#################################################
#Function to interface between the user and XPEB
#Exported Function: run.xpeb()
#################################################

#################################################
#Description
#################################################
#Loads and checks the target and base-GWAS files, containing GWAS statistics of type p-value for the target and base.
#Converts target and base GWAS p-values into chi-squares.
#Applies genomic control correction if user-specified.
#Runs 2 rounds of XPEB including trimming for LD.
#Returns: XPEB locfdr and estimated overlap.

#################################################
#Usage example
#################################################
# results <- run.xpeb(path.target="inputData/target.10K.assoc.linearWITHNA",path.base="inputData/base.10K.metaWITHNA",gc.target=T,gc.base=T,n.target=10000,n.base=10000,n.iter=1e5)
# results$overlap

#################################################
#Input options
#################################################
#path.target and path.base: string indicating the path to the GWAS result files in the target and base populations.
			#Requirements for these files:			
			#1 file for base and 1 file for target
			#Target file contains 4 mandatory columns with header: P SNP BP CHR
			#Target file contains 2 mandatory columns with header: P SNP
			#1 marker data per row, and no duplicates marker names. 
				#Note: If using PLINK and there are covariates, suppress the covariates in the association output files using the PLINK option: --hide-covar
			#white space delimited (1 or more space, tab, newlines or carriage return)
#gc.target and gc.base: T or F to indicate whether a genomic control correction should be applied (T) or not (F) 
			#to the target and base-GWAS satistics. (defaults to T)
#n.target and n.base: median sample size in target and base.
#n.iter: number of iterations for the MCMC.

#################################################
#Ouput
#################################################
#A list containing 2 elements:
#1. data.frame: dimmensions: number of markers in target x 4 columns
				#columns are: SNP,CHR,BP,LOCFDR
  
#2.the estimated overlap


#################################################
#Checks performed by this function on the input files and parameters
#################################################
#Check that the {SNP P CHR BP} and {SNP P} columns are present and not duplicated in the target and base-GWAS files, respectively.
#Check that no marker name is NA.
#Check that BP are numeric and not NA.
#Check that CHR are strings or numeric and not NA. 
#Check if some SNP IDs are duplicated : terminate if there are duplicates
	#Another option is to take the first instance and remove the others and issue a warning how many are removed- stop if more than 10% are removed
#Check that there is a min number of markers in the target and base >=minMarker (1e5)
#Check P values are all numeric, with P between 0 and 1 or NA.
#Check that, out of the markers with data in base and target, 
	#the number of markers in target that are genome-wide significant (p<p.Threshold 5e-8) is >=minSign (3)
	#the number of markers in base that are genome-wide significant (p<p.Threshold 5e-8) is >=minSign (3)
#Check that n.base and n.target are numeric and >=minSample (100)
#Check that n.iter is numeric and >=minNiter (1e5)
#Check that gc.base and gc.target are T or F


#################################################
#run.XPEB is the only exported function
#################################################
run.xpeb <- function(path.target,path.base,n.target,n.base,gc.target=T,gc.base=T,n.iter=1e6){
	##################	
	#Print start time
	##################
	print(paste("Starting XPEB on ",Sys.time(),sep=""))
		
	##################	
	#process the user input data
	##################
	processed.Input <- xpeb.importData(path.target=path.target,path.base=path.base,gc.target=gc.target,gc.base=gc.base, n.target=n.target,n.base=n.base,n.iter=n.iter)
	#processed.Input is a list with:
	#data with no NA: processed.Input$ceb.input.noNA  chi.stats and keep.markers
	#full user data, chi squares: processed.Input$ceb.stats.full
	#full user data, marker info : processed.Input$marker.info  
	
	##################
	#Run XPEB with pruning for LD
	##################
	print("----------------------------------------------------------------------")
	print("Running XPEB")
	print("----------------------------------------------------------------------")
	
	results <- run.xpeb.internal(dat.noNA=processed.Input$ceb.input.noNA,stats.full=processed.Input$ceb.stats.full, n.base=n.base, n.target=n.target,n.iter=n.iter,marker.info.all=processed.Input$marker.info)

	##################	
	#Print end time
	##################
	print(paste("XPEB finished on ",Sys.time(),sep=""))
	print(sprintf("Estimated overlap is:  %.3f.",results$overlap))
	print(sprintf( "%i markers have an XPEB locfdr<0.05, out of %i markers with P-value in the target.",sum(results$locfdr$LOCFDR<0.05,na.rm=T),sum(!is.na(results$locfdr$LOCFDR))))

	##################
	#Return the results
	##################
	return(list(locfdr=results$locfdr,overlap=results$overlap))

}


#################################################
#run.xpeb.internal() runs XPEB with LD pruning
#################################################
#one global variable is defined:
#max overlap
p11.constraint=.9

run.xpeb.internal <- function(dat.noNA,stats.full, n.base, n.target, n.iter,marker.info.all){
	##################
	#Specify parameters
	##################
	#for the ld trimming:
	ld.trim.prob=0.5;
	ld.trim.bp=2e5;
	#grid length
	grid.length=500
	
	##################
	#run EB round 1
	##################
	print("Running EB round 1...")
	ceb.obj.shared0=analyzeArun5(dat=dat.noNA$chi.stats, n.base=n.base, n.target=n.target, n.iter=n.iter,grid.length=grid.length)

	##################
	# pruning for LD the "suggestive loci"
	##################
	post.base0=ceb.obj.shared0$post.alt.base
	post.target0=ceb.obj.shared0$post.alt.target
	post.xpeb0=ceb.obj.shared0$eb$post.prob.alt.target
	markermax=apply(cbind(post.base0, post.target0, post.xpeb0), 1, max) #9997
	marker.info=marker.info.all[dat.noNA$keep.markers,] #9997
	dummy=(markermax>ld.trim.prob) #9997    - 679 T rest F

	#identifying windows (200kb separation)
	cluster0=makeCluster(dummy=(markermax>ld.trim.prob), marker.info=marker.info, ld.trim.bp=ld.trim.bp, flank=1e4)
	window.mat=cluster0[[1]]; #610 6
	rm.indx=unlist(cluster0[[2]]) #list with 610 elements ->735 indeces to rm

	#making the prunned dataset
	newdat=dat.noNA$chi.stats[nrow(window.mat),];
	for (j in 1:nrow(window.mat)) {
	  zz=which.max(post.xpeb0[window.mat[j,1]: window.mat[j,2]]);
	  newdat[j,]=as.numeric(dat.noNA$chi.stats[window.mat[j,1]-1+zz,])
	}
	#the new dataset has 1 SNP in each of 610 windows  + the windows with no sign finding
	#this is the chi square data trimmed for ld in the suggestive windows
	pruned.chi.dat=rbind(dat.noNA$chi.stats[-rm.indx,], newdat); #this is  9872    2
	
	##################
	#run EB round 2 on the ld trimmed dataset
	##################
	print("Running EB round 2...")
	ceb.obj.shared1=analyzeArun5(dat=pruned.chi.dat, n.base=n.base, n.target=n.target, n.iter=n.iter,grid.length=grid.length)
		
	##################
	#fix up the mixture parameters
	##################
	ceb.obj.shared.fixed=ceb.obj.shared1
	ceb.obj.shared.fixed$eb[[1]]=c(ceb.obj.shared1$eb[[1]][1], max(0, min(-ceb.obj.shared1$eb[[1]][1]+sum(ceb.obj.shared0$eb[[1]]), -ceb.obj.shared1$eb[[1]][1]+ilogistic(p11.constraint))))

	
	##################
	#use the new parameters to predict locfdr for all markers
	##################
	posterior.full=eb.predict.postprob(dat=stats.full, res=ceb.obj.shared.fixed)
	betas.final=ceb.obj.shared.fixed$eb[[1]]

	##################		
	#calculate the overlap:
	##################
	overlap <- exp(betas.final[1]+betas.final[2])/(1+exp(betas.final[1]+betas.final[2]))

	##################
	#assemble the results:
	##################
		
	locfdr <- cbind(marker.info.all,1-posterior.full)
	names(locfdr)<- c("SNP","CHR","BP","LOCFDR")
	
	return(list(locfdr=locfdr,overlap=overlap))
}





#################################################
#xpeb.importData() import users p-value and outputs XPEB input
#################################################
xpeb.importData <- function(path.target,path.base,gc.target=T,gc.base=T,n.target,n.base,n.iter){
	
	###Set some constants as bounds for the input checks
	minSample <- 100     #requires that n.target and n.base be >100
	minMarker <- 1e5     # require at least this many markers in the target and base files
	minSign <- 3         # require at least this many GW sign markers in base and target
	p.Threshold <- 5e-8  # GW significance threshold in term of p value	
	minNiter <- 1e5

	###Columns to retrieve from base and target
	col.target <- c("SNP","P","CHR","BP")
	col.base <- c("SNP","P")		

	
	###First check that the function arguments are in the correct format
	check.arguments.importData(gc.target=gc.target,gc.base=gc.base,n.target=n.target,n.base=n.base,minSample=minSample,n.iter=n.iter,minNiter=minNiter)
	
	###Print to the log the options chosen by the user.
	print("----------------------------------------------------------------------")
	print("User input parameters for xpeb.importData")	
	print("----------------------------------------------------------------------")
	print(paste("Target-GWAS file is: ",path.target,sep=""))
	print(paste("Sample size for the target-GWAS is: ",n.target,".",sep=""))	
	print(paste("Sample size for the base-GWAS is: ",n.base,".",sep=""))	
	print(paste("The number of individuals in the base-GWAS is: ",n.iter,".",sep=""))
	if(gc.target){
		print("Applying genomic control correction to the target-GWAS.")
	} else(print("No genomic control correction to the target-GWAS."))
	print(paste("Base-GWAS file is: ",path.base,sep=""))
	if(gc.base){
		print("Applying genomic control correction to the base-GWAS.")
	} else(print("No genomic control correction to the base-GWAS."))
	print(paste("The number of iterations for the MCMC is: ",n.iter,".",sep=""))
	print("----------------------------------------------------------------------")
	print("Preparing target and base-GWAS data for XPEB")
	print("----------------------------------------------------------------------")

	
	###Load and check the target-GWAS file
	target <- read.table(path.target,as.is=T,header=T)
	#make an index so we can recover the order of the user markers in target
	target$index.user <- 1:dim(target)[1]
	#check that the target file has the right columns, and that the stats look like proper stats.
	check.colNames(col=col.target,header=names(target),gwas="target")
	check.markers(markers=target[,col.base[1]],gwas="target",minMarker=minMarker)
	check.stat(stat=target[,col.target[2]],gwas="target")
	check.chr(chr=target[,col.target[3]],gwas="target")
	check.bp(bp=target[,col.target[4]],gwas="target")
	
	#Feedback to user that the target was successfully read
	print(paste("The file for the target GWAS was read. It contains P-values for ",dim(target)[1]," markers.",sep=""))

	###Load and check the base-GWAS file
	base <- read.table(path.base,as.is=T,header=T)
	check.colNames(col=col.base,header=names(base),gwas="base")
	check.markers(markers=base[,col.base[1]],gwas="base",minMarker=minMarker)
	check.stat(stat=base[,col.base[2]],gwas="base")
	#Feedback to user that the target was successfully read
	print(paste("The file for the base GWAS was read. It contains P-values for ",dim(base)[1]," markers.",sep=""))
	

	###Merge the target and base file, keeping all target markers (this is useful since we can then sort by the initial index in target, which is coming from a user file, and keep this order throughout)
	print("Merging target and base-GWAS files. Keeping all markers in target-GWAS.")
	data <- merge(target[,c(col.target,"index.user")],base[,col.base],by.x=col.target[1],by.y=col.base[1],all.x=T)
	names(data) <- c("SNP","P.target","CHR","BP","index.user","P.base")
	
	#re-sort by the order in target and reorganize the columns
	data <- data[order(data$index.user),]
	data <- data[,c("P.target","P.base","SNP","index.user","CHR","BP")]

	#Give the user feedback on the merge
	print(sprintf("Out of %i markers,",dim(data)[1]))
	print(sprintf("%i have non missing values (NA) in the base and target-GWAS,",sum(!is.na(data[,"P.target"])&!is.na(data[,"P.base"]))))
	print(sprintf("%i have missing values in the target-GWAS only,",sum(is.na(data[,"P.target"])&!is.na(data[,"P.base"]))))
	print(sprintf("%i have missing values in the base-GWAS only,",sum(!is.na(data[,"P.target"])&is.na(data[,"P.base"]))))
	print(sprintf("and %i have missing values in both the target and base-GWAS.",sum(is.na(data[,"P.target"])&is.na(data[,"P.base"]))))

	#Check that at least minMarker have data in the base and the target:
	check.markerNbr(dat=data,minMarker=minMarker)
	
	###Prepare the data as chi-square for XPEB:
	ceb.input.noNA <- prepData(dat=data,gc.target=gc.target,gc.base=gc.base,na.action=-1)
	ceb.input.full <- prepData(dat=data,gc.target=gc.target,gc.base=gc.base,na.action=0)
		
	###Print the lambda GC:
	print(paste("For the ",sum(!is.na(data[,"P.target"])&!is.na(data[,"P.base"]))," markers with statistics available in both the target and base-GWAS,",sep="")	)
	print(sprintf("Lambda GC in the target-GWAS is %.3f.",ceb.input.noNA$lambda[1]))
	print(sprintf("Lambda GC in the base-GWAS is %.3f.",ceb.input.noNA$lambda[2]))
	
	###Check that there is enough signal in the target and in the base
	check.signal(chi.vector=ceb.input.noNA$chi.stats[,1],gwas="target",minSign=minSign,p.Threshold=p.Threshold)
	check.signal(chi.vector=ceb.input.noNA$chi.stats[,2],gwas="base",minSign=minSign,p.Threshold=p.Threshold)

	#keep what we need
	ceb.input.noNA <- ceb.input.noNA[c("chi.stats","keep.markers")]

	###Output the data to be used to run XPEB:
	return(list(ceb.input.noNA=ceb.input.noNA,ceb.stats.full=ceb.input.full$chi.stats,marker.info=data[,c(3,5,6)]))
	
}


#############################################
#This function checks the gc correction switch option is a boolean T or F
##############################################
check.arguments.importData <- function(gc.target,gc.base,n.target,n.base,minSample,n.iter,minNiter){
	#gc can only be T or F
	if(!is.logical(gc.target)|is.na(gc.target)){stop("GC correction (gc.target) must be specified as a boolean: T or F.")}
	if(!is.logical(gc.base)|is.na(gc.base)){stop("GC correction (gc.base) must be specified as a boolean: T or F.")}
	if(!is.numeric(n.target)|is.na(n.target)){stop("Sample size in target-GWAS must be specified as a number.")}
	if(!is.numeric(n.base)|is.na(n.base)){stop("Sample size in base-GWAS must be specified as a number.")}
	if(n.target<minSample){stop(paste("Sample size in target-GWAS is not sufficient. It needs to be greater than or equal to",minSample,".",sep=""))}
	if(n.base<minSample){stop(paste("Sample size in base-GWAS is not sufficient. It needs to be greater than or equal to ",minSample,".",sep=""))}
	if(!is.numeric(n.iter)|is.na(n.iter)){stop("The option n.iter must be specified as a number.")}
	if(n.iter<minNiter){stop(paste("The number of iterations (n.iter) must be greater than or equal to ",minNiter,".",sep=""))}

}

#############################################
#This function checks that the column names exist and are not duplicated.
#############################################
check.colNames <- function(col,header,gwas){
	#check if the column names are present:
	if (sum(is.element(col,header))!=length(col)){
		#find which col is missing:
		tmp <- col[!is.element(col,header)] 
		#give error message
		stop(paste("Column(s) \"",capture.output(cat(tmp)),"\" cannot be found in the ",gwas,"-GWAS.",sep=""),call.=F)
	}
	#make sure they are no duplicated column names:
	if(length(header[is.element(header,col)])!=length(col)){
		tmp <- header[is.element(header,col)]
		dup <-  names(table(tmp)[table(tmp)>1])
		stop(paste("There are multiple instances of column(s): \"",capture.output(cat(dup)),"\" in the ",gwas,"-GWAS.",sep=""),call.=F)
	}
}

#############################################
#This function checks that there are no NA markers, no duplicated marker names, and the min number of markers .
#############################################
check.markers <- function(markers,gwas,minMarker){
	if(sum(is.na(markers))>0){
		stop(paste("There are marker(s) with names \"NA\" in the ",gwas,"-GWAS file.\nRemove these marker(s) before running XPEB.",sep=""),call.=F)
	}
	if(length(unique(markers))!=length(markers)){
		stop(paste("There are duplicated marker names in the ",gwas,"-GWAS file.\nRemove duplicated markers before running XPEB.",sep=""),call.=F)
	}	
	if(length(markers)<minMarker){
		stop(paste("There are less than ",minMarker," markers in the ",gwas,"-GWAS file.\nThis is too small of a dataset to run XPEB.",sep=""),call.=F)
	}
}

#############################################
#This function checks that there at least minMarker with data in both the target and base-GWAS.
#############################################
check.markerNbr <- function(dat,minMarker){
	markersWithData <- sum(!is.na(dat[,1])&!is.na(dat[,2]))
	if(markersWithData<minMarker){
		stop(paste("There are less than ",minMarker," markers with non missing data in both the target and base-GWAS.\nThis is too small of a dataset to run XPEB.",sep=""),call.=F)
	}
}


#############################################
#This function checks that the statistics are numeric and between 0 and 1 for the p-values.
##############################################
check.stat <- function(stat,gwas){
	#whether it is  numeric
	if(class(stat)!="numeric"){
		stop(paste("Statistics in ",gwas,"-GWAS file are not all numeric.",sep=""),call.=F)
	}
	#for a p-value, check that it is between 0 and 1
	if(min(stat,na.rm=T)<0|max(stat,na.rm=T)>1){
			stop(paste("Some P-values in ",gwas,"-GWAS file are less than 0 or greater than 1.",sep=""),call.=F)
	}
	
}


#############################################
#This function checks that the chr information is numeric or strings and that none are NA
##############################################
check.chr <- function(chr,gwas){
	if(!(class(chr)=="numeric"|class(chr)=="character"|class(chr)=="integer")){
		stop(paste("The chromosome column should be of type \"number\" or \"character\" but instead it is a ",class(chr)," in the ",gwas,"-GWAS file.",sep=""))
	}
	if(sum(is.na(chr))>0){
		stop(paste("There are some NA values for the chromosome column in the ",gwas,"-GWAS file.",sep=""))
	}
}


#############################################
#This function checks that the position information are numeric and that none are NA
##############################################
check.bp <- function(bp,gwas){
	if(!(class(bp)=="numeric"|class(bp)=="integer")){
		stop(paste("Position in ",gwas,"-GWAS file are not all numeric.",sep=""),call.=F)
	}
	if(sum(is.na(bp))>0){
		stop(paste("There are some NA values for the position column in the ",gwas,"-GWAS file.",sep=""))
	}
}


#############################################
#This function checks that there is enough signal in the base and target: 4 or more are p<5e-8 <=> chi2.df1>29.71679
##############################################
check.signal <- function(chi.vector,gwas,minSign,p.Threshold){
	chiThreshold <- qchisq(p.Threshold, df=1, lower.tail=F) #29.71679 for p 5e-8
	if(sum(chi.vector>=chiThreshold,na.rm=T)<minSign+1){
		stop(paste("Less than ",minSign+1," markers are genome-wide significant (p<",p.Threshold,") in ",gwas,"-GWAS file.\nThere is not enough information to run XPEB.",sep=""),call.=F)		
	}
}

#############################################
#Functions that transforms a p value into the XPEB input chi-square
#############################################
p2chi=function(pval, df=1, thresh=1e-15, lower=F) {
	pval=pmin(pval, 1-thresh) #this replaces p values of 1 by 1-1e-15 so that the min chisquare is not 0 but  1.568286e-30
	chi.val=qchisq(pval, df=df, lower.tail=lower)#p=NA is chi=NA; p=0 is chi=Inf,p=1e-300 chi is 1373, p=1 is now chi= 1.568286e-30	
	#as far as I can tell the next 3 lines could be removed
	nas=which(!is.finite(chi.val))
	if (length(nas))
		chi.val[nas]=qchisq(1-pval[nas], df=df, lower.tail=1-lower)
	##end of part that could be removed	
	return(chi.val);
}


#############################################
#Function that implements the GC correction of chi-squareand calculates lambda GC
#############################################
gc.scale.chisq=function(chi,gc.scale,df=1) {
  lambda=pmax(median(chi,na.rm=T)/qchisq(.5,df=df,lower.tail=FALSE) ,1)
  if(gc.scale){	
	 chi1=chi / lambda
	 return(list(chi1=chi1, lambda=lambda))
  } else {
  	 return(list(chi1=chi, lambda=lambda))
  }
}

#############################################
#Function that transforms input data (P values) into chi-square statistics for the base and target-GWAS.
#It outputs a list with chi.stats, lambda and keep.markers.
#############################################
prepData=function(dat,type.target,gc.target,type.base,gc.base,na.action) {
	
	###Get the index of markers to keep: zz
	if (na.action==-1) {
		zz=which( (!is.na(dat[,1])) &  (!is.na(dat[,2])))  #keep markers with info in both the target and the base 
	} else if (na.action==1) {
		zz=which( (!is.na(dat[,1])) )  #keep markers with info in the target
	} else if (na.action==2) {
		 zz=which( (!is.na(dat[,2])) )  #keep markers with info in the base
	} 
	else if(na.action==0){
		 zz=1:nrow(dat); ## keep all
	 }

	###subset the data on the markers to keep        
	dat=dat[zz,]

	###Calculate the chi-square,
	#for the target
	chisq.target=p2chi(dat[,1])
	#for the base
	chisq.base=p2chi(dat[,2])

	#Do GC-correction if the user wants to otherwise still return the lambdaGC for information.
	chisq.target= gc.scale.chisq(chisq.target,gc.scale=gc.target)    
	chisq.base= gc.scale.chisq(chisq.base,gc.scale=gc.base)    

	return(list(chi.stats=data.frame(chisq.target=chisq.target[[1]], chisq.base=chisq.base[[1]]), keep.markers=zz, lambda=c(chisq.target[[2]],chisq.base[[2]])));
}












