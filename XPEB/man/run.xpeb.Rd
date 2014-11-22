\name{run.xpeb}
\alias{run.xpeb}
\title{Calculates the XPEB locfdr and overlap for a target and base-GWAS}
\description{
  This function runs XPEB. It takes as input the target and base-GWAS files, containing the GWAS p-values for the target and base.  It applies genomic control correction per user specification. Then it runs XPEB including trimming for LD based on base pair position. The function returns the XPEB locfdr and estimated overlap.
}
\usage{
run.xpeb(path.target,path.base,n.target,n.base,gc.target=T,gc.base=T,n.iter=1e6)
}

\arguments{
  \item{path.target}{string indicating the path to the target-GWAS file.}
  
  \item{path.base}{string indicating the path to the base-GWAS file.} 

  \item{n.base}{numeric, number of individuals in the base-GWAS}

  \item{n.target}{numeric, number of individuals in the target-GWAS.}

  \item{gc.target}{logical. Should genomic controls be applied to the target-GWAS? True by default.}

  \item{gc.base}{logical. Should genomic controls be applied to the base-GWAS? True by default.}

  \item{n.iter}{numeric: number of iterations for MCMC.}
}

\value{
  \item{locfdr}{a data frame giving for each marker, the locfdr.}
  
  \item{overlap}{numeric: the overlap in genetic architecture between the two populations.}
}

\references{Marc A Coram, Sophie I Candille, Qing Duan, Katie Chan, Yun Li, Charles Kooperberg, Alex P Reiner, Hua Tang. An Ethnicity-Aware Approach for Mapping Complex Traits Using Multi-ethnic Data. American Journal of Human Genetics (submitted)}

\examples{
	#Unzip the example files from the package and retrieve their path
	#target-GWAS
	path.target <- system.file("extdata", "target.gwas.txt.zip", package="XPEB")
	unzip(path.target,exdir="TMPinput")
	path.target <- "TMPinput/target.gwas.txt"
	#base-GWAS
	path.base <- system.file("extdata", "base.gwas.txt.zip", package="XPEB")
	unzip(path.base,exdir="TMPinput")
	path.base <- "TMPinput/base.gwas.txt"
	
	#Run XPEB on the example files
	res <- run.xpeb(path.target=path.target,path.base=path.base,n.target=1e4,n.base=1e5)
	
	#print the overlap
	print("The estimated overlap in genetic architecture is:", res$overlap, sep="")
	
	#Save the locfdr calculation results to a text file
	write.table(res$locfdr, file="locfdrResults.txt", sep="\t", quote=FALSE, row.names=FALSE) 
}
 
