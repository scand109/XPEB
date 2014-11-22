#############################
#Instructions to install and run the XPEB package
#############################

#to install package, at the command prompt run:
#>R CMD INSTALL XPEB_0.1.tar.gz


library(XPEB)

#Retrieve input example from the package
path.target <- system.file("extdata", "target.gwas.txt.zip", package="XPEB")
unzip(path.target,exdir="TMPinput")
path.target <- "TMPinput/target.gwas.txt"

path.base <- system.file("extdata", "base.gwas.txt.zip", package="XPEB")
unzip(path.base,exdir="TMPinput")
path.base <- "TMPinput/base.gwas.txt"

#run XPEB
res <- run.xpeb(path.target=path.target,path.base=path.base,n.target=10000,n.base=100000,gc.target=T,gc.base=T,n.iter=1e5)

#locfdr results: 
#res$locfdr$LOCFDR
#overlap: 
#res$overlap





