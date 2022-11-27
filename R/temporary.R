
printf <- function(...) invisible(print(sprintf(...)))


resetInferN0 <- function(){
	detach("package:InferN0", unload=TRUE)
	library(InferN0)
}

testA <-function(){
	library(InferN0)
scptr<- infern0LoadFile(path="/lustre/scratch117/cellgen/team218/lh20/inferN0output/frozenhigh.ifn.scp")
	genelist = c("BCL11B","CUX1","CUX2","ETV1","FEZF2","NR4A2","POU3F2","RELN","RORB","SATB2","SOX5","CNTNAP2","CALB1","CALB2","CCK","GAD1","GAD2","NPY","SST","TH","ERBB4","MEF2C","NPAS3","GRIA1","NLGN1","DLX6-AS1","VIM","APOE","APP","PSEN1","STMN2","ROBO1","TOP2A","HIST1H4C","H3F3B","MALAT1","RPLP0","RPL41","RPL12","MT-ND2")
	
return(infern0ComputeCovariance(scptr,genelist,nb.threads=4))}

testB <-function(){
	library(InferN0)
scptr <- infern0LoadFile(path="/lustre/scratch117/cellgen/team218/lh20/inferN0output/prefrontalcortex.ifn.scp")
	genelist = c("BCL11B","CUX1","CUX2","ETV1","FEZF2","NR4A2","POU3F2","RELN","RORB","SATB2","SOX5","CNTNAP2","CALB1","CALB2","CCK","GAD1","GAD2","NPY","SST","TH","ERBB4","MEF2C","NPAS3","GRIA1","NLGN1","DLX6-AS1","VIM","APOE","APP","PSEN1","STMN2","ROBO1","TOP2A","HIST1H4C","H3F3B","MALAT1","RPLP0","RPL41","RPL12","MT-ND2")
	output <- infern0ComputeCovariance(scptr,genelist,nb.threads=1)
return(output)}



testC <- function(){
	library(InferN0)
scptr <- infern0LoadFile("inferN0output/prefrontalcortex.ifn.scp"); 
genelist = c("BCL11B","CUX1","CUX2","ETV1","FEZF2","NR4A2","POU3F2","RELN","RORB","SATB2","SOX5","CNTNAP2","CALB1","CALB2","CCK","GAD1","GAD2","NPY","SST","TH","ERBB4","MEF2C","NPAS3","GRIA1","NLGN1","DLX6-AS1","VIM","APOE","APP","PSEN1","STMN2","ROBO1","TOP2A","HIST1H4C","H3F3B","MALAT1","RPLP0","RPL41","RPL12","MT-ND2");
#outlist <- infern0ComputeCovariance(scptr, genelist, method="Partial")
outlist <- infern0ComputeCovariance(scptr, genelist, method="Partial")
outlist <- infern0ComputeCovariance(scptr, genelist, method="EM")
#outlist <- infern0ComputeCovariance(scptr, genelist, method="EM", nb.crosseval.partition = 8)

outlist2 <- infern0IdentifyNetwork(scptr, genelist, method="EMM")


runrun <- Infern0IdentifyNetwork(scptr, genelist, method="Partial")

}


testD <- function(){
	library(InferN0)
scptr <- infern0LoadFile("inferN0output/prefrontalcortex.ifn.scp"); 
genelist = c("BCL11B","CUX1","CUX2","ETV1","FEZF2","NR4A2","POU3F2","RELN","RORB","SATB2","SOX5","CNTNAP2","CALB1","CALB2","CCK","GAD1","GAD2","NPY","SST","TH","ERBB4");
#outlist <- infern0ComputeCovariance(scptr, genelist, method="Partial")
outlist <- infern0ComputeCovariance(scptr, genelist, method="Partial")
outlist <- infern0ComputeCovariance(scptr, genelist, method="EM")

}

errornorm <- function(X,Y){
	nb_edge = sum(Y) / 2;
	nb_nonedge = sum((-diag(dim(Y)[1]) +1) - Y)/2
	count =matrix(0,2);
	d = dim(X);
	for(i in 2:d[1])
	for(j in 1:(i-1))
	if (X[i,j] != Y[i,j]) {
		if (X[i,j] < Y[i,j]) count[1]= count[1]+1
		else count[2]= count[2]+1
	}
	return(c(count[1]/ nb_edge, (nb_nonedge -count[2])/ nb_nonedge))
}

IdentifyNetworkQUIC <- function(findcov.out, M=-1, lambdas=(1:50)/200, excludezeros=T) {
    require('QUIC')
    S <- findcov.out$Covar
    n <- dim(S)[2]
    Z <- c()
    ms <- c()
    Z.list <- list()
    for (j in 1:length(lambdas)) {
        quicout <- QUIC(S, lambdas[j], msg=0)$X
				Z <- 1*(abs(quicout)>1e-6)-diag(n)
        Z.list[[length(Z.list)+1]] <- Z
    }
    return( Z.list )
}



testE <- function(){
	library(InferN0)
	Z <- GetGoldStandardPluripotencyNetwork();
	X <- infern0GenerateSyntheticNetworkData(Z, p.min=.5, p.max=.95, 1000);
	out <- infern0ComputeCovariance(X, method="CorrectedPartial",nb.crosseval.partition=8)
	outnold <- infern0IdentifyNetworkOlderVersion(out,  evaltarget=Z,do.manually.select.nb.edges=T,check.nb=37)
	dev.new()
	infern0PlotNetwork(outnold)
}


testF <- function(Z = GetGoldStandardPluripotencyNetwork(), rep = 250, n = 1000){

	#out2 <- infern0ComputeCovariance(Matrix(X,sparse=T), method="")
	#out2 <- infern0ComputeCovariance(Matrix(X,sparse=T), method="EMM")
	#outn <- infern0IdentifyNetwork(X, method=out,  evaltarget=Z)
	#outeval <- infern0FindConstrainedCovariance(outn$PrecisionMatrix !=0, out)
	#outn$PrecisionMatrix
	#outeval$PrecisionMatrix
	#outn <- infern0IdentifyNetwork(X, method=out,  evaltarget=Z,do.manually.select.nb.edges=T,check.nb=37)
	#outeval <- infern0FindConstrainedCovariance(outn$PrecisionMatrix !=0, out)

	
	maxcov_stat <- matrix(0,2,50)
	quic_stat <- matrix(0,2,50)
	quic_nb <-0
	

	for( i in 1:rep){
		X <- infern0GenerateSyntheticNetworkData(Z, n);
		out <- infern0ComputeCovariance(X, method="EMM",nb.crosseval.partition=8)

    mat_arr <- try(IdentifyNetworkQUIC(out,excludezeros=T))
		if (class(mat_arr) != "try-error"){
			quic_nb %+=% 1
			for(j in 1:50){
			  quic_stat[,j] %+=% t(errornorm(Tnet,mat_arr[[j]]))
			}
		}

	}

	quic_stat %/=% quic_stat

	plot(quic_stat[,1], quic_stat[,2])

}

testG <- function(){
	library(InferN0)
	source("~/work/readFromNew.R")
	scptr <- infern0LoadFile("/lustre/scratch117/cellgen/team218/lh20/inferN0output/frozenhigh.ifn.scp");subcell <- readRDS("tmp.rds");
	out <- infern0ComputeCovariance(scptr,gene.list= c("BCL11B","CUX1","CUX2", "ETV1","FEZF2", "NR4A2","POU32F2","RELN","RORB","SATB2","SOX5","GRIA1","CALB2","ROBO1","ROBO2","GAD2","DLX5","CNTNAP2","SIX3","MEIS2","MEF2C","NPAS3","GRIK2","CCND2","SOX4"), cell.list=subcell,nb.crosseval.partition=8, method="EM")
	outnold <- infern0IdentifyNetworkOlderVersion(out,check.nb=37)
}

ExtractNetworkfromCoordinates <- function(table, nbedges=0, n=0){
	if (n==0) n = max(table[,1:2])+1
	fout <- diag(n)
	if (table[1,1] == table[1,2]) table <- table[2:(dim(table)[1]),]
	if (nbedges==0) nbedges = dim(table)[1]
	for(i in 1:nbedges){
		fout[table[i,2],table[i,1]] <- 1
	  fout[table[i,1],table[i,2]] <- 1
	}
return(fout)}




#Return a network of 15 genes from the core pluripotency network. If the argument ns contains gene names, only those genes present in ns will be included.
GetGoldStandardPluripotencyNetwork <- function(ns=c()) {
    if (length(ns)==0) {
        ns <- c("Esrrb", "Jarid2", "Klf4", "Tbx3", "Nanog", "Pou5f1", "Sox2", "Tcf3", "Stat3", "Zfp281", "Myc", "Nr0b1", "Zfx", "Sall4", "Zfp42")
    }
    x <- matrix(rep(0, length(ns)^2), length(ns))
    rownames(x) <- ns
    colnames(x) <- ns
    x[which(ns=="Esrrb"), which(ns=="Jarid2")] <- 1
    x[which(ns=="Esrrb"), which(ns=="Klf4")] <- 1
    x[which(ns=="Esrrb"), which(ns=="Tbx3")] <- 1
    x[which(ns=="Esrrb"), which(ns=="Nanog")] <- 1
    x[which(ns=="Esrrb"), which(ns=="Pou5f1")] <- 1
    x[which(ns=="Jarid2"), which(ns=="Sox2")] <- 1
    x[which(ns=="Jarid2"), which(ns=="Tcf3")] <- 1
    x[which(ns=="Jarid2"), which(ns=="Sox2")] <- 1
    x[which(ns=="Stat3"), which(ns=="Sox2")] <- 1
    x[which(ns=="Stat3"), which(ns=="Zfp281")] <- 1
    x[which(ns=="Stat3"), which(ns=="Myc")] <- 1
    x[which(ns=="Klf4"), which(ns=="Myc")] <- 1
    x[which(ns=="Klf4"), which(ns=="Nanog")] <- 1
    x[which(ns=="Klf4"), which(ns=="Pou5f1")] <- 1
    x[which(ns=="Klf4"), which(ns=="Tbx3")] <- 1
    x[which(ns=="Klf4"), which(ns=="Zfp42")] <- 1
    x[which(ns=="Tbx3"), which(ns=="Sox2")] <- 1
    x[which(ns=="Tbx3"), which(ns=="Pou5f1")] <- 1
    x[which(ns=="Tbx3"), which(ns=="Nr0b1")] <- 1
    x[which(ns=="Tbx3"), which(ns=="Tcf3")] <- 1
    x[which(ns=="Tbx3"), which(ns=="Zfx")] <- 1
    x[which(ns=="Zfp42"), which(ns=="Sox2")] <- 1
    x[which(ns=="Zfp42"), which(ns=="Pou5f1")] <- 1
    x[which(ns=="Zfp42"), which(ns=="Nr0b1")] <- 1
    x[which(ns=="Zfp42"), which(ns=="Zfx")] <- 1
    x[which(ns=="Tcf3"), which(ns=="Pou5f1")] <- 1
    x[which(ns=="Tcf3"), which(ns=="Sox2")] <- 1
    x[which(ns=="Tcf3"), which(ns=="Nanog")] <- 1
    x[which(ns=="Tcf3"), which(ns=="Sall4")] <- 1
    x[which(ns=="Tcf3"), which(ns=="Nr0b1")] <- 1
    x[which(ns=="Nr0b1"), which(ns=="Sox2")] <- 1
    x[which(ns=="Nr0b1"), which(ns=="Nanog")] <- 1
    x[which(ns=="Nr0b1"), which(ns=="Sall4")] <- 1
    x[which(ns=="Nr0b1"), which(ns=="Pou5f1")] <- 1
    x[which(ns=="Sall4"), which(ns=="Pou5f1")] <- 1
    x[which(ns=="Nanog"), which(ns=="Pou5f1")] <- 1
    x[which(ns=="Sox2"), which(ns=="Pou5f1")] <- 1
    x[which(ns=="Sox2"), which(ns=="Nanog")] <- 1
    x <- x + t(x)
    return( x )
}

#Similar as before but a sparser hematopoiesis network
GetGoldStandardHematopoieticNetwork <- function(ns=c()) {
    if (length(ns)==0) {
        ns <- c("Bmp4", "Smad1", "Runx1", "Smad6", "Erg", "Eng", "Lmo2", "Lyl1", "Fli1", "Scl", "Runx3", "Hex", "Gata2", "Pu.1", "Elf1")
    }
    x <- matrix(rep(0, length(ns)^2), length(ns))
    rownames(x) <- ns
    colnames(x) <- ns
    x[which(ns=="Elf1"), which(ns=="Lmo2")] <- 1
    x[which(ns=="Elf1"), which(ns=="Scl")] <- 1
    x[which(ns=="Pu.1"), which(ns=="Lyl1")] <- 1
    x[which(ns=="Pu.1"), which(ns=="Eng")] <- 1
    x[which(ns=="Pu.1"), which(ns=="Lmo2")] <- 1
    x[which(ns=="Gata2"), which(ns=="Hex")] <- 1
    x[which(ns=="Gata2"), which(ns=="Scl")] <- 1
    x[which(ns=="Gata2"), which(ns=="Lyl1")] <- 1
    x[which(ns=="Gata2"), which(ns=="Eng")] <- 1
    x[which(ns=="Gata2"), which(ns=="Runx1")] <- 1
    x[which(ns=="Gata2"), which(ns=="Smad6")] <- 1
    x[which(ns=="Gata2"), which(ns=="Fli1")] <- 1
    x[which(ns=="Scl"), which(ns=="Hex")] <- 1
    x[which(ns=="Scl"), which(ns=="Runx3")] <- 1
    x[which(ns=="Scl"), which(ns=="Lyl1")] <- 1
    x[which(ns=="Scl"), which(ns=="Lmo2")] <- 1
    x[which(ns=="Fli1"), which(ns=="Hex")] <- 1
    x[which(ns=="Fli1"), which(ns=="Lyl1")] <- 1
    x[which(ns=="Fli1"), which(ns=="Lmo2")] <- 1
    x[which(ns=="Fli1"), which(ns=="Eng")] <- 1
    x[which(ns=="Erg"), which(ns=="Scl")] <- 1
    x[which(ns=="Erg"), which(ns=="Lyl1")] <- 1
    x[which(ns=="Erg"), which(ns=="Eng")] <- 1
    x[which(ns=="Bmp4"), which(ns=="Smad6")] <- 1
    x[which(ns=="Bmp4"), which(ns=="Smad1")] <- 1
    x[which(ns=="Smad6"), which(ns=="Smad1")] <- 1
    x[which(ns=="Smad6"), which(ns=="Runx1")] <- 1
    x[which(ns=="Smad1"), which(ns=="Runx1")] <- 1
    x <- x + t(x)
    return( x )
}





