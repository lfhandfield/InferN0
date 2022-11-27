


#' @export
InferN0TestCpp <- function(){
return(.infernalTest())}


#' Example function to analyse data using provided function. This function serves as a guideline on the intended use of other functions.
#'
#' @param Seurat Object/SingleCellExperiment/data matrix with gene/feature as rows and cells as collunms
#' @param network_genelist list of genes for network inference
#' @param nbgene.to.detect number of co-expressed genes to detect, ("default" is a number equal to the length of 'network_genelist')
#'
#' @export
InferN0RunExample <- function(data, network_genelist, nbgene.to.detect="use.default",nb.threads=4){
	out <- list()
	out$infscp <- InferN0Init(data)
	out$ModelData <- InferN0ModelData(out$infscp,nb.hidden.row=8,nb.hidden.col=5,nb.step=25,nb.thread=nb.threads);
	out$IdentifyNetwork <- InferN0IdentifyNetworkPipeline(out$infscp,network_genelist, nb.threads=nb.threads);
return(out)}

#' Displays a summary of a InferN0 scope
#'
#' @param path InferN0 scope file path
#' @return void
#' @export
InferN0Summary <- function(ifrn0){
	.infernal_show(ifrn0@ptr);
}

#' Unsupported function for internal testing purposes
#'
#' @param data matrix with gene/feature as rows and cells as collunms
#'
#' @export
.InferN0Test<- function(scope){.infernal_runTestFunction(scope@ptr)}
#' Fit Negative Binomial to data
#'
#' @param InferN0 scope
#'
#' @export
.InferN0FitNegativeBinomial <- function(scope){
	.infernal_fitNegativeBinomial(scope@ptr)
}

#' (...)
#'
#' @param InferN0 scope
#' @param path prefix for output files (ex: "~/outdir/output")
#'
#' @export
InferN0SaveHierarchical <- function(scope, pathprefix){.infernal_saveHierarchical(scope@ptr,pathprefix)}

#' Execute HierarchicalClustering
#'
#' @param InferN0 scope
#'
#' @export
InferN0HierarchicalClusteringWithInput <- function(ifrn0, cells.group){
#	out <- internal_InferN0_hierarchicalClustering_CustomGroups(ifrn0$data, cells.group, max(cells.group));
return(c(ifrn0, list(gene.tree=out)))}

#' Execute HierarchicalClustering
#'
#' @param InferN0 scope
#' @param cluster.cells , TRUE to  cluster cell, FALSE to cluster GENES
#' @param (over)writes in InferN0 scope the ordering as the new hierarchical clustering indicates
#'
#' @export
InferN0HierarchicalClustering <- function(ifrn0, cluster.cells=TRUE, write.ordering=TRUE){
out <- .infernal_HierarchicalClustering(ifrn0@ptr, list(cluster_cells=cluster.cells, write_ordering=write.ordering));
return(out)}

#' Find optimal Covariance with associated constrained Precision matrix
#'
#' @param symmetric sparse matrix constraint
#' @param list equivalent to what is outputed by the 'InferN0ComputeCovariance' function (if used, genes must agree)
#' @return list
#'
#' @export
InferN0FindConstrainedCovariance <- function(constraint, covarstr){
	outlist <- .infernal_IdentifyConstrainedCovar(covarstr$Input, list(covarstr$Covar,covarstr$Mean,partition=covarstr$Testset.Partition,constraint=constraint), covarstr$Training.Covar ,  covarstr$Training.Mean)


	tryCatch({
		outlist$LLDerivative = solve(outlist$Precision) - covarstr$Covar
		outlist$LLDerivative[outlist$Precision == 0] <- 0
	},error=function(cond){print("Warning, Obtained Precision Matrix is Singular!")})

	outlist$LLtable <- matrix(0,3+length(covarstr$Training.Mean),4)
	outlist$LLtable[1,1] <- 1
	outlist$LLtable[1,3] <- determinant(outlist$Precision,logarithm =T)$modulus
	outlist$LLtable[1,4] <- sum(outlist$Precision * covarstr$Covar)
	outlist$LLtable[1,2] <- (outlist$LLtable[1,3] - outlist$LLtable[1,4] - log(pi * 2.0)) * 0.5

	outlist$LLtable[2,1] <- length(covarstr$Testset.Partition)
	outlist$LLtable[2,3] <- determinant(outlist$Precision,logarithm =T)$modulus
	outlist$LLtable[2,4] <- sum(outlist$Precision * outlist$SuffitientMatrix) / outlist$LLtable[2,1]
	outlist$LLtable[2,2] <- (outlist$LLtable[2,3] - outlist$LLtable[2,4] - log(pi * 2.0)) * 0.5
	for(i in 1:length(covarstr$Training.Mean)){ 
		outlist$LLtable[3+i,1] <- sum(covarstr$Testset.Partition == (i-1))
		outlist$LLtable[3+i,3] <- determinant(outlist$TrainingPrecisionMatrices[[i]],logarithm =T)$modulus
		outlist$LLtable[3+i,4] <- sum(outlist$TrainingPrecisionMatrices[[i]] * outlist$TestSetSuffitientMatrices[[i]]) / outlist$LLtable[3+i,1]
		outlist$LLtable[3+i,2] <- (outlist$LLtable[3+i,3] - outlist$LLtable[3+i,4] - log(pi * 2.0)) * 0.5
	}

	outlist$LLtable[3,] <- colSums(outlist$LLtable[4:(length(covarstr$Training.Mean)+3),]) / length(covarstr$Training.Mean)
	colnames(outlist$LLtable) <- c("n", "LL over n", "LogDeterminent", "TestSetDeviationTerm")
	rownames(outlist$LLtable) <- c("Covariance Data", "Whole Data", "Testset average", paste("Testset", as.character(1:length(covarstr$Training.Mean)) ,sep=""))
return(outlist)}

#' Export hierarchical clustering(s) to target directory
#'
#' @param InferN0 scope
#' @param path to directory for output tree files
#'
#' @export
InferN0_SaveAsCDTfile <- function(ifrn0, path){out <- infernal_saveHierarchical(ifrn0@ptr, path);}

#' Internal: filter maximimal numb of row/col to retain percentage of data
#'
#' @param SingleCellExperiement Object to load
#'
#' @export
InferN0internal_FilterRowCol <- function(input, fraction = 0.9, cost_ratio = -1, do_print = FALSE){
    dmat <- as.matrix(input) != 0
    tsum <- sum(dmat)
    target <- floor(tsum * (1 - fraction))
    idim <- dim(dmat)
    rowf <- rep(TRUE, idim[1])
    colf <- rep(TRUE, idim[2])
    rite <- 1
    cite <- 1
    nbfilt <- 0
    ctar <- 0

    if (cost_ratio < 0) cost_ratio <- idim[1] / idim[0] # illegal here, use aspect ratio instead!

    if (target >= tsum){
        return(list(row_filter= rep(FALSE, idim[1]), col_filter = rep(FALSE, idim[2])))
    }

    while(target > ctar){
        ctar <- (ctar + target+1) /2
        csum <- colSums(dmat)
        rsum <- rowSums(dmat)
        cord <- order(csum)
        rord <- order(rsum)
        while(nbfilt < ctar){
            if (csum[cord[cite]] < rsum[rord[rite]] * cost_ratio){
                colf[cord[cite]] = FALSE;
                nbfilt <- nbfilt + sum(dmat[,cord[cite]])
                dmat[,cord[cite]] <- 0
                cite <- cite + 1
            }else{
                rowf[rord[rite]] = FALSE;
                nbfilt <- nbfilt + sum(dmat[rord[rite],])
                dmat[rord[rite],] <- 0
                rite <- rite + 1
            }
        }

    }
    if (do_print){
        print(paste(sum(rowf), "/",idim[1]," rows & ", sum(colf), "/",idim[2]," collumns remains", sep=""))
    }

    return (list(row_filter= rowf, col_filter = colf))
}

#' Recover Raw Counts
#'
#' Recovers the 'raw.data' using 'data' assuming cells had their expression normalized (no transcript-length/gene specific normalization)
#'
#' @param path or list of paths to tab delimited text table
#' @param pathprefix: output prefix for the 3 output files
#'
#' @export
InferN0ReverseSizeFactor <- function(infrn0, do.log=F){
	if (do.log) flag <- 1
	else flag<-0
	.infernal_reverseTCM(infrn0@ptr, list=list(flag=flag))
return}

#' Convert file list of dense matrices to 3 files for gene name, cell names and sparse data MTX format
#'
#' @param path or list of paths to tab delimited text table
#' @param pathprefix: output prefix for the 3 output files
#'
#' @export
.InferN0ReverseTCMinMTX <- function(inprefix, outprefix, do.log=F){
	if (do.log) flag <- 1
	else flag<-0
	.infernal_reverseTCM_inMTX(list(inprefix=inprefix,outprefix=outprefix,flag=flag))
return}

#' Function that generates synthetic read counts, assuming there are "nb.hidden.col" cell types with "nb.hidden.row" subset of differentially expressed genes(or features)
#'
#' @param ifrn0 InferN0 scope
#' @param nb.row number of features
#' @param nb.col number of cells
#' @param nb.hidden.row number of differentially expressed features groups
#' @param nb.hidden.col number of cell types
#'
#' @export
InferN0ExportPrincipalComponents <- function(ifrn0,nb.pc=50, nb.thread=4){
	return(.infernal_exportpcs(ifrn0@ptr, list(nb_pc=nb.pc,nb_thread=nb.thread,nbstep= 10, flags=1)))
}

#' Merge sparse matrices using matching row/column names
#'
#' @param matrixA: first matrix to merge
#' @param matrixB: second matrix to merge
#' @param mode: operation "o" overwrite non-zero, "O" overwrite
InferN0SparseMatrixMerge <- function(matrixA, matrixB, mode= "o"){
	return(.infernal_mergeSparseMatrices(matrixA,matrixB))
}

#' Compute Wilcox test adjusted for missing data (Hypergeometric-Normal mixture for rank sum test)
#'
#' @param input.matrix: sparse matrix or InferN0 Scope with transposed matrix...
#' @param list.pos: list of rows for Negative Partition, or T/F vector matching dimention
#' @param list.neg: list of rows for Positive Partition, or T/F vector matching dimention
InferN0ZeroCorrectedWilcox <- function(input, list.pos, list.neg, do.transpose=T, do.cmp.mean=F, use.clusterID=F,purewilcox=F,print.progress=T){
	if (do.cmp.mean) flag <- 2
	else flag <- 0
	if (print.progress) flag <- flag + 8
	if (class(input) == "InferN0Scope"){
		nbcell <- length(input$gene.names)
		if (use.clusterID){
			clustID <- input$gene.cluster
			list.neg <- clustID %in% list.neg
			list.pos <- clustID %in% list.pos
		}		
		if (length(list.neg) == nbcell) list.neg <- (1:nbcell)[list.neg];
		if (length(list.pos) == nbcell) list.pos <- (1:nbcell)[list.pos];
		if (purewilcox) flag <- flag + 16
		return(.infernal_wilcox_scope(input@ptr,list(list.pos=list.pos,list.neg=list.neg,do.transpose=flag)))
	}else{
		if (do.transpose) {flag <- flag + 1; nbcell <- dim(input)[2]}
		else{ nbcell <- dim(input)[2]}


		if (length(list.neg) == nbcell) list.neg <- (1:nbcell)[list.neg];
		if (length(list.pos) == nbcell) list.pos <- (1:nbcell)[list.pos];

		return(.infernal_wilcox_matrix(input@ptr,list(list.pos=list.pos,list.neg=list.neg,do.transpose=flag)))
	}
}

#' Compute Wilcox test adjusted for missing data (Hypergeometric-Normal mixture for rank sum test)
#'
#' @param input.matrix: InferN0 Scope with transposed matrix...
#' @param list.pos: list of rows for Negative Partition, or T/F vector matching dimention
#' @param list.neg: list of rows for Positive Partition, or T/F vector matching dimention
#' @param total_count: cells read depth, used to partition into quartiles and normalize fold change
#'
#' @export
InferN0ZeroCorrectedWilcox2 <- function(input, list.pos, list.neg, total_count, nb.threads=4, list.partition= c(), nbpart=4, use.clusterID=F,print.progress=T,do.quartile.average=T, gene.use=c(), do.hypergeo.correct=T,do.plotZscores=F, do.downsample=F, do.output.overwritematrix=F, do.fdr.correct=T){
	if (print.progress) flag <- 1
	else flag <- 0
	ordering <- order(total_count)
	if (do.hypergeo.correct) flag <- flag + 2
	if (do.downsample) flag <- flag + 64
	if (do.output.overwritematrix) flag <- flag + 128
	if (do.fdr.correct) flag <- flag + 256
	if (class(input) == "InferN0Scope"){
		nbcell <- length(input$gene.names)
		if (use.clusterID){
			clustID <- input$gene.cluster
			list.neg <- clustID %in% list.neg
			list.pos <- clustID %in% list.pos
		}		
		if (length(list.neg) == nbcell) list.neg <- (1:nbcell)[list.neg];
		if (length(list.pos) == nbcell) list.pos <- (1:nbcell)[list.pos];
		if (is.null(gene.use)) gene.use = 1:length(input$cell.names)
		if (is.null(list.partition)) list.partition = c(1)

		output <- .infernal_find_DE(input@ptr,list(list.pos=list.pos,list.neg=list.neg,order=total_count,partition= list.partition, nbpart=nbpart,nbthread=nb.threads, do.transpose=flag,glist=gene.use))
		if (!do.quartile.average) {
			output$Zscore <- output$Zscore[,((1:nbpart)*2 -1)];
			output$LogitAuroc <- output$LogitAuroc[,((1:nbpart)*2 -1)];
			return(output)
		}
		denP <- mean(total_count[list.pos])
		denN <- mean(total_count[list.neg]) 
		if (do.plotZscores != 0){
			#print(plot(output$Zscore[,do.plotZscores*2-1],output$LogitAuroc[,do.plotZscores*2-1]))
			hist(output$Zscore, breaks=100)
			for(i in 1:ncol(output$Zscore)) print(paste(var(output$Zscore[,i])))
		}
		col1 <- rowSums(output$Zscore[,((1:nbpart)*2 -1)] * sqrt(output$Weight[,((1:nbpart)*2 -1)])   ,na.rm=T)
		col2 <- rowSums(output$Weight[,((1:nbpart)*2 -1)])
		#col2 <- rowMeans(output$Zscore[,((1:nbpart)*2)],na.rm=T)
		output$Zscore <- col1 / sqrt(col2)
	
		col1 <- rowSums(output$CoverageEnrichment[,((1:nbpart)*2 -1)] * output$Weight[,((1:nbpart)*2 -1)] ,na.rm=T)
		#col2 <- rowMeans(output$LogitAuroc[,((1:nbpart)*2)],na.rm=T)
		output$CoverageEnrichment <- col1 / col2

		col1 <- rowSums(output$LogitAuroc[,((1:nbpart)*2 -1)] * output$Weight[,((1:nbpart)*2 -1)],na.rm=T)
		#col2 <- rowMeans(output$LogitAuroc[,((1:nbpart)*2)],na.rm=T)
		output$LogitAuroc <- col1 / col2 #cbind(col1,col2)
		output$Weight <- col2
		col1 <- rowMeans(output$Mean[,((1:nbpart)*2 -1)],na.rm=T) / denP
		col2 <- rowMeans(output$Mean[,((1:nbpart)*2)],na.rm=T) / denN
			
		output$MeanTPM <- (col1 + col2) * 500000
		output$log2FC <- log2(col1) - log2(col2)

		return(output)
	}else stop("InferN0 scope expected")
}

#' Compute Wilcox test adjusted for missing data (Hypergeometric-Normal mixture for rank sum test)
#'
#' @param value: Matrix with measured feature that is to be ranked
#' @param pvalue: pvalue associated with observation
#'
#' @export
InferN0DElistCompare <- function(value, pvalue, pvaluethreshold=0.05){
	value[pvalue < pvaluethreshold] <- 0;
	input <- Matrix(value,sparse=T)

}

#' Compute Wilcox test adjusted for missing data (Hypergeometric-Normal mixture for rank sum test)
#'
#' @param input.matrix: InferN0 Scope with transposed matrix...
#' @param partition: list of positive integer with cluster identity of cells
#' @param ordering: cells ordered by major counfounding factor, used to partition into quartiles
#' @param is.ordering.a.partition: defines the partition directly
#' @param use.pure.wilcox: use pure wilcox test
#'
#' @export
InferN0FindMarkers <- function(input, nbpart=1, partition=c(), ordering=c(), print.progress=T,do.quartile.average=T, is.ordering.a.partition=F, use.pure.wilcox =F, do.downsample=F, cell.use = c()){
	if (print.progress) flag <- 1
	else flag <- 0
	if (use.pure.wilcox) flag <- flag + 2
	if (do.downsample) flag <- flag + 64
	if (is.null(cell.use)) cell.use = rep(T, length(input$gene.names))
	if (class(input) == "InferN0Scope"){
		nbcell <- length(input$gene.names)
		if (is.null(partition)) partition = input$gene.cluster;
		output <- .infernal_find_markers(input@ptr,list(partition=partition,order=ordering,nbpart=nbpart,do.transpose=flag, cell.use));
		nbcluster <- dim(output$Zscore)[2] / nbpart;
		if (nbpart != 1){
			swap <- matrix(0,dim(output$Zscore)[1], nbcluster)
			rownames(swap) <- rownames(output$LogitAuroc)
			colnames(swap) <- colnames(output$LogitAuroc)[1:nbcluster]
			for(i in 1:nbcluster) swap[,i] <- rowMeans(output$LogitAuroc[,((0:(nbpart-1))*nbcluster +i)],na.rm=T)
		}else swap <-output$LogitAuroc 
		output$MarkPartition <- array(0,dim(output$LogitAuroc)[1])
		for(i in 1:length(output$MarkPartition)) {
			tmp <- which(swap[i,] == max(swap[i,]))
			if (length(tmp) == 1) output$MarkPartition[i] <- tmp
			else output$MarkPartition[i] <- tmp[1]
		}
		if ((!do.quartile.average)||(nbpart == 1)) return(output)
		output$LogitAuroc <- swap;
		for(i in 1:nbcluster) swap[,i] <- rowMeans(output$Zscore[,((0:(nbpart-1))*nbcluster +i)],na.rm=T)
		output$Zscore <- swap;

		for(i in 1:nbcluster) swap[,i] <- rowMeans(output$Mean[,((0:(nbpart-1))*nbcluster +i)],na.rm=T)
		output$Mean <- swap;
		return(output)
	}else stop("InferN0 scope expected")
}

#' Function that generates synthetic read counts, assuming there are "nb.hidden.col" cell types with "nb.hidden.row" subset of differentially expressed genes(or features)
#'
#' @param nb.row number of features
#' @param nb.col number of cells
#' @param nb.hidden.row number of differentially expressed features groups
#' @param nb.hidden.col number of cell types
#' @return InferN0 scope
#'
#' @export
InferN0GenerateSyntheticData <- function(nb.row=320,nb.col=200,nb.hidden.row=8,nb.hidden.col=5,mixedstate_weight=1){
	ptr <- .infernal_createScope();
	.infernal_genSynthetic(ptr, list(nbrow=nb.row,nbcol=nb.col,nbhrow=nb.hidden.row,nbhcol =nb.hidden.col, mixed= mixedstate_weight));

	return(new("InferN0Scope",ptr=ptr))
}

#' Generate Synthetic Data with a given Network
#'
#' @param Z: Symetrix matrix as network
#' @param n: number of generated cells (collumns) 
#' @param mu.min,mu.max: mean expression of transcripts
#' @param p.min,p.max: dropout rate transcripts
#' @param s.min,s.max: background noise for expression of transcripts
InferN0GenerateSyntheticNetworkData <- function(Z, n=100, mu.min = 10, mu.max = 10, p.min=.3, p.max=.3, s.min=3, s.max=3, noise.min=0, noise.max=0, muvec = c(), pvec = c(), svec = c(),noisevec=c(),  do.log.transform=T, is.log.dependency=T) {
    require('MASS')
    #Create the covariance matrix
    Z <- (Z != 0)
    sum <- rowSums(abs(Z))
    S <- diag(rowSums(abs(Z))) + Z


    #Generate synthetic data
    cond <- 0
    d <- dim(S)[1]

    for(i in 1:d){
        S[,i] <- S[,i] / sqrt(sum[i])
        S[i,] <- S[i,] / sqrt(sum[i])
    }
    S <- solve(S) # got to inverse to get Covariance parameter
    if (length(muvec) != d) muvec <- sample(runif(d) * (mu.max- mu.min) + mu.min)
    if (length(pvec) != d) pvec <- sample(runif(d) * (p.max- p.min) + p.min):w
    if (length(svec) != d) svec <- sample(runif(d) * (s.max- s.min) + s.min)
    if (length(noisevec) != d) noisevec <- sample(runif(d) * (noise.max- noise.min) + noise.min)
    
    if (is.log.dependency) {
	muvec <- log2(muvec)
        svec <- log2(svec)
    }

    S <- svec * S * svec
    x <- mvrnorm(n, muvec, S) + mvrnorm(n, rep(0, d), (noisevec %*% t(noisevec)) +diag(noisevec))
    #Remove random entries as dropouts
    for (i in 1:dim(x)[1]) {
        for (j in 1:dim(x)[2]) {
            if (runif(1)<pvec[j])  x[i, j] <- 0
            else if (is.log.dependency){
                x[i, j] <- ceiling(2^x[i,j]-1)
            }else if (x[i,j] <= 0) x[i,j] = 0
	    else x[i,j] <- ceiling(x[i,j])

        }
    }
	if (do.log.transform) return( t(log2(1+x)) )
	else return( t(x) )
}

# library(InferN0); scptr <- InferN0LoadFile("InferN0output/prefrontalcortex.ifn.scp"); genelist = c("BCL11B","CUX1","CUX2","ETV1","FEZF2","NR4A2","POU3F2","RELN","RORB","SATB2","SOX5","CNTNAP2","CALB1","CALB2","CCK","GAD1","GAD2","NPY","SST","TH","ERBB4","MEF2C","NPAS3","GRIA1","NLGN1","DLX6-AS1","VIM","APOE","APP","PSEN1","STMN2","ROBO1","TOP2A","HIST1H4C","H3F3B","MALAT1","RPLP0","RPL41","RPL12");cov <- InferN0ComputeCovariance(scptr,gene.list=genelist)
# library(InferN0); sro <- readRDS("frozenmerged3.sro.rds"); InferN0SparseMatrixMerge(sro@raw.data[1:100,1:100], sro@raw.data[51:150,51:150]*10) 


#' Example function to analyse data using provided function. This function serves as a guideline on the intended use of other functions.
#'
#' @param 2-collumn data points 
#' @param network_genelist list of genes for network inference
#' @param nbgene.to.detect number of co-expressed genes to detect, ("default" is a number equal to the length of 'network_genelist')
#'
#' @export
InferN0GetDensity <- function(coor, weight=c(),data=c(), covariance=c(), bandwidth= 0.25, mapsize=256, rect=c(), do.normalize.bandwidth=T){
	if ((is.null(dim(coor)))||(dim(coor)[2] != 2)) stop("Expected input is a matrix with tw collumn")
	if (is.null(rect)){
		rect <- rep(0,4);
		tmp <- range(coor[,1],na.rm=T);	rect[1] <- tmp[1];rect[3] <- tmp[2];
		tmp <- range(coor[,2],na.rm=T);	rect[2] <- tmp[1];rect[4] <- tmp[2];
	}else{
		
	}
	if (is.null(covariance)){covariance = cov(coor)}
	covariance = covariance * (bandwidth*bandwidth)

	if (!is.null(weight)) {
		coor <- cbind(coor,weight)
		flag = 1
	}else flag =0
	if (!is.null(data)){
		coor <- cbind(coor,data)
		flag = flag +2
	}
	if (do.normalize.bandwidth) flag <= flag +4 
return(.infernalDensity(coor,list(flag=flag,dims= c(mapsize,mapsize), rect=rect, covar=c(covariance[1,1], covariance[1,2],covariance[2,2]) )))}

callCumulant <- function(op, arga, argb, quartile, x, fake = 0,order=16){
return(.infernalCumulant(list(args= c(arga, argb, quartile, x), opsel = op, fake =fake, order=order)))
}



#' Example function to analyse data using provided function. This function serves as a guideline on the intended use of other functions.
#'
#' @param 2-collumn data points 
#' @param network_genelist list of genes for network inference
#' @param nbgene.to.detect number of co-expressed genes to detect, ("default" is a number equal to the length of 'network_genelist')
#'
#' @export
InferN0Deconvolve2D <- function(X, Y, celltypes, donors, genenames, cell_filter = c(), gene_filter = c(),minimum_countpergene = 100, nbstep = 100, save_path="", do_load=F, use.old=T, use.intersept=F, use.oldbulkdeconv=F, do.gpmodel=F){
	Y <- as(Y,"dgCMatrix")
	if (type(Y) == "double"){
		print("fractionnal counts detected! sampling the fractionnal counts away to get integers")
		Y@x = as.integer(Y@x) + as.integer(runif(length(Y@x)) < (Y@x %% 1.0))
		print(type(Y@x))
	}
	if (is.null(colnames(Y))) colnames(Y) <- paste("bulk", 1:ncol(Y))

	if (is.null(gene_filter)) gene_filter <- rep(1, nrow(Y))

	tcount <- Matrix::rowSums(X)
	gene_filter[tcount < minimum_countpergene] = 0;

	if (!is.null(cell_filter)){
		X <- as(X,"dgCMatrix")
		X <- X[, cell_filter]
		print(paste("filtering", sum(!cell_filter),"cells"))
		celltypes <- droplevels(as.factor(celltypes[cell_filter]))
		donors <- droplevels(as.factor(donors[cell_filter]))
	}else X <- as(X,"dgCMatrix") 
	if (do_load) flag <- 1
	else flag <- 0
	if (use.old) flag <- flag + 2
	if (use.intersept) flag <- flag +4
	if (use.oldbulkdeconv) flag <- flag + 8
	if (do.gpmodel) flag <- flag + 16
	return( .infernalDeconvolve2D(X,Y,list(ctlist = celltypes@.Data -1, dnlist = donors@.Data - 1, ctnames = levels(celltypes), dnnames = levels(donors), genenames = rownames(X), gene_filter=gene_filter, nbstep=nbstep, save_path=save_path, flag=flag )))
}

#' Example function to analyse data using provided function. This function serves as a guideline on the intended use of other functions.
#'
#' @param SingleCellExperiement 
#' @param colData celltype key
#' @param colData sameple key
#'
#' @export
InferN0DualDE <- function(sce, obskey_celltypes, obskey_samples){
	library(DESeq2)
	
}



