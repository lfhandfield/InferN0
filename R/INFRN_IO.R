
#' Load read counts using files produced by 'cellranger count (...)'
#'
#' @param directory path to a cellranger output folder, which expected that contains 3 files named barcodes.tsv genes.tsv matrix.mtx respectively. (ex: "./data/")
#' @param colsubset: if set, getCells with matching name only
#' @return List containing the following: (data, filtered.data, etc)
#'
#' @export
InferN0ReadRangerMatrixCPP <- function(path, threshold.low.UMI=0, threshold.low.GENE=0, threshold.genes.ignored=c("thisgenedoesnotexistsurely"), do.soup.regression=F,do.print.progress=T,do.plot.soup=F, plot.soup.title=c(), threshold.soup.fraction=1.0,threshold.soup.logUMIslope=0,is.floatingpoint=F, do.TCM.reverse=F, do.output.data.only=F, nb.soup.samples=32){
	if ((do.plot.soup)&&(is.null(plot.soup.title))) plot.soup.title <- ""
	if (do.print.progress) flag <- 1
	else flag <- 0	
	if (do.TCM.reverse) flag <- flag + 8
	if (do.output.data.only) flag <- flag + 16
	else if (do.soup.regression) flag <- flag + 2
	if (is.floatingpoint) flag <- flag + 4
	print(flag)
	if (nb.soup.samples == 0) soupsamarg = c(10,0,0)
	else soupsamarg =c(5,25,nb.soup.samples)
	output <- .infernal_readMatrixFolder(list(path=path,low_UMI_threshold=threshold.low.UMI,low_UMI_threshold=threshold.low.GENE,threshold_genes_ignored=threshold.genes.ignored,threshold.soup=threshold.soup.fraction,threshold.soup.logUMIslope=threshold.soup.logUMIslope, souprange=soupsamarg, do.print.progress=flag))
	if (!("data" %in% names(output))){ # no data or all cell filtered!
		return(list());
	}
	if ((do.soup.regression)&&(!is.null(plot.soup.title))&&(!do.output.data.only)){
		layout(matrix(c(1,2),nrow=1,ncol=2,byrow=T))
 		tra <- log10(Matrix::colSums(output$data))
		trb <- log10(Matrix::colSums(output$filtered.data))
		tra2 <- log10(Matrix::colSums(output$data != 0))
		trb2 <- log10(Matrix::colSums(output$filtered.data != 0))
    		xlim = range(output$soupfraction_cell, output$soupfraction_filteredcell)
    		ylim = range(tra,trb)
		ylim2 = range(tra2,trb2)
		plot(output$soupfraction_cell, tra, ,xlim=xlim,ylim=ylim,col="red",xlab="",ylab="" )
		par(new=TRUE)
		plot(output$soupfraction_filteredcell, trb ,xlim=xlim,ylim=ylim, col="gray",xlab="Soup Fraction",ylab="Log10 Total cell UMI",main=plot.soup.title)
		plot(output$soupfraction_cell, tra2, ,xlim=xlim,ylim=ylim2,col="red",xlab="",ylab="" )
		par(new=TRUE)
		plot(output$soupfraction_filteredcell, trb2 ,xlim=xlim,ylim=ylim2, col="gray",xlab="Soup Fraction",ylab="Log10 Gene Coverage",main=plot.soup.title)
	}
return(output)}

#' Convert file list of dense matrices to 3 files for gene name, cell names and sparse data MTX format
#'
#' @param path or list of paths to tab delimited text table
#' @param pathprefix: output prefix for the 3 output files
#' @param cellnameprefix: string appended to cellnames from their respective file
#' @param is_floatingpoint: if set, expects floating point values
#' @param do.unzip: if set, expects the files to have an addionnnal .gz suffix, and temporalily uncompress them 
#'
#' @export
InferN0ConvertFilesToSparse <- function(pathlist, pathprefix, cellnameprefix=c(), is_floatingpoint=F, do.unzip=F,is_commaseparated=F){
	if (is_floatingpoint) flag <- 1
	else flag <- 0
	if (is_commaseparated) flag <- flag + 2

	if (length(cellnameprefix) != length(pathlist)) cellnameprefix = rep("",length(pathlist)) 
	if (do.unzip){
		library(R.utils)
		for(item in pathlist) gunzip(paste(item,".gz",sep=""))
	}
	.infernal_convertDenseTextMatrixToSparseFormat(list(pathlist=pathlist, pathprefix=pathprefix,cellnameprefix=cellnameprefix,flag=flag))
	if (do.unzip){
		for(item in pathlist) gzip(item)				
	}
return}



#' Recover Raw Counts
#'
#' Recovers the 'raw.data' using 'data' assuming cells had their expression normalized (no transcript-length/gene specific normalization)
#'
#' @param path or list of paths to tab delimited text table
#' @param pathprefix: output prefix for the 3 output files
#'
#' @export
InferN0ExportDataAsMtx <- function(infrn0, pathprefix, matname= c(), matrix.filename="matrix.mtx", gene_filename = "genes.tsv", cell_filename = "barcodes.tsv",  use.normalized=F, do.gzip = F){
	if (use.normalized) flag <- 1
	else flag<-0
	if (is.null(matname)){
		if (use.normalized) matname <- "Expression"
		else matname <- "Counts"
	}
	.infernal_exportDataAsMtx(infrn0@ptr, list=list(pathprefix=pathprefix,matname=matname,flag=flag))
	if (do.gzip){
		for(item in paste(pathprefix, c(matrix.filename, gene_filename,cell_filename),sep="")) gzip(item)
	}
return}


#' Load InferN0 Scope
#'
#' Loads InferN0 scope from file.
#'
#' @param path InferN0 scope file path
#' @return InferN0 scope (can only be saved to a file using \code{InferN0SaveFile})
#' @export
InferN0LoadFileOld <- function(path = "./ifrn0.scp"){
	ptr <- .infernal_readScope(path)
	return(new("InferN0Scope",ptr=ptr))
}

#' Load InferN0 Scope
#'
#' Loads InferN0 scope from file.
#'
#' @param path InferN0 scope file path
#' @return InferN0 scope (can only be saved to a file using \code{InferN0SaveFile})
#' @export
InferN0LoadFile <- function(path = "./ifrn0.scp"){
	out <- .infernal_readScope2(path)
	if (is.null(out)) stop()
	outlist <- unserialize(out$serial)
	if (!"gene.meta" %in% names(outlist)) outlist$gene.meta <- data.frame()
	if (!"gene.coor" %in% names(outlist)) outlist$gene.coor <- list()
	if (!"cell.coor" %in% names(outlist)) outlist$cell.coor <- list()
	return(tryCatch(
		new("InferN0Scope",ptr=out$ptr,misc=outlist$misc,meta.data= outlist$meta.data,display=outlist$color,gene.meta.data=outlist$gene.meta, gene.coor=outlist$gene.coor, cell.coor=outlist$cell.coor),
                error=function(cond){print("serialized data mismatch, returning them separated instead.");return(list(scp =new("InferN0Scope",ptr=out$ptr), metalist = outlist))}))
}

#' Save InferN0 Scope
#'
#' Saves InferN0 scope to file. (Needed Previously generated by "loadInferN0file")
#'
#' @param path InferN0 scope file path
#' @param version R serialize version (see ?serialize)
#' @return void
#' @export
InferN0SaveFile <- function(ifrn0, path = "./ifrn0.scp", version = NULL){
	if (class(ifrn0) != "InferN0Scope") stop(paste("InferN0Scope object expected! got", class(ifrn0), "instead!"))
	ptr <- .infernal_writeScope2(ifrn0@ptr, serialize(list(meta.data = ifrn0@meta.data,gene.meta = ifrn0@gene.meta.data, gene.coor= ifrn0@gene.coor, cell.coor= ifrn0@cell.coor,  misc=ifrn0@misc, color=ifrn0@display),NULL, version=version), path)
}

#' Retrieve Zscores associated with observed counts under the learned model
#'
#' Saves InferN0 scope to file. (Needed Previously generated by "loadInferN0file")
#'
#' @param infrn0 InferN0 scope object
#' @param genenames list of genes to output (names or indexes)
#' @param cellnames list of cell to output (names or indexes)
#' @param size.warning get a warning if the retrieved size of the matrix is very large
#' @return dense(genes x cells) matrix of Zscores
#' @export
InferN0GetDeviationMatrix <- function(ifrn0, genenames=c(), cellnames =c(),  size.warning=T){
	if (is.null(cellnames)) {cellnames = ""; nbentry <- length(ifrn0$cell.names)  }
	else nbentry <- length(cellnames)
	if (is.null(genenames)) {genenames = ""; nbentry <- nbentry * length(ifrn0$gene.names)  }
	else nbentry <- nbentry * length(genenames)

	if (size.warning) {
		if (nbentry > 10000000) stop("requested object is very large, retry with 'size.warning=F' to proceed")
	}
return(.infernal_getFullDeviations(ifrn0@ptr, list(genenames=genenames, cellnames=cellnames)))}


#' Convert InferN0 scope to Seurat Object, including custom normalization
#'
#' @param path InferN0 scope file path
#' @param do.hidden.as.pca: uses cell hidden.state to populate sro@dr$pca
#' @return seurat Object
#'
#' @export
InferN0ConvertToSeurat <- function(ifrn0, do.compute.TSNE=F){
	if (class(ifrn0) != "InferN0Scope") stop(paste("InferN0Scope object expected! got", class(ifrn0), "instead!"))
	if (!requireNamespace("Seurat", quietly = TRUE)) {
        	warning("The Seurat package must be installed to use this functionality")
	       	return(NULL)
    	}
	library(Seurat)
	data <- InferN0Get(ifrn0,"data");
	print(class(InferN0Get(ifrn0,"cell.state")))
	hid <- cbind(Matrix::t(InferN0Get(ifrn0,"cell.state")), InferN0Get(ifrn0,"gene.state.deviation"));
	sro <- CreateSeuratObject(InferN0Get(ifrn0,"raw.data"));
	if ("assays" %in% names(attributes(sro))) sro@assays$RNA@data <- data
	else sro@data <- data
	
	#sro@scale.data <- rbind(hid,data);
	sro <- AddMetaData(sro, hid);

	if ("cell.scale" %in% names(ifrn0)) sro <- AddMetaData(sro, metadata=data.frame(InferN0Get(ifrn0,"cell.scale")));
	if ("gene.scale" %in% names(ifrn0)) sro@misc$gene.meta.data <- AddMetaData(sro, metadata=data.frame(InferN0Get(ifrn0,"gene.scale")));
	if ("cell.state" %in% names(ifrn0)){
		drinput <- Matrix::t(ifrn0$cell.state)
		rownames(drinput) <- rownames(sro@meta.data)
		colnames(drinput) <- paste("PC", 1:ncol(drinput),sep="")
		if ("dr" %in% names(sro)) sro@dr[["ifrn0state"]] <- new("dim.reduction",cell.embeddings=as.matrix(drinput),key = "PC")
		else sro@reductions[["ifrn0state"]] <- new("DimReduc",cell.embeddings=as.matrix(drinput),key = "PC",assay.used=DefaultAssay(sro))
		if (do.compute.TSNE){
			sro <- RunTSNE(sro, reduction = "ifrn0state", reduction.name="infrTSNE", dims = 1:ncol(drinput),dim.embed=2)
		}
	}
	sro@misc <- ifrn0@misc
return(sro)}

#' Convert InferN0 scope to Singlecell Experiement Object, including custom normalization
#'
#' @param path InferN0 scope file path
#' @param do.hidden.as.pca: uses cell hidden.state to populate sro@dr$pca
#' @return seurat Object
#'
#' @export
InferN0ConvertToSingleCellExperiment <- function(ifrn0){
	if (class(ifrn0) != "InferN0Scope") stop(paste("InferN0Scope object expected! got", class(ifrn0), "instead!"))
	if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
        	warning("The SingleCellExperiment package must be installed to use this functionality")
	       	return(NULL)
    	}
	library(SingleCellExperiment)
	sce <- SingleCellExperiment(assays=list(counts=ifrn0$rawdata),reducedDims= SimpleList(ifrn0@cell.coor), colData= ifrn0@cell.meta , rowData = ifrn0@gene.meta.data, metadata = ifrn0@misc)
return(sce)}





