
#' Initialization of InferN0Scope
#'
#' Initialize Infern0 scope from input data, which can be a matrix, a singlecelleperiment or Seurat object
#'
#' @return InferN0 scope object (can only be saved to a file using \code{infern0SaveFile})
#' @param input versatile input data, matrix with gene/feature as rows and cells as collunms, seurat object or SingleCellExperiment object
#' @param use.normalized.data selects if the input populates the rawdata (integer counts) or the data slot for normalized data (floating-point value)
#' @param meta.data meta data associated with each cell
#' @export
infern0Init <- function(input, meta.data= c(), use.normalized.data=F){
	ptr <- .infernal_createScope()
	if (class(input) == "SingleCellExperiment"){
		library("SingleCellExperiment")
		if (class(counts(input)) == "matrix") .infernal_readDenseMatrix(ptr,counts(input), rawdata)
		else if (isS4(counts(input)) & is(counts(input),'sparseMatrix') ) .infernal_readSparseMatrix(ptr,as(counts(input),"dgCMatrix"),rawdata)
		scp <- new("InferN0Scope",ptr=ptr)
	}else{
	scp <- new("InferN0Scope",ptr=ptr)
	 if ((class(input) == "seurat")|(class(input) == "Seurat")){
		if ("dr" %in% attributes(input)) drname <- "dr"
		else drname <- "reductions"
		for(i in names(slot(sro,"reductions"))) {
			scp@cell.coor[[i]] <- slot(sro,"reductions")[[i]]@cell.embeddings
		}
		scp@misc <- sro@misc
		if (is.null(meta.data)) meta.data <- input@meta.data
		if (!use.normalized.data) {
			if ("raw.data" %in% names(attributes(input))) input <- input@raw.data
			else if ("RNA" %in% names(input@assays)) input <- input@assays[["RNA"]]@counts
			else input <- input@assays[[input@active.assay]]@counts
		}else {
			if ("data" %in% names(attributes(input))) input <- input@data
			else if ("RNA" %in% names(input@assays)) input <- input@assays[["RNA"]]@data
			else input <- input@assays[[input@active.assay]]@data
		}
	} 
	if (class(input) == "data.frame") input <- as.matrix(input)
	if (is.null(rownames(input))) rownames(input) <- paste("GENE", 1:nrow(input),sep="_")
	if (is.null(colnames(input))) colnames(input) <- paste("CELL", 1:ncol(input),sep="_")

	if (class(input) == "matrix") {if (use.normalized.data) scp$data <- input else scp$raw.data <- input}
	else if (isS4(input) & is(input,'sparseMatrix') ) {if (use.normalized.data) scp$data <- as(input,"dgCMatrix") else scp$raw.data <- as(input,"dgCMatrix") }
	else {
		stop(paste("'",class(input) ,"' is not a recognised type",sep=""))
	}
	}

	scp@meta.data <- data.frame(row.names= scp$cell.names)
	scp@gene.meta.data <- data.frame(row.names= scp$gene.names)
	scp@cell.coor <- list()
	scp@gene.coor <- list()

	if (!is.null(meta.data)) {
		if (is.null(rownames(meta.data))){
			if (nrow(meta.data) != length(scp$cell.names)) map <- c()
			else map <- 1:length(scp$cell.names)
		}else map <- match(scp$cell.names, gsub("-", ".", gsub(" ", ".",rownames(meta.data))))
		if (is.null(colnames(meta.data))) colnames(meta.data) <- paste("meta.data", 1:ncol(meta.data))
		if ((!is.null(map))&(sum(!is.na(map)) > 0)&(ncol(meta.data) >0)){
			for(i in 1:ncol(meta.data)) {
				if (class(meta.data[[i]]) == "factor"){
					scp@meta.data[[colnames(meta.data)[i]]] <- factor(NA, levels = levels(meta.data[[i]]))
					scp@meta.data[[colnames(meta.data)[i]]]@.Data[!is.na(map)] <- meta.data[[i]]@.Data[map[!is.na(map)]]
				}else{
					scp@meta.data[[colnames(meta.data)[i]]][!is.na(map)] <- meta.data[map[!is.na(map)],i]
				}
			}
		}else{
			warning("Cell names mismatch in provided metadata. No metadata is imported")
		}
	}
	return(scp)
}



#' Alias of infern0Init for conversion operation
#'
#' @param versatile input data, matrix with gene/feature as rows and cells as collunms, seurat object or SingleCellExperiment object
#' @return InferN0 scope
#' @export
as.InferN0Scope <- function(object){return(infern0Init(object))}

