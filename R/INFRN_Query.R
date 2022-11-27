
#' Query variable within InferN0Scope
#'
#' @param InferN0 scope
#' @param InferN0 variable name, one of ("cell.cluster", "gene.cluster", "cell.scale", "cell.order","gene.scale","gene.order")
#' @return requested attribute of InferN0 scope
#'
#' @export
InferN0Get <- function(scope, varname){ 
	if (is.null(scope)) stop("InferN0 scope is not initialized!")
	if (is.null(scope@ptr)) stop("InferN0 scope is not initialized!")
	if (varname == "cell.cluster") return(.infernal_getParameter(scope@ptr, 0))
	if (varname == "gene.cluster") return(.infernal_getParameter(scope@ptr, 1))
	if (varname == "cell.scale") return(.infernal_getParameter(scope@ptr, 2))
	if (varname == "gene.scale") return(.infernal_getParameter(scope@ptr, 3))
	if (varname == "cell.state") return(.infernal_getParameter(scope@ptr, 4))
	if (varname == "gene.state") return(.infernal_getParameter(scope@ptr, 5))
	if (varname == "cell.order") return(.infernal_getParameter(scope@ptr, 8))
	if (varname == "gene.order") return(.infernal_getParameter(scope@ptr, 9))
	if (varname == "cell.names") return(.infernal_getParameter(scope@ptr, 10))
	if (varname == "gene.names") return(.infernal_getParameter(scope@ptr, 11))
	if (varname == "cell.nbcounts") return(.infernal_getParameter(scope@ptr, 21))
	if (varname == "gene.nbcounts") return(.infernal_getParameter(scope@ptr, 22))


	if (varname == "raw.data") return(.infernal_getParameter(scope@ptr, 6))

	if (varname == "data") return(.infernal_getParameter(scope@ptr, 7))

	if (varname == "raw.data") return(.infernal_getParameter(scope@ptr, 6))
	if (varname == "data") return(.infernal_getParameter(scope@ptr, 7))
	if (varname == "cell.state.deviation") return(.infernal_getParameter(scope@ptr, 12))
	if (varname == "gene.state.deviation") return(.infernal_getParameter(scope@ptr, 13))
	if (varname == "state.matrices") return(.infernal_getParameter(scope@ptr, 14))
	if (varname == "cell.coverage") return(.infernal_getParameter(scope@ptr, 15))
	if (varname == "gene.coverage") return(.infernal_getParameter(scope@ptr, 16))
	if (varname == "pval.data") return(.infernal_getParameter(scope@ptr, 17))
	if (varname == "dropprob.data") return(.infernal_getParameter(scope@ptr, 18))
	if (varname == "cell.state.transition") return(.infernal_getParameter(scope@ptr, 23))
	if (varname == "gene.state.transition") return(.infernal_getParameter(scope@ptr, 24))

	stop(paste(varname, "is not an attribute of an InferN0 scope"))
}

#' (...)
#'
#' @param InferN0 scope
#' @param list/matrix with clusterID as values
#'		
#' @export
InferN0SetClusterIdentity <- function(scope, assosication){.infernal_setParameter(scope@ptr,0,assosication)}

#' Set variable within InferN0Scope
#'
#' @param InferN0 scope
#' @param InferN0 variable name, one of ("cell.cluster", "gene.cluster", "cell.scale", "cell.order","gene.scale","gene.order")
#' @param value to set, its type must match the output of 'InferN0Get' with matching variable name
#' @return requested attribute of InferN0 scope
#'
#' @export
InferN0Set <- function(scope, varname, value){ 
	expectedclass <- switch(varname
, "cell.state" = "readonly"
, "gene.state" = "readonly"
, "cell.scale" = "readonly"
, "gene.scale" = "readonly"
, "state.matrices" = "readonly"
, "cell.coverage"= "readonly"
, "gene.coverage"= "readonly"
, "cell.nbcounts"= "readonly"
, "gene.nbcounts"= "readonly"
, "gene.state.transition" = "readonly"
, "cell.state.transition" = "readonly"
, "cell.cluster" = c("list", "integer")
, "gene.cluster" = c("list", "integer")
, "cell.order" = c("list", "integer")
, "gene.order" = c("list", "integer")
, "cell.names" = c("list", "character")
, "gene.names" = c("list", "character")
, "data"= "custom"
, "raw.data"= "custom"
, "unknown"
)
	if (length(expectedclass) == 1)	{
		if (expectedclass == "readonly") stop(paste(varname, "cannot be altered"))
		else if (expectedclass == "unknown") stop(paste(varname, "is not a part of InferN0 scope"))
	}
	if (is.null(scope)) stop("InferN0 scope is not initialized!")
	if (is.null(scope@ptr)) stop("InferN0 scope is not initialized!")
	if (varname == "raw.data") {
		if (class(value) == "matrix") {
			if (length(rownames(value)) != dim(value)[1]) rownames(value) <- paste( "gene", 1:dim(value)[1])
			else rownames(value) <- make.names(rownames(value),unique=T)
			if (length(colnames(value)) != dim(value)[2]) colnames(value) <- paste( "cell", 1:dim(value)[2])
			else colnames(value) <- make.names(colnames(value),unique=T)
			outcode <- .infernal_readDenseMatrix(scope@ptr,value, T)
                }else if (isS4(value) & is(value,'sparseMatrix') ) {
			if (length(rownames(value)) != dim(value)[1]) rownames(value) <- paste( "gene", 1:dim(value)[1])
			else rownames(value) <- make.names(rownames(value),unique=T)
			if (length(colnames(value)) != dim(value)[2]) colnames(value) <- paste( "cell", 1:dim(value)[2])
			else colnames(value) <- make.names(colnames(value),unique=T)

			outcode <- .infernal_readSparseMatrix(scope@ptr,as(value,"dgCMatrix"),T)
		}
		else stop("assignment expects a matrix or sparse Matrix")
#		if (outcode == 0){
#			scope@cell.names <- return(.infernal_getParameter(scope@ptr, 10))
#			scope@gene.names <- return(.infernal_getParameter(scope@ptr, 11))
#		}
	return()
	}
	if (varname == "data") {
                if (class(value) == "matrix") {
			if (length(rownames(value)) != dim(value)[1]) rownames(value) <- paste( "gene", 1:dim(value)[1])
			else rownames(value) <- make.names(rownames(value),unique=T)
			if (length(colnames(value)) != dim(value)[2]) colnames(value) <- paste( "cell", 1:dim(value)[2])
			else colnames(value) <- make.names(colnames(value),unique=T)

			outcode <- .infernal_readDenseMatrix(scope@ptr,value, F)
                }else if (isS4(value) & is(value,'sparseMatrix') ) {
			if (length(rownames(value)) != dim(value)[1]) rownames(value) <- paste( "gene", 1:dim(value)[1])
			else rownames(value) <- make.names(rownames(value),unique=T)
			if (length(colnames(value)) != dim(value)[2]) colnames(value) <- paste( "cell", 1:dim(value)[2])
			else colnames(value) <- make.names(colnames(value),unique=T)

			outcode <- .infernal_readSparseMatrix(scope@ptr,as(value,"dgCMatrix"),F)
		}
		else stop("assignment expects a matrix or sparse Matrix")
#		if (outcode == 0){
#			scope@cell.names <- return(.infernal_getParameter(scope@ptr, 10))
#			scope@gene.names <- return(.infernal_getParameter(scope@ptr, 11))
#		}
	return()
	}

	err <- T
	if (expectedclass[1] == 'list'){
		if (length(value) > 1){
			if (class(value) == expectedclass[2]) err <-F
			else tryCatch({value <- as(value, expectedclass[2]); err <- F},error=function(cond){print(paste("Could not convert", class(value), "to", expectedclass[2]))})
		}else{
			if (class(value[[1]]) == expectedclass[2]) err <-F
		}
	}

	if (err){
		stop(switch(varname
, "cell.cluster" = "Assignment expects a list of integers"
, "gene.cluster" = "Assignment expects a list of integers"
, "cell.order" = "Assignment expects a list of integers"
, "gene.order" = "Assignment expects a list of integers"
, "cell.names" = "Assignment expects a list of strings"
, "gene.names" = "Assignment expects a list of strings"
))
	}else{
		
	outcode <- switch(varname,
"cell.cluster" = .infernal_setParameter(scope@ptr, 0,value),
"gene.cluster" = .infernal_setParameter(scope@ptr, 1,value),
"cell.order" = .infernal_setParameter(scope@ptr, 8,value),
"gene.order" = .infernal_setParameter(scope@ptr, 9,value),
"cell.names" = .infernal_setParameter(scope@ptr, 10,value),
"gene.names" = .infernal_setParameter(scope@ptr, 11,value),
"state.matrices" = .infernal_setParameter(scope@ptr, 14,value),
	print(paste(varname, "is not an attribute of an InferN0 scope"))
)
	}
	return()
}


#' @export
#setMethod("data","InferN0Scope", function(scope){return(InferN0Get(scope, "data"))})
#' @export
#setReplaceMethod("data","InferN0Scope", function(scope,value){return(InferN0Set(scope, "data",value))})

#' @export
#setMethod("raw.data","InferN0Scope", function(scope){return(InferN0Get(scope, "raw.data"))})
#' @export
#setReplaceMethod("raw.data","InferN0Scope", function(scope,value){return(InferN0Set(scope, "raw.data",value))})

#' @export
setMethod(f = "$", signature = "InferN0Scope", function(x, name) {return(InferN0Get(x,name))})

#' @export
setMethod(f = "$<-", signature = "InferN0Scope", function(x, name, value) {InferN0Set(x,name,value);return(x)})

#' @export
setMethod("names","InferN0Scope", function(x){return(.infernal_getParameterList(x@ptr))})

#' @export
# setReplaceMethod("[", signature(x="InferN0Scope"), function(x, i, j) {return(InferN0Set(x,i,j))})


#' @export
#setMethod("print","InferN0Scope", function(scope){return(print("hello"))})

#' @export
#setMethod("rownames","InferN0Scope",function(scope){return(InferN0Get(scope, "gene.names"))})
#' @export
#setReplaceMethod("rownames","InferN0Scope",function(scope,value){return(InferN0Set(scope, "gene.names",value))})
#' @export
#setMethod("colnames","InferN0Scope",function(scope){return(InferN0Get(scope, "cell.names"))})
#' @export
#setReplaceMethod("colnames","InferN0Scope",function(scope,value){return(InferN0Set(scope, "cell.names",value))})

