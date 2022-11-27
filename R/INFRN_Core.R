

#' Effectively execute 'InferN0FindCoexpressed', 'InferN0IdentifyNetwork' and 'InferN0PlotNetwork' while returning the output of 'InferN0IdentifyNetwork' with default parameters
#'
#' Sequentially execute 'InferN0FindCoexpressed', 'InferN0IdentifyNetwork', 'InferN0PlotNetwork'.
#' @param i0scp InferN0 scope
#'
#' @export
InferN0IdentifyNetworkPipeline <- function(i0scp, gene.list, cell.list=as.character(c()), cell.cluster="",method="EMM", nb.threads=4,max.cyclic.component=10,label.cex=0.65,edge.width=10, do.plot=T,vertex.size.base=8,vertex.size.incr=1.6, vertex.height=5,check.nb=10,check.min=0,color="#FF4400",extendedcolor="#FFFFFFFF",nb.output=0,individual.coexp.max=0, mincell=5){
	outp <-	InferN0FindCoexpressed(i0scp,  gene.list=gene.list,  cell.list=cell.list, cell.cluster=cell.cluster, nb.threads=nb.threads, nb.output=nb.output,individual.coexp.max=individual.coexp.max, mincell=mincell)
	out <- InferN0IdentifyNetwork(i0scp, c(colnames(outp$correlation), outp$gene.list), cell.cluster=cell.cluster,cell.list=cell.list,method=method,nb.threads=nb.threads,check.min=check.min,check.nb=check.nb)
	if (length(color) == 1) color <- rep(color, dim(outp$correlation)[1])
	else color <- color[match(colnames(outp$correlation), gene.list)]
	InferN0PlotNetwork(out,vertex.size.base=vertex.size.base,vertex.size.incr=vertex.size.incr,edge.width=edge.width,color= c(color, rep(extendedcolor, length(outp$gene.list) )))
return(out)}


#' Find Coexpressed
#'
#' Compute Partial correlation assoiated with a subset of all transcripts
#' @param i0scp InferN0 scope
#' @param gene.list List of transcript name for which a network is to be found (<50 gene expected)
#' @param cell.list Filter to listed cells
#' @param nb.output number of coexpressed transcript to find (default: same as size of 'gene.list')
#' @param cell.cluster uses matching clusterID as 'cell.list', ignored if custom cell.list is also provided
#' @param nb.threads number of threads required for compututations (>=1)
#' @param mincell minimum amount of cell with transcript pairs needed to evalutate the correlation.
#' @param do.call.findnetwork list of arguments supplied to a call to 'InferN0IdentifyNetwork' done uppon completion, which overwrites output.
#'
#' @export
InferN0FindCoexpressed <- function(i0scp,  gene.list,  cell.list=as.character(c()), cell.cluster="", nb.threads=4, nb.output=0,individual.coexp.max=0, mincell=5,do.append=F, call.findnetwork=c()){
	if (nb.threads < 1) stop("illegal nb.threads: needs >= 1 threads")
	if (nb.output == 0) nb.output = length(gene.list)
        if ((length(cell.list) == 0)&&(cell.cluster !="")){
                cell.clust <- InferN0Get(i0scp,"cell.cluster");
                cell.name <- InferN0Get(i0scp,"cell.names");
                cell.list <- cell.name[cell.clust == cell.cluster];
        }
	if (length(cell.list) == 0){
		nbcellscap <- length(InferN0Get(i0scp,"cell.names")) 
	}else nbcellscap <- length(cell.list)
	thr = rep(0,nbcellscap-2)
	for(n in 1:(nbcellscap-2)){
		tmptmp <- qt(0.999,df=n)^2;
		thr[n] <- tmptmp / (n + tmptmp)
	}

	tmpout <- .infernal_computePartialCorrelation(i0scp@ptr, list(gene_list=gene.list, cell_list= cell.list, nb_threads= nb.threads,nb_output=nb.output,individual_coexp_max=individual.coexp.max,mincell=mincell,minthr=thr))
	if (!is.null(call.findnetwork)) return(InferN0IdentifyNetwork(i0scp, c(colnames(tmpout$correlation), tmpout$gene.list), cell.cluster=cell.cluster,cell.list=cell.list,nb.threads=nb.threads,color=c(rep("#FF4400", dim(tmpout$correlation)[2]),rep("#0088FF", dim(tmpout$correlation)[2]))))
	if (do.append) tmpout$gene.list <- c(colnames(tmpout$correlation), tmpout$gene.list)
return(tmpout)}


#' Estimate Covariance from sparse observations
#'
#' Estimate Covariance from sparse observations, which needs to account for missing information that can induce singular matrices if partial correlations are considered. Missing information can be inputed, or eigen values can be altered. Finally, one can infer the probability of biological dropout by modeling observations, and then normalize observation by reporting deviation to expectation instead.
#'
#' @param i0scp InferN0 scope (or Input matrix gene as Rows, cells as Collumns)
#' @param method One for the following:
#' \enumerate{
#'   \item "Zero" computes covariance with 0 for missing values
#'   \item "Partial" computes partial-correlation matrix (can be singular!)
#'   \item "CorrectedPartial" computes partial-correlation matrix and ensures matrices are positive definite
#'   \item "EM" use inputation to find covariance
#'   \item "EMM" use inputation for conditionnal parameters to find covariance
#'   \item "Model" uses infered hidden states to normalize values (which assigns a non-trivial values to zero-counts)
#' }
#' @param gene.list List of transcript name for which a network is to be found (<50 gene expected)
#' @param cell.list Filter to listed cells
#' @param cell.meta Filter based on meta data associated to cells
#' @param cell.cluster uses matching clusterID as 'cell.list', ignored if custom cell.list is also provided
#' @param nb.threads number of threads required for compututations (>=1)
#' @param nb.crosseval.partition additionnal output which are multivariate sample covariance for (k-1)/k of the observations
#' @param seed seed for randomly assigning cells into training and test set partitions
#' @param outlier_frequency fraction of the data ignored as outliers (method in {"EM","EMM"} only)
#' @param EM.itemax number of EM iteration (method in {"EM","EMM"} only)
#' @return list(Covar = matrix, Mean = vector, Input = matrix, Training.Covar = list(matrix()), Training.Mean = list(vector()), Testset.Partition = list(partition_index))
#'
#' @export
InferN0ComputeCovariance <- function(i0scp, gene.list, method="EM", cell.list=c(), cell.meta=c(), cell.cluster="", nb.threads=4, nb.crosseval.partition =0,outlier_frequency =0.0,EM.itemax=10, seed = -1){
	if (nb.threads < 1) stop("illegal nb.threads: needs >= 1 threads")
	if (class(i0scp) == "InferN0Scope"){
		if (is.null(gene.list)) stop("Query 'gene.list' is required!")
		if ((length(cell.list) == 0)&&(cell.cluster !="")){
	      cell.clust <- InferN0Get(i0scp,"cell.cluster");
	      cell.name <- InferN0Get(i0scp,"cell.names");
	      cell.list <- cell.name[cell.clust == cell.cluster];
		}
		if (length(cell.list) == 0){
			nbcellscap <- length(InferN0Get(i0scp,"cell.names")) 
		}else nbcellscap <- length(cell.list)
	}else{
		if (class(i0scp) == "data.frame") i0scp <- as.matrix(i0scp)
		if (class(i0scp) == "matrix") i0scp <- Matrix(i0scp,sparse=T)
		if (isS4(i0scp) & is(i0scp,'sparseMatrix') ) i0scp <- as(i0scp,"dgCMatrix")
		if (class(i0scp) != "dgCMatrix") stop(paste(class(i0scp),"is not a valid input class; data.frame,matrix are allowed"))
		if (is.null(rownames(i0scp))) rownames(i0scp) <- paste("Gene", 1:dim(i0scp)[1])
		if (is.null(colnames(i0scp))) colnames(i0scp) <- paste("Cell", 1:dim(i0scp)[2])
		gene.list = rownames(i0scp)
		cell.list = colnames(i0scp)
	}
	if (method == "EM")	{
		listarg <- list(gene_list=gene.list, cell_list= cell.list, nb_threads=nb.threads, nbpartition=nb.crosseval.partition,nbiteration=EM.itemax, outlier_frequency=outlier_frequency, include.uncertainty=0,seed=seed)
		if (class(i0scp) == "InferN0Scope") output <- .infernal_cmpVariance(i0scp@ptr, listarg)
		else output <- .infernal_cmpVariance(i0scp, listarg)
		output$method <- method
		return(output)
	}else if (method =="EMM"){
		listarg <- list(gene_list=gene.list, cell_list= cell.list, nb_threads=nb.threads, nbpartition=nb.crosseval.partition,nbiteration=EM.itemax, outlier_frequency=outlier_frequency,include.uncertainty=1,seed=seed)
		if (class(i0scp) == "InferN0Scope") output <- .infernal_cmpVariance(i0scp@ptr, listarg)
		else output <- .infernal_cmpVariance(i0scp, listarg)
		output$method <- method
		return(output)
	}else if (method == "Model"){
		if (class(i0scp) != "InferN0Scope") stop("InferN0 scope with infered hidden states and parameter is required for this method")
		if (!"cell.state" %in% names(i0scp)) stop("InferN0 scope with infered hidden states and parameter is required for this method")
		if (is.null(cell.list)) cell.list = ""
		output <- .infernal_cmpModeledVariance(i0scp@ptr,list(gene_list=gene.list, cell_list= cell.list, nb_threads=nb.threads, nbpartition=nb.crosseval.partition,seed=seed,method=0))
		if (is.null(output)) stop("Failed")
		else{
			output$method= method
			output$Input <- Matrix(output$Input, sparse=T)
		}
		return(output)
	}else if (method =="Partial"){
		if (class(i0scp) == "InferN0Scope") output <- .infernal_computeCovar(i0scp@ptr, list(gene_list=gene.list, cell_list= cell.list, nbcross=nb.crosseval.partition,is.partial=1,seed=seed))
		else output <- .infernal_computeCovar(i0scp, list(gene_list=rownames(i0scp), cell_list= colnames(i0scp), nbcross=nb.crosseval.partition,is.partial=1,seed=seed))
		output$method <- method
		return(output)
	}else if (method =="CorrectedPartial"){
		if (class(i0scp) == "InferN0Scope") output <- .infernal_computeCovar(i0scp@ptr, list(gene_list=gene.list, cell_list= cell.list, nbcross=nb.crosseval.partition,is.partial=2,seed=seed))
		else output <- .infernal_computeCovar(i0scp, list(gene_list=rownames(i0scp), cell_list= colnames(i0scp), nbcross=nb.crosseval.partition,is.partial=2,seed=seed))
		output$method <- method
		return(output)
	}else if (method =="Zero"){
		if (class(i0scp) == "InferN0Scope") output <- .infernal_computeCovar(i0scp@ptr, list(gene_list=gene.list, cell_list= cell.list, nbcross=nb.crosseval.partition,is.partial=0,seed=seed))
		else output <- .infernal_computeCovar(i0scp, list(gene_list=rownames(i0scp), cell_list= colnames(i0scp),nbcross=nb.crosseval.partition,is.partial=0,seed=seed))
		output$method <- method
		return(output)
	}else stop("unknown method")
return(tmpout)}


#' Execute Network inference
#'
#' @param i0scp InferN0 scope, a raw input matrix or a list. If it is a matrix, 'InferN0ComputeCovariance' will be called with (or an allowed alternative, a list equivalent to what is outputed by the 'InferN0ComputeCovariance' function (if used, ignores any provided arguments that shared with that function))
#' @param list of transcript name for which a network is to be found (<50 gene expected)
#' @param method string argument needed to compute covariance ("EM", "EMM","Zero", "Partial", "CorrectedPartial","Model", see the 'InferN0ComputeCovariance' function)
#' @param evaltarget optionnal 'true' Network argument to report on true/false prediction rates
#' @param max.cyclic.component maximal size for a cyclic compnent (cycles in network), search might be slower if a more larger threshold size
#' @param cell.list list of cell names, list of cell index or T/F filter vector
#' @param cell.cluster uses matching clusterID as 'cell.list', ignored if custom cell.list is also provided
#' @param nb.threads number of threads required for compututations (>=1)
#' @param check.nb number of additionnal edges considered after an optimum network is found
#' @param check.min lowerbound on number of edges for the search depth. (if =0, sets value as size of network instead)
#' @param do.manually.select.nb.edges uses 'check.nb' as forced number of edges in output graph
#'
#' @export
InferN0IdentifyNetwork <- function(i0scp, gene.list=c(""), method= "EM", evaltarget= c(),nb.crosseval.partition=8, cell.list=as.character(c()), cell.cluster="",nb.threads=4,max.cyclic.component=10,label.cex=0.75,edge.width=10, do.plot=T,vertex.size.base=8,vertex.size.incr=1.6, vertex.height=5,check.nb=10,check.min=0,color="#FF4400",do.manually.select.nb.edges=T){
	if (nb.threads < 1) nb.threads = 1;
	if ((length(cell.list) == 0)&&(cell.cluster !="")){
		cell.clust <- InferN0Get(i0scp,"cell.cluster");
		cell.name <- InferN0Get(i0scp,"cell.names");
		cell.list <- cell.name[cell.clust == cell.cluster];
	}

	if (is.null(evaltarget)) evaltarget <- matrix(1,1,1) # will be ignored

	if (class(i0scp) == "list") {
		if (4 != sum(c("Covar", "Mean", "Input", "method") %in% names(i0scp))) stop("Expects output structure equivalent to what outputed by the 'InferN0ComputeCovariance' function")
		if (3 != sum(c("Training.Covar", "Training.Mean", "Testset.Partition") %in% names(i0scp))) stop("Expects training and set sets (obtained with nb.crosseval.partition != 0 for 'InferN0ComputeCovariance')")
		covarstr <- i0scp
	}else if (class(i0scp) == "InferN0Scope"){
		if (nb.crosseval.partition < 2) nb.crosseval.partition = 2 # argument need to be at least 2
		covarstr <- InferN0ComputeCovariance(i0scp, gene.list, cell.list=cell.list, method=method, nb.crosseval.partition=nb.crosseval.partition)
	}else stop(paste(class(i0scp), " class was used as input, InferN0Scope expected! (or list with entries matching the output of 'InferN0ComputeCovariance')"))
	flag = 0;
	if (do.manually.select.nb.edges) flag <- flag +1; 


	out <- .infernal_IdentifyNetwork(covarstr$Input, list(covarstr$Covar,covarstr$Mean,partition=covarstr$Testset.Partition, max_cyclic_component=max.cyclic.component,nb_threads=nb.threads,check_nb=check.nb,check_min=check.min, target= evaltarget,flag=flag),
		covarstr$Training.Covar,covarstr$Training.Mean);

	# verify that solution is optimum

	tryCatch({
		out$LLDerivative = solve(out$PrecisionMatrix) - covarstr$Covar
		out$LLDerivative[out$PrecisionMatrix == 0] <- 0
	},error=function(cond){print("Warning, Obtained Precision Matrix is Singular!")})

	# if the LLDerivative is a matrix filled with zeroes, than the solution is a local maximum likelihood


	if (do.plot){
		if (dim(evaltarget)[1] == dim(covarstr$Input)[1]){
			# plot ROC curve!
			plot(out$LL_incrs[,9],out$LL_incrs[,8], main = paste("ROC curve"), xlim= c(0,1), ylim= c(0,1), type="l", col='red',xlab="False Positive Rate",ylab="True Positive Rate");
			par(new=T)
			plot(out$LL_incrs[(out$nbedges+1),9],out$LL_incrs[(out$nbedges+1),8], xlim= c(0,1), ylim= c(0,1),pch=42,col="red",cex=5,xlab="",ylab="")
			.addROCCmpLines(covarstr$Covar, evaltarget)
		}else if ("EdgeListCoor" %in% names(out)){
			if (length(color) != 1) color <-color[match(colnames(out$covar), gene.list)]
			InferN0PlotNetwork(out,vertex.size.base=vertex.size.base,vertex.size.incr=vertex.size.incr,edge.width=edge.width,color=color);
	}	}
return(out)}




#' Execute Network inference Older Version...
#'
#' @param i0scp InferN0 scope (or list equivalent to what is outputed by the 'InferN0ComputeCovariance' function (if used, ignores any provided arguments that shared with that function)
#' @param list of transcript name for which a network is to be found (<50 gene expected)
#' @param method: string argument needed to compute covariance ("EM", "EMM","Zero", "Partial", "CorrectedPartial","Model" see the 'InferN0ComputeCovariance' function)
#' @param evaltarget: optionnal 'true' Network argument to report on true/false prediction rates
#' @param max.cyclic.component: maximal size for a cyclic compnent (cycles in network), search might be slower if a more larger threshold size
#' @param cell.list: list of cell names, list of cell index or T/F filter vector
#' @param cell.cluster: uses matching clusterID as 'cell.list', ignored if custom cell.list is also provided
#' @param nb.threads: number of threads required for compututations (>=1)
#' @param check.nb: upperbound on number of edges for the search depth, beyond the optimum reported. 
#' @param check.min: lowerbound on number of edges for the search depth. (if =0, sets value as size of network instead)
#' @param do.manually.select.nb.edges: uses 'check.nb' as forced number of edges in output graph
#'
#' @export
InferN0IdentifyNetworkOlderVersion <- function(i0scp, gene.list=c(""), method= "EM", evaltarget= c(),nb.crosseval.partition=8, cell.list=as.character(c()), cell.cluster="",nb.threads=4,max.cyclic.component=10,label.cex=0.75,edge.width=10, do.plot=T,vertex.size.base=8,vertex.size.incr=1.6, vertex.height=5,check.nb=10,check.min=0,color="#FF4400",do.manually.select.nb.edges=T){
	if (nb.threads < 1) nb.threads = 1;
	if ((length(cell.list) == 0)&&(cell.cluster !="")){
		cell.clust <- InferN0Get(i0scp,"cell.cluster");
		cell.name <- InferN0Get(i0scp,"cell.names");
		cell.list <- cell.name[cell.clust == cell.cluster];
	}

	if (is.null(evaltarget)) evaltarget <- matrix(1,1,1) # will be ignored

	if (class(i0scp) == "list") {
		if (4 != sum(c("Covar", "Mean", "Input", "method") %in% names(i0scp))) stop("Expects output structure equivalent to what outputed by the 'InferN0ComputeCovariance' function")
		if (3 != sum(c("Training.Covar", "Training.Mean", "Testset.Partition") %in% names(i0scp))) stop("Expects training and test sets (obtained with nb.crosseval.partition > 1 for 'InferN0ComputeCovariance')")
		covarstr <- i0scp
	}else if (class(i0scp) == "InferN0Scope"){
		if (nb.crosseval.partition < 2) nb.crosseval.partition = 2 # argument need to be at least 2
		covarstr <- InferN0ComputeCovariance(i0scp, gene.list, cell.list=cell.list, method=method, nb.crosseval.partition=nb.crosseval.partition)
	}else stop(paste(class(i0scp), " class was used as input, InferN0Scope expected! (or list with entries matching the output of 'InferN0ComputeCovariance')"))
	flag = 0;
	if (do.manually.select.nb.edges) flag <- flag +1



	out <- InferN0_identifyL0NetworkGoldResurrected(covarstr$Input, list(covarstr$Covar,covarstr$Mean,partition=covarstr$Testset.Partition, max_cyclic_component=max.cyclic.component,nb_threads=nb.threads,check_nb=check.nb,check_min=check.min, target= evaltarget,flag=flag),
		covarstr$Training.Covar,covarstr$Training.Mean);

	# verify that solution is optimum
	# LL = logdet(P) - tr(PS)
	# dLL/dP = P^{-1} - S   which should be true for sparse entries
	
	tryCatch({
		out$LLDerivative = solve(out$PrecisionMatrix) - covarstr$Covar
		out$LLDerivative[out$PrecisionMatrix == 0] <- 0
	},error=function(cond){print("Warning, Obtained Precision Matrix is Singular!")})

	print(paste("number of edges ",out$nbedges))
	# if the LLDerivative is a matrix filled with zeroes, than the solution is a local maximum likelihood

	
	if (do.plot){
		if (dim(evaltarget)[1] == dim(covarstr$Input)[1]){
			# plot ROC curve!
			plot(out$LL_incrs[,9],out$LL_incrs[,8], main = "ROC curve", xlim= c(0,1), ylim= c(0,1), type="l", col='red',xlab="False Positive Rate",ylab="True Positive Rate");
			par(new=T)
			plot(out$LL_incrs[(out$nbedges+1),9],out$LL_incrs[(out$nbedges+1),8], xlim= c(0,1), ylim= c(0,1),pch=42,col="red",cex=5,xlab="",ylab="")
			.addROCCmpLines(covarstr$Covar, evaltarget)

		}else if ("EdgeListCoor" %in% names(out)){
			if (length(color) != 1) color <-color[match(colnames(out$covar), gene.list)]
			InferN0PlotNetwork(out,vertex.size.base=vertex.size.base,vertex.size.incr=vertex.size.incr,edge.width=edge.width,color=color);
		}

	}
return(out)}

# library(InferN0); ips_inf<- InferN0LoadFile(path="/lustre/scratch117/cellgen/team218/lh20/Nguyen_Powell_ips.ifn.scp"); path = "/lustre/scratch117/cellgen/team218/vk8/scRNAseqData/Pluripotency/"; tf = as.character(unlist(read.delim(file = paste(path,"ESC_network.txt",sep="")))); inf_cov <- InferN0ComputeCovariance(i0scp = ips_inf, gene.list=tf,  method="Partial"); inf_cov <- InferN0ComputeCovariance(i0scp = ips_inf, gene.list=tf,  method="EM", nb.crosseval.partition=8, cell.list=colnames(inf_cov$Input)[colSums(inf_cov$Input) < 10], EM.itemax=3) ;inf_nw = InferN0IdentifyNetworkOlderVersion(i0scp = ips_inf, method= inf_cov,check.nb=32 )

#' Sparse Factor based Negative-Binomial modeling
#'
#' Function that normalizes the data by recovering the major coexpression factors accross cells (aka cell-types) and genes
#'
#' @param i0scp InferN0 scope
#' @param nb.hidden.row number of differentially expressed features groups
#' @param nb.hidden.col number of cell types
#' @param nb.step number of EM steps to perform
#' @param nb.thread number of threads
#' @param cell.cluster (optionnal) start search using input clusters
#' @param state.labels pair of names or name prefix used for outputs (stored within i0scp)
#' @param do.runTSNE T/F flag for runing Rtsne of hidden states
#'
#' @export
InferN0ModelData <- function(i0scp, nb.hidden.row=8,nb.hidden.col=5,nb.step=200, nb.thread=4, cell.cluster= c(), do.output.tiff=FALSE, version= 1, state.labels = c("CELL_STATE", "GENE_STATE"), do.runTSNE=T, use.heavy.prior=F, batch.map =c()){
	if (nb.hidden.row <= 1) stop("nb.hidden.row should be more than 1")
	if (nb.hidden.col <= 1) stop("nb.hidden.col should be more than 1")

	if (!is.null(cell.cluster)){
		if (length(cell.cluster) != i0scp$cell.names) stop("provided clustering does not match of its cell number")
		if (!is.factor(cell.cluster)) cell.cluster <- as.factor(cell.cluster)
		if (length(levels(cell.cluster)) > nb.hidden.col) stop("nb.hidden.col should be greater or equal to the number of clusters provided")
		.infernal_initModelWithClustering(i0scp, list(nbrow=nb.hidden.row, nbcol=nb.hidden.col, cell.cluster@.Data))
	}


	if (version == 0) fout <- .infernal_modelhidden(i0scp@ptr, list(nbhrow=nb.hidden.row,nbhcol =nb.hidden.col,nb_step=nb.step,nb_thread=nb.thread,do_output_tiff=do.output.tiff))
	else {
		if (use.heavy.prior) flag = 1
		else flag =0
		fout <- .infernal_modelhidden_v2(i0scp@ptr, list(nbhrow=nb.hidden.row,nbhcol =nb.hidden.col,nb_step=nb.step,nb_thread=nb.thread,flag=flag, batch = batch.map))
	}
	if (is.null(fout)) stop("Error")
	
	
	gc <- i0scp$cell.nbcounts; freq <-  i0scp$cell.state %*% matrix(gc,length(gc),1);	freq = diag(as.vector(freq / sum(freq)))
	tmp <- i0scp$cell.state
	i0scp@meta.data[[paste(state.labels[1],"maximum",sep="_")]] <- 0
	for(i in 1:ncol(tmp)){i0scp@meta.data[[paste(state.labels[1],"maximum",sep="_")]][i] <- which.max(tmp[,i])}
	i0scp@meta.data[[paste(state.labels[1],"maximum",sep="_")]] <- as.factor(i0scp@meta.data[[paste(state.labels[1],"maximum",sep="_")]])
	i0scp@meta.data$nb.cell.states <- Matrix::colSums(tmp != 0)
	i0scp@cell.coor[[state.labels[1]]] <- matrix(0,length(gc),ncol(i0scp$state.matrices$R))
	rownames(i0scp@cell.coor[[state.labels[1]]]) <- scp$cell.names
	colnames(i0scp@cell.coor[[state.labels[1]]]) <- paste(state.labels[1], 1:ncol(i0scp$state.matrices$R), sep="_")

	
	if (do.runTSNE){
		if (!requireNamespace("Rtsne", quietly = TRUE)) {
        		warning("The Rtsne package must be installed to get TSNEs")
	        	do.runTSNE <- F
    		}else{
		library(Rtsne);
		which <- (i0scp@meta.data$nb.cell.states > 1)
		which1 <- (i0scp@meta.data$nb.cell.states == 1)
		mmm<- Rtsne(t(as.matrix(rbind(i0scp$state.matrices$M %*% freq,i0scp$state.matrices$R %*% freq) %*% i0scp$cell.state))[which,] , check_duplicates=F)
		ncoor <- matrix(0,length(gc),2)
		ncoor[which,] <- mmm$Y
		rownames(ncoor) <- scp$cell.names
		for(i in 1:ncol(i0scp$state.matrices$R)) {
		   center <- (t(mmm$Y) %*% (i0scp$cell.state[i,which]^2))  / sum(i0scp$cell.state[i,which]^2)
		   ncoor[(i0scp$cell.state[i,] != 0.0)&(which1), ] <- rep(center, each=sum((i0scp$cell.state[i,] != 0.0)&(which1)))
		   i0scp@cell.coor[[state.labels[1]]][,i] <- i0scp$cell.state[i,]
		}
		i0scp@cell.coor$TSNE <- ncoor;
		}
	}


	gc <- i0scp$gene.nbcounts;freq <-  i0scp$gene.state %*% matrix(gc,length(gc),1);	freq = diag(as.vector(freq / sum(freq)))
	tmp <- i0scp$gene.state 
	i0scp@gene.meta.data[[paste(state.labels[2],"maximum",sep="_")]] <- 0
	for(i in 1:ncol(tmp)){i0scp@gene.meta.data[[paste(state.labels[2],"maximum",sep="_")]][i] <- which.max(tmp[,i])}
	i0scp@gene.meta.data[[paste(state.labels[2],"maximum",sep="_")]] <- as.factor(i0scp@gene.meta.data[[paste(state.labels[2],"maximum",sep="_")]])
	i0scp@gene.meta.data$nb.gene.states <- Matrix::colSums(tmp != 0)
	i0scp@gene.coor[[state.labels[2]]] <- matrix(0,length(gc),nrow(i0scp$state.matrices$R))
	rownames(i0scp@gene.coor[[state.labels[2]]]) <- scp$gene.names
	colnames(i0scp@gene.coor[[state.labels[2]]]) <- paste(state.labels[2], 1:nrow(i0scp$state.matrices$R), sep="_")

	if (do.runTSNE){	
		which <- (i0scp@gene.meta.data$nb.gene.states > 1)
		which1 <- (i0scp@gene.meta.data$nb.gene.states == 1)
		mmm<- Rtsne(t(as.matrix(t(cbind(freq %*% i0scp$state.matrices$M , freq %*% i0scp$state.matrices$R)) %*% i0scp$gene.state))[which,], check_duplicates=F)
		ncoor <- matrix(0,length(gc),2)
		ncoor[which,] <- mmm$Y
		rownames(ncoor) <- scp$gene.names
		for(i in 1:nrow(i0scp$state.matrices$R)) {
		   center <- (t(mmm$Y) %*% (i0scp$gene.state[i,which]^2))  / sum(i0scp$gene.state[i,which]^2)
		   ncoor[(i0scp$gene.state[i,] != 0.0)&(which1), ] <- rep(center, each=sum((i0scp$gene.state[i,] != 0.0)&(which1)))
		   i0scp@gene.coor[[state.labels[2]]][,i] <- i0scp$gene.state[i,]
	       	}
		i0scp@gene.coor$TSNE <- ncoor;
	}
return(i0scp)}


