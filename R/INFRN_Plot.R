
printf <- function(...) invisible(print(sprintf(...)))



#' Plot Inferred Network
#'
#' @param networkoutput: list generated by 'InferN0IdentifyNetwork'
#'
#' @export
InferN0PlotNetwork <- function(out, layout= c(),nb.edges=0, vertex.size.base=12,vertex.size.incr=2.4,edge.width=10,edge.color="#000000FF",edge.color.low="#222222FF", vertex.height=5,label.cex=0.75, incl_isolates=F,color = "#FF4400", do.make.interactive=T){
	if (!requireNamespace("igraph", quietly = TRUE)) {
	        warning("The igraph package must be installed to use this functionality")
	        return(NULL)
    	}
	library(igraph);
	if (!"PrecisionMatrix" %in% names(out)) stop("expects the output from 'InferN0IdentifyNetwork' as input")
	if (nb.edges == 0) {nb.edges = out$nbedges
	}else if (nb.edges > dim(out$EdgeListCoor)[1]) nb.edges = dim(out$EdgeListCoor)[1]

	if (nb.edges == 0) return()
	range <- 3:(nb.edges*2+2);
      lone = setdiff(1:dim(out$PrecisionMatrix)[1], (unique(as.vector(t(out$EdgeListCoor))[range])+1));
	if (incl_isolates){
		
		gr <- graph(edges=rownames(out$PrecisionMatrix)[ as.vector(t(out$EdgeListCoor + 1))[range] ], isolates = rownames(out$PrecisionMatrix)[lone], directed = F); 
	}else{
		gr <- graph(edges=rownames(out$PrecisionMatrix)[ as.vector(t(out$EdgeListCoor + 1))[range] ], directed = F); 
	}
	#gr <- permute(gr, match(1:length(names(V(gr))),order(match(names(V(gr)), colnames(out$PrecisionMatrix)))))
	#print(names(V(gr)))	
	nbnode = length(V(gr))
	if (is.null(layout)) layout <- layout_nicely(gr)
	else{
		map <- match(names(V(gr)),rownames(layout))
		layout <- layout[map,]
	}

	vmap <- match(names(V(gr)) , colnames(out$PrecisionMatrix))
	if (length(color) != 1){
		tmpcol <- rep("#AAAAFF",nbnode)
		tmpcol[!is.na(vmap)] <- color[vmap[!is.na(vmap)]]
		color <- tmpcol
	}

	xrange <- c( min(layout[,1]), max(layout[,1]))
	yrange <- c( min(layout[,2]), max(layout[,2]))

	xfact <- xrange[2] - xrange[1];
	yfact <- yrange[2] - yrange[1];

	diff <- out$LL_incrs[2:(nb.edges+1),1] - out$LL_incrs[1:(nb.edges),1];
	diff <- ((diff - min(diff)) / (max(diff) - min(diff)))

	par(mar=c(1,1,1,1) + 0.1)
	V(gr)$label.cex <- 0.0001
	col <- rep("", nb.edges)
	wid <- rep(0, nb.edges)
	for(i in 1:nb.edges){
		tmp <- -out$PrecisionMatrix[out$EdgeListCoor[i+1,1]+1,out$EdgeListCoor[i+1,2]+1] / sqrt(out$PrecisionMatrix[out$EdgeListCoor[i+1,1]+1,out$EdgeListCoor[i+1,1]+1]* out$PrecisionMatrix[out$EdgeListCoor[i+1,2]+1,out$EdgeListCoor[i+1,2]+1])
		if (tmp >0) col[i] <- rgb(tmp/2, 0.5, tmp/2)
		else col[i] <- rgb(1+tmp*0.25,-0.5 *tmp, -0.5 *tmp)
		wid <- edge.width * (1.0 + abs(tmp))/2
	}
	E(gr)$color <- col
	E(gr)$width <- wid

	print(length(E(gr)$color))
		logv <- log(diag(out$PrecisionMatrix)[vmap])
	logvrange <- range(logv)
	logv <- (logv - logvrange[1]) / (logvrange[2] - logvrange[1])
	V(gr)$color <- sapply(logv, function(x) rgb(0.4,1.0-x,1))
	print(V(gr)$color)
	print(E(gr)$color)
	V(gr)$size <- (nchar(names(V(gr)) )*vertex.size.incr + vertex.size.base)
	if (do.make.interactive){
		return(list(tk=tkplot(gr), nam= names(V(gr))))
	}else{
		print(col)
		plot.igraph(gr,layout=layout,rescale=F,asp=0,xlim=xrange,ylim=yrange,vertex.size= 0.1);
		par(new=T)
		V(gr)$label.cex <- label.cex
		bpsi <- (nchar(names(V(gr)) )*vertex.size.incr + vertex.size.base) * 0.66 * xfact
		plot.igraph(gr,layout=layout,rescale=F,asp=0,xlim=xrange,ylim=yrange, vertex.size= bpsi,
			vertex.shape="rectangle" ,vertex.size2=rep(vertex.height * yfact,nbnode),edge.width=0.001, label.font=4);
	}
}

#' Plot a consensus from a list of inferred networks
#'
#' @param networkoutput: list of lists generated by 'InferN0IdentifyNetwork'
#'
#' @export
InferN0PlotNetworkConsensus <- function(glist, layout= c(), vertex.size.base=12,vertex.size.incr=2.4,edge.width=10,edge.color="#000000FF",edge.color.low="#222222FF", vertex.height=5,label.cex=0.75, incl_isolates=F,color = "#FF4400", do.make.interactive=T){
	if (!requireNamespace("igraph", quietly = TRUE)) {
	        warning("The igraph package must be installed to use this functionality")
	        return(NULL)
    	}
	library(igraph);
	if ((class(glist) != "list")|| (!"PrecisionMatrix" %in% names(glist[[1]]))) stop("expects a list of output from 'InferN0IdentifyNetwork' as input")

	vnam <- colnames(glist[[1]]$PrecisionMatrix)

	edges <- glist[[1]]$EdgeListCoor[2:(glist[[1]]$nbedges+1),]
	edges <- cbind(edges, rep(1,nrow(edges)), 1+edges[,1] + (edges[,2] * length(vnam)))
	edges <- cbind(edges, glist[[1]]$PrecisionMatrix[edges[,4]])
	for( i in 2:length(glist)){
		print(dim(edges))
		eval <- 1 + glist[[i]]$EdgeListCoor[2:(glist[[i]]$nbedges+1),1] + (glist[[i]]$EdgeListCoor[2:(glist[[i]]$nbedges+1),2] * length(vnam));
		map <- match(eval, edges[,4])
		edges[map[!is.na(map)],3] <- edges[map[!is.na(map)],3] +1
		edges[map[!is.na(map)],5] <- edges[map[!is.na(map)],5] + glist[[i]]$PrecisionMatrix[edges[map[!is.na(map)],4]]

		edges <- rbind(edges , cbind( glist[[i]]$EdgeListCoor[2:(glist[[i]]$nbedges+1),][is.na(map),] , rep(1, sum(is.na(map))), eval[is.na(map)], glist[[i]]$PrecisionMatrix[eval[is.na(map)]] ))
	}
	edges <- cbind(edges[,2], edges[,1],edges[,3:5])
	edges <- edges[order(edges[,3]),]
	lone = setdiff(1:length(vnam), (unique(as.vector(edges[,1:2])+1)));
	if (incl_isolates){
		gr <- graph(edges=vnam[ as.vector(t(edges[,1:2]) + 1)], isolates = vnam[lone], directed = F); 
	}else{
		gr <- graph(edges=vnam[ as.vector(t(edges[,1:2]) + 1)], directed = F); 
	}
	#gr <- permute(gr, match(1:length(names(V(gr))),order(match(names(V(gr)), colnames(out$PrecisionMatrix)))))
	print(names(V(gr)))
	gr <- set_edge_attr(gr,"weight", value = edges[,3])	
	nbnode = length(V(gr))
	if (is.null(layout)) layout <- layout_nicely(gr)
	else{
		if (!is.null(rownames(layout))){
			map <- match(names(V(gr)),rownames(layout))
			layout <- layout[map,]
		}
	}

	xrange <- c( min(layout[,1]), max(layout[,1]))
	yrange <- c( min(layout[,2]), max(layout[,2]))
	xfact <- xrange[2] - xrange[1];
	yfact <- yrange[2] - yrange[1];

	par(mar=c(1,1,1,1) + 0.1)
	V(gr)$label.cex <- 0.0001
	V(gr)$color <- sapply(names(V(gr)), function(x){return(ifelse(x %in% c("MEG3", "FCN1", "DYNLT3", "MAGEH1", "EMC4"), "#00FFFF", "#FF00FF"))})

	E(gr)$color <- mapply(function(x,y){if (x != length(glist)) return(ifelse(y<0,"#00DD00", "#FF0000")); return(ifelse(y<0,"#006600", "#880000"));}, edges[,3] ,edges[,5])
	E(gr)$width <- edge.width * (edges[,3] / length(glist))
	
#	logv <- log(diag(out$PrecisionMatrix)[vmap])
#	logvrange <- range(logv)
#	logv <- (logv - logvrange[1]) / (logvrange[2] - logvrange[1])
#	V(gr)$color <- sapply(logv, function(x) rgb(0.4,1.0-x,1))
#	print(V(gr)$color)
#	print(E(gr)$color)
	V(gr)$size <- (nchar(names(V(gr)) )*vertex.size.incr + vertex.size.base)
	if (do.make.interactive){
		tk <- tkplot(gr)
		coor <- tk_coords(tk)
		rownames(coor) <- names(V(gr))
		return(list(tk=tk, nam= names(V(gr)), coor=coor))
	}else{
		plot.igraph(gr,layout=layout,rescale=F,asp=0,xlim=xrange,ylim=yrange,vertex.size= 0.1);
		par(new=T)
		V(gr)$label.cex <- label.cex
		bpsi <- (nchar(names(V(gr)) )*vertex.size.incr + vertex.size.base) * 0.66 * xfact
		plot.igraph(gr,layout=layout,rescale=F,asp=0,xlim=xrange,ylim=yrange, vertex.size= bpsi,
			vertex.shape="rectangle" ,vertex.size2=rep(vertex.height * yfact,nbnode),edge.width=0.001, label.font=4);
	}

}

.addROCCmpLines <- function(Covar, target){
	par(new=T)
	m <- matrix(1,2,2)
  m[1,] <- m[1,] * 0
	plot(m, type="l", xlim= c(0,1), ylim= c(0,1), col= "gray",xlab="",ylab="")
	n <- dim(Covar)[1]
	daord <- order(abs(Covar),decreasing=T);
	damat <- matrix(0,1 + (n * (n -1))/2, 2)
	damat[1,1] =0; damat[1,2] =0;
	j=2;
	newrow <- c(0,0)
	for(i in daord){
		if (((i-1) %% n) < floor((i-1) / n)){
			if (target[(i %% n), (floor((i-1) / n)+1)] != 0) newrow[2] <- newrow[2] +1
			else newrow[1] <- newrow[1] +1
			damat[j,] <- newrow
			j <- j + 1
		}
	}
	damat[,1] <- damat[,1] / damat[(j-1),1]
	damat[,2] <- damat[,2] / damat[(j-1),2]
	par(new=T)
	plot(damat, type="l", xlim= c(0,1), ylim= c(0,1), col= "blue",xlab="",ylab="")

	tryCatch({
		library('QUIC')
		tryCatch({
			quicout <- QUIC(out$Covar, 0 , msg=0)
			daord <- order(abs(quicout$X),decreasing=T);
			damat <- matrix(0,1 + (n * (n -1))/2, 2)
			damat[1,1] =0; damat[1,2] =0;
			j=2;
			newrow <- c(0,0)
			for(i in daord){
				if (((i-1) %% n) < floor((i-1) / n)){
					if (target[(i %% n), (floor((i-1) / n)+1)] != 0) newrow[2] <- newrow[2] +1
					else newrow[1] <- newrow[1] +1
					damat[j,] <- newrow
					j <- j + 1
				}
			}
			damat[,1] <- damat[,1] / damat[(j-1),1]
			damat[,2] <- damat[,2] / damat[(j-1),2]
			par(new=T)
			plot(damat, type="l", xlim= c(0,1), ylim= c(0,1), col= "green",xlab="",ylab="")
		},error=function(cond){print("Warning, QUIC package failed to produre l1 constrained Precision matrix!")})
	},error=function(cond){print("Warning, QUIC package is not installed for some benchmarking!")})
}

#' Function used to easily populate ggplot related options stored in a list, that can be used by ploting function "overlay" as plot.attrib argument
#'
#' @param infrnscp: InferN0 scope
#' @param query: meta.data collumn or gene
#' @param transform: string representing a transformation of input field ("Z" Zscore, "N", q-value, "E" log2(expectation+1) under model, "V" log2(variance) under model) 
#' @param coor: 1 to 3 strings representing x,y,z coordinates for overlay. If a single string is provided, a pair with "_0" and "_1" is used instead.
#' @param filter: T/F vector to select cells or genes to render
#' @param filter.meta: meta.data/gene.meta filter criterions to select cells/genes to render
#' @param plot.attribs: additionnal flags to change ggplot layout, use the output from InferN0plotAttribs
#' @param color.palette: list of colors for overlays if the query is a numerical value
#' @param pt.size: size of drawn points
#' @return ggplot2 object or rgl overlay
#' @examples
#' scp <- InferN0Init(datamatrix)
#' InferN0Overlay(scp, "geneA", coor= c("geneB", "geneC"))
#' @export
InferN0Overlay <- function(infrnscp, query, transform=c(), coor= "TSNE", filter=c(), filter.meta=c(),plot.attribs=c(), color.palette="fireice",na.alpha = 0.1, na.color ="#888888", pt.size=1){
	use.shape <- c()
	use.alpha <- c()

	if (length(query) != 1){
		catequery <- query		
		query <- query[1]
	}else catequery <- c()


	if (query %in% names(infrnscp@meta.data)) {
		mode <- "cells"
		filter <- .celllistquery(infrnscp, filter.meta, filter)
		if (sum(filter) == 0) {warning("Current cell filter filters every cell")}
		data <- infrnscp@meta.data[filter, query]
		if (class(data) == "characters") data <- as.factor(data)
		if (class(data) == "factor"){
			if (query %in% names(infrnscp@display)) {
				if (class(infrnscp@display[[query]]) == "list") color <- infrnscp@display[[query]]$color
				else color <- infrnscp@display[[query]]
			}else{
				color <- .myrainbow(length(levels(data)))
				names(color) <- levels(data)
			}
		}else color <- c()
	}else if (query %in% infrnscp$gene.names) {
		mode <- "cells"
		filter <- .celllistquery(infrnscp, filter.meta, filter)
		if (sum(filter) == 0) {warning("Current cell filter filters every cell")}
		if (is.null(transform)){

			transform = "log2 Counts +1"
			data <- log2(.infernal_getDataRow(infrnscp@ptr, list(G=query,C="", 0))[filter]+1)
			data[ data == 0] <- NA
		}
		else if (transform == "N" ) {data <- pnorm(as.vector(.infernal_getDataRow(infrnscp@ptr, list(G=query,C="", 1)))[filter]); transform = "Erf(Zscore)"}
		else if (transform == "n" ) data <- pnorm(as.vector(.infernal_getDataRow(infrnscp@ptr, list(G=query,C="", 6)))[filter])
		else if (transform == "E" ) {data <- log2(as.vector(.infernal_getDataRow(infrnscp@ptr, list(G=query,C="", 2)))[filter]+1); transform = "Log2(Expectation+1)"}
		else if (transform == "V" ) {data <- log2(as.vector(.infernal_getDataRow(infrnscp@ptr, list(G=query,C="", 3)))[filter]); transform = "Log2 Variance"}
		else if (transform == "R" ) data <- as.vector(.infernal_getDataRow(infrnscp@ptr, list(G=query,C="", 4)))[filter]
		else if (transform == "P" ) data <- as.vector(.infernal_getDataRow(infrnscp@ptr, list(G=query,C="", 5)))[filter]
#  infrnscp$data[query, filter]
		else {
			transform = "log2 Counts +1"
			data <- log2(.infernal_getDataRow(infrnscp@ptr, list(G=query,C="", 0))[filter]+1)
			data[ data == 0] <- NA
		}
		#infrnscp$raw.data[query, filter]
		
		color <- c()
	}else if (query %in% infrnscp$cell.names) {
		mode <- "genes"
		filter <- .genelistquery(infrnscp, filter.meta, filter)
		if (is.null(transform)){
			transform = "log2 Counts +1"
			data <- log2(.infernal_getDataRow(infrnscp@ptr, list(G=query,C="", 0))[filter]+1)
			data[ data == 0] <- NA
		}
		else if (transform == "N" ) {data <- pnorm(as.vector(.infernal_getDataRow(infrnscp@ptr, list(G="",C=query, 1)))[filter]); transform = "Erf(Zscore)"}
		else if (transform == "n" ) data <- pnorm(as.vector(.infernal_getDataRow(infrnscp@ptr, list(G="",C=query, 6)))[filter])
		else if (transform == "E" ) {data <- log2(as.vector(.infernal_getDataRow(infrnscp@ptr, list(G="",C=query, 2)))[filter] + 1); transform = "Log2(Expectation+1)"}
		else if (transform == "V" ) {data <- log2(as.vector(.infernal_getDataRow(infrnscp@ptr, list(G="",C=query, 3)))[filter] ); transform = "Log2 Variance"}
		else if (transform == "R" ) data <- as.vector(.infernal_getDataRow(infrnscp@ptr, list(G="",C=query, 4)))[filter]
		else if (transform == "P" ) data <- as.vector(.infernal_getDataRow(infrnscp@ptr, list(G="",C=query, 5)))[filter]
		#<- infrnscp$data[filter,query]
		else{
			transform = "log2 Counts +1"
			data <- log2(.infernal_getDataRow(infrnscp@ptr, list(G="",C=query, 0))[filter]+1)
			data[ data == 0] <- NA
		}
		# data <- infrnscp$raw.data[, query]
		color <- c()
	}else if (query %in% names(infrnscp@gene.meta.data)) {
		mode <- "genes"
		filter <- .genelistquery(infrnscp, filter.meta, filter)
		data <- infrnscp@gene.meta.data[filter, query]
		if (class(data) == "characters") data <- as.factor(data)
		if (class(data) == "factor"){
			if (query %in% names(infrnscp@colors)) color <- infrnscp@colors[[query]]
			else{
				color <- .myrainbow(length(level(data)))
				names(color) <- levels(data)
			}
		}else color <- c()
	}else if (query %in% names(infrnscp@gene.coor)) {
		mode <- "genes"
		filter <- .genelistquery(infrnscp, filter.meta, filter)
		data <- infrnscp@gene.coor[filter, query]
		color <- c()
	}else if (query %in% names(infrnscp@cell.coor)) {
		mode <- "cells"
		filter <- .celllistquery(infrnscp, filter.meta, filter)
		if (sum(filter) == 0) {warning("Current cell filter filters every cell")}
		data <- infrnscp@cell.coor[filter, query]
		color <- c()
	}else{
		warning("Query does not match any metadata or gene")
		return();
	}

	if (!(is.null(catequery))){
		data <- grepl(catequery[2], data);
	}



	if (mode == "cells"){
		if (length(coor) == 1){
			if (!coor %in%names(infrnscp@cell.coor)) stop(paste(coor, "is not a cell coordinate"))
			if (ncol(infrnscp@cell.coor[[coor]]) >= 3) {
				gdata <- data.frame(infrnscp@cell.coor[[coor]][filter,1:3])
				coor <- paste(coor, 1:3,sep = "_")
			}
			else {
				gdata <- data.frame(infrnscp@cell.coor[[coor]][filter,1:2])
				coor <- paste(coor, 1:2,sep = "_")
			}
		}else{
			for(i in 1:min(3,length(coor))){
				coorname <- sub("_[^_]", "", coor[i])
				col <- as.integer(sub(".*_", "", coor[i]))
				if (!coorname %in% names(infrnscp@cell.coor)){
					stop(paste("did not find coordinate",coorname,"data for cell"))
					return;
				}
				if ((col == 0)||(col > ncol(infrnscp@cell.coor[[coorname]]))){
					stop(paste("invalic coordinate column selected: ",col))
					return;
				}
				if (i == 1) {gdata <- data.frame(infrnscp@cell.coor[[coorname]][filter, col]); colnames(gdata) <- coor[1];}
				else gdata[[coor[i]]] <- infrnscp@cell.coor[[coorname]][filter, col] 
			}
		}
	}else{
		if (length(coor) == 1){
			if (!coor %in%names(infrnscp@gene.coor)) stop(paste(coor, "is not a gene coordinate"))
			if (ncol(infrnscp@gene.coor[[coor]]) >= 3) {
				gdata <- data.frame(infrnscp@gene.coor[[coor]][filter,1:3])
				coor <- paste(coor, 1:3,sep = "_")
			}
			else {
				gdata <- data.frame(infrnscp@gene.coor[[coor]][filter,1:2])
				coor <- paste(coor, 1:2,sep = "_")
			}
		}else{
			for(i in 1:min(3,length(coor))){
				if (!coor[i] %in% names(infrnscp@gene.coor[[coorname]])){
					warning("did not find coordinate data for gene")
					return;
				}
			}
			gdata <- data.frame(infrnscp@gene.coor[filter, coor])
		}
	}


	if (is.null(color)){
		# is numeric since no color is assigned yet
		if (!is.null(transform)){
		if (transform == "log2p1") data <- log2(data + 1)
		else if (transform == "log2") {data <- log2(data); plotdata[is.infinite(data)] <- NA}
		else if (transform == "log10") {data <- log10(data); plotdata[is.infinite(data)] <- NA}
		else if (transform == "log") {data <- log(data); plotdata[is.infinite(data)] <- NA}
		} else transform = "value"
	}else{


	}

	if (length(coor) == 2){
		
		colnames(gdata) <- c("X","Y")
		gdata[["A"]] <- data
		if (is.null(color)){
			p <- ggplot(gdata, aes(x=X,y=Y,color=A, alpha=A))
			daccrange <- range(gdata[["A"]],na.rm=T)
			if (is.infinite(daccrange[1])) daccrange = 1:41
			else{
				if (daccrange[1] * daccrange[2] < 0) {
					if (-daccrange[1] > daccrange[2]) daccrange = 1:(21+ floor(-20 *daccrange[2] /daccrange[1]))
					else daccrange = (21- floor(-20 *daccrange[1] /daccrange[2])):41
				}else daccrange = 1:41
			}
			if ((length(color.palette) != 1)||(color.palette=="fireice")) colpal <- colorRampPalette(c("#00FFFF", "#00B0FF","#0079FF","#0000E8", "#000074","#000000","#4B0000","#960000","#E10000","#FF8000","#FFD600"))(41)[daccrange]
			else if (color.palette == "fire") colpal <- colorRampPalette(c("#000000","#4B0000","#960000","#E10000","#FF8000","#FFD600"))(41)[daccrange]
	 	       	else if (color.palette == "ice") colpal <- colorRampPalette(rev(c("#00FFFF", "#00B0FF","#0079FF","#0000E8", "#000074","#000000")))(41)[daccrange]
			else colpal <- colorRampPalette(color.palette)(41)[daccrange]
	
			p <- p + scale_color_gradientn(name=transform,colours=colpal, na.value= na.color)
			p <- p + scale_alpha_continuous(position=NULL,guide="none", na.value=na.alpha, range = c(1, 1))
		}else{
			p <- ggplot(gdata, aes(x=X,y=Y,color=A))
			p <- p + scale_color_manual(name=query,values=color,label=names(color),drop = FALSE)
			if (!is.null(use.shape)) p <- p + scale_shape_manual(name=query,values=use.shape,label=names(use.shape),drop = FALSE)
			if (!is.null(use.alpha)) p <- p + scale_alpha_manual(name=query,values=use.alpha,label=names(use.alpha),drop = FALSE)
			p <- p + guides(colour = guide_legend(override.aes = list(shape = 15,size = 5)))
		}
		p <- p + xlab(coor[1]) + ylab(coor[2])
	}else{
		if (!requireNamespace("rgl", quietly = TRUE)) {
        		warning("The rgl package must be installed to render 3D plots")
	        	return(NULL)
    		}
		library(rgl)
		evalcolor = rep("#000000", sum(filter))
		if (!is.null(color)){
			evalcolor <- color[data]
		}else{
			aur <- range(data,na.rm=T)
			color <- (data - aur[1]) / (aur[2] - aur[1])
			flt = !is.na(color)
			if ((length(color.palette) != 1)||(color.palette=="fireice")) evalcolor <- colorRampPalette(c("#00FFFF", "#00B0FF","#0079FF","#0000E8", "#000074","#000000","#4B0000","#960000","#E10000","#FF8000","#FFD600"))(41)[color[filter]*40]
			else if (color.palette == "fire") evalcolor <- colorRampPalette(c("#000000","#4B0000","#960000","#E10000","#FF8000","#FFD600"))(41)[color[filter]*40]
	 	       	else if (color.palette == "ice") evalcolor <- colorRampPalette(rev(c("#00FFFF", "#00B0FF","#0079FF","#0000E8", "#000074","#000000")))(41)[color[filter]*40]
			else evalcolor <- colorRampPalette(color.palette)(41)[color[filter]*40]
			evalcolor[!flt] = "#AAAAAA"
		}
		plot3d(x = gdata[,1], y= gdata[,2], z= gdata[,3] , col=evalcolor)
		return;
	}
return(.changeStyle(p+ geom_point(size=pt.size),plot.attribs))}

#' Plot for proportion/frequency for pair for meta data entries
#'
#' @param ggplot.or.plot.attrib: Either a ggplot object or a default list of attributes that is to be modified (effectively equivalent to setting the "baseAttribs" instead). If a ggplot is provided, outputs an ggplot object with the applied modifications instead.
#' @param baseAttribs: list of attributes that is to be modified
#'
#' @export
InferN0FrequencyBarPlot <- function(xannot,hannot, wannot=c(), do.horizscale=T, use.colors=c(), show.label.lowerbound.fraction= 0.05, do.resize=c(), is.vertical=F, plot.attribs=c(),return.counts.instead=F,make.piechart=F){
	if ("flags" %in% names(plot.attribs)) flags.plot <- plot.attribs[["flags"]]
	else flags.plot <- c()

	if (!is.null(do.resize)) setPlotSize(do.resize[1],do.resize[2])
	if (length(xannot) != length(hannot)) {
		print(table(xannot))
		print(table(hannot))
		 stop("annotation length mistaches")
	}
	if (class(xannot) == "factor") {
		vlabs = levels(xannot)
		if ("y.flt.empty" %in% flags.plot){
			tmp <- table(xannot)
			print(tmp)
			print(tmp[vlabs])
			vlabs <- vlabs[tmp[vlabs] > 0] # remove 0
		}
	}else vlabs = unique(xannot)
	if (is.vertical == FALSE) vlabs <- rev(vlabs)

	if (class(hannot) == "factor") {
		hlabs = levels(hannot)
		if ("x.flt.empty" %in% flags.plot){
			tmp <- table(hannot)
			hlabs <- hlabs[tmp[hlabs] > 0] # remove 0
		}
	}else hlabs = unique(hannot)

	if ("x.remap.class" %in% names(plot.attribs)){
		remap <- match(levels(xannot), names(plot.attribs$x.remap.class))
		filter <- !is.na(remap)
		levels(xannot)[filter] <- as.character(plot.attribs$x.remap.class[remap[filter]])
	}

	if ("x.rev.order" %in% flags.plot) hlabs <- rev(hlabs)
	if ("y.rev.order" %in% flags.plot) vlabs <- rev(vlabs)
	
	if (is.null(use.colors)) use.colors <- myrainbow(length(hlabs))
	counts <- matrix(0,length(hlabs), length(vlabs))


	if (length(wannot) == length(xannot)){
		for(i in 1:length(xannot)) {
			pa <- match(hannot[i], hlabs);
			pb <- match(xannot[i], vlabs);
			counts[pa, pb] = counts[ pa, pb] + wannot[i]
		}
	}else{
		for(i in 1:length(xannot)) {
			pa <- match(hannot[i], hlabs);
			pb <- match(xannot[i], vlabs);
			counts[pa, pb] = counts[ pa, pb] + 1
		}
	}



	if (return.counts.instead){
		rownames(counts) <- hlabs
		colnames(counts) <- vlabs
		return(counts)
	}

	tt = structure(counts, .Dim = c(length(hlabs), length(vlabs)), .Dimnames = structure(list(Cluster = hlabs, verti = vlabs), .Names = c("Cluster", "verti")), class = "table")
	tt <- as.data.frame(tt)
	tt[["Frequency"]] <- tt[["Freq"]]
	if (do.horizscale) {
		denums <- colSums(counts)
		names(denums) <- vlabs
		for(i in 1:nrow(tt)) tt[i, "Frequency"] <- tt[i, "Frequency"] / denums[tt[i, "verti"]]
	}
	tt[["Label"]] <- as.character(tt[["Freq"]])
	for(i in 1:nrow(tt)) {
		if (is.na(tt[i,"Frequency"])||(tt[i,"Frequency"] <= show.label.lowerbound.fraction)) tt[[i,"Label"]]  <- ""
	}

	library(ggplot2)
	p <- ggplot(as.data.frame(tt),aes(x=factor(verti),y=Frequency,fill=Cluster))
	if ("no.xnamedticks" %in% flags.plot) p <-p + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
	else{
		if ("no.xticks" %in% flags.plot) p <-p + theme(axis.ticks.x=element_blank())
		if ("no.xnames" %in% flags.plot) p <-p + theme(axis.text.x=element_blank())
		else if (!(("rot.xnames" %in% flags.plot)||("bars.rot.xnames" %in% flags.plot))) p <- p + theme(axis.text.x=element_text(angle = 90, hjust = 1))
	}

	if (is.vertical){
	#if ("x.rev.order" %in% flags.plot)  p <- p + geom_bar(stat="identity",position=position_stack(reverse = TRUE))
	#else
	p <- p + geom_bar(stat="identity",position=position_stack(reverse = FALSE))
	if ("make.pie" %in% flags.plot) p <- p + coord_polar(theta="y")
	else{
	p <- p + geom_text(aes(label=Label),position=position_stack(vjust=0.5,reverse = FALSE))
	}
	
	p <- p +  scale_fill_manual(values=use.colors)
	}else{
	#if ("x.rev.order" %in% flags.plot) p <- p + geom_bar(stat="identity",position=position_stack(reverse = FALSE))
	#else
	p <- p + geom_bar(stat="identity",position=position_stack(reverse = TRUE))
	if ("make.pie" %in% flags.plot ) {
		p <- p + coord_polar(theta="y")
	 	p <- p + geom_text(aes(label=Label),position=position_stack(vjust=0.5,reverse = TRUE))
		p <- p +  scale_fill_manual(values=use.colors)
	}else{
  p <- p + geom_text(aes(label=Label),position=position_stack(vjust=0.5,reverse = TRUE))+
  scale_fill_manual(values=use.colors)+
  coord_flip()
	}
	}
	p <- p + guides(fill=guide_legend(ncol=1))
return(changeStyle(p,plot.attribs,"bar"))}


#' Function used to easily populate ggplot related options stored in a list, that can be used by ploting function "Overlay" as plot.attrib argument be applied on a input ggplot.
#'
#' @param ggplot.or.plot.attrib: Either a ggplot object or a default list of attributes that is to be modified (effectively equivalent to setting the "baseAttribs" instead). If a ggplot is provided, outputs an ggplot object with the applied modifications instead.
#' @param baseAttribs: list of attributes that is to be modified
#'
#' @export
InferN0plotAttribs <-function(..., ggplot.or.baseAttribs=c(), baseAttribs=c(),no.legend=F,no.jitter=F, no.xticks=F, no.xlabel=F, no.yticks=F, no.ylabel=F, title=c(), xtitle= c(), ytitle=c(),xlabel=c(), ylabel=c(), xnames.bold=F,xnames.color=c(), xnames.size=c(), ynames.bold=F,ynames.color=c(), ynames.size=c(),legnames.size=c(),no.xnamedticks=F,no.ynamedticks=F,nb.col.legend=c()){
	if ((is.null(ggplot.or.baseAttribs))||(inherits(ggplot.or.baseAttribs, "ggplot"))) {
		if (is.null(baseAttribs)) fout <- list(flags=c())
		else fout <- baseAttribs
	}else fout <- ggplot.or.baseAttribs

	fout[["flags"]]  <- c(fout[["flags"]], ...)
	if (no.jitter) fout[["flags"]] <- c(fout[["flags"]], "no.jitter")
	if (no.legend) fout[["flags"]] <- c(fout[["flags"]], "no.legend") 
	if (no.xlabel) fout[["flags"]] <- c(fout[["flags"]], "no.xlabel")
	if (no.xticks) fout[["flags"]] <- c(fout[["flags"]], "no.xticks")
	if (no.xnamedticks) fout[["flags"]] <- c(fout[["flags"]], "no.xnamedticks")
	if (no.ylabel) fout[["flags"]] <- c(fout[["flags"]], "no.ylabel")
	if (no.yticks) fout[["flags"]] <- c(fout[["flags"]], "no.yticks")
	if (no.ynamedticks) fout[["flags"]] <- c(fout[["flags"]], "no.ynamedticks")

	if (xnames.bold) fout[["flags"]] <- c(fout[["flags"]], "xnames.bold")
	if (ynames.bold) fout[["flags"]] <- c(fout[["flags"]], "ynames.bold")


	if (!is.null(xnames.color)) fout[["xnames.color"]] = xnames.color
	if (!is.null(xnames.size)) fout[["xnames.size"]] = xnames.size
	if (!is.null(ynames.color)) fout[["ynames.color"]] = ynames.color
	if (!is.null(ynames.size)) fout[["ynames.size"]] = ynames.size
	if (!is.null(legnames.size)) fout[["legnames.size"]] = legnames.size

	if (!is.null(title)) fout[["title"]] = title
	if (!is.null(xlabel)) fout[["xlabel"]] = xlabel
	if (!is.null(ylabel)) fout[["ylabel"]] = ylabel
	if (!is.null(nb.col.legend)) fout[["nb.col.legend"]] = nb.col.legend


	if (inherits(ggplot.or.baseAttribs, "ggplot")) return(.changeStyle(.ggplot.or.baseAttribs, fout))
return(fout)}


