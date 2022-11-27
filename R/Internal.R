


.celllistquery <- function(infrscp, metaquery=c(), listquery=c(),default=T,getnames=F,negate=F){
	if (!is.null(listquery)){
		if (length(listquery) != nrow(infrscp@meta.data)) {
			listquery <- !is.null(match(listquery, rownames(infrscp@meta.data)))
		}else{
			listquery[is.na(listquery)] <- F
			listquery <- listquery
		}
		if (is.null(metaquery)) return(listquery)
	}else{
		if (is.null(metaquery)) return(rep(default,nrow(infrscp@meta.data)))
		listquery <- rep(T,nrow(infrscp@meta.data))
	}
	if (class(metaquery) == "list"){
		dalist <- rep(T,nrow(infrscp@meta.data))
		for(i in names(metaquery)){
			if (!i %in% colnames(infrscp@meta.data)) stop(paste("Missing", i, "in meta.data"))
			listquery = listquery & (infrscp@meta.data[[i]] %in% c(metaquery[[i]]))
		}
	}else{
	if (length(metaquery) == 1) {
		if (!metaquery %in% colnames(infrscp@meta.data)) stop(paste("Missing", metaquery, "in meta.data"))
		dalist <- infrscp@meta.data[[metaquery]]
	}else{
		if (!metaquery[1] %in% colnames(infrscp@meta.data)) {
			if (substr(metaquery[1],1,1) %in% c(">", "<")){
				if (substr(metaquery[1],2,2) == "=") {
					if (substr(metaquery[1],1,1) == "<") dalist <- infrscp@meta.data[[substr(metaquery[1],3,10000)]] <= as.numeric(metaquery[2])
					else dalist <- infrscp@meta.data[[substr(metaquery[1],3,10000)]] >= as.numeric(metaquery[2])
				}else{
					if (substr(metaquery[1],1,1) == "<") dalist <- infrscp@meta.data[[substr(metaquery[1],2,10000)]] < as.numeric(metaquery[2])
					else dalist <- infrscp@meta.data[[substr(metaquery[1],2,10000)]] > as.numeric(metaquery[2])
				}
				dalist[is.na(dalist)] <- F
			}else{
				if (substr(metaquery[1],1,1) == "!"){
					negate = T
					metaquery[1] <- substr(metaquery[1],2,10000)
				}
				if (!metaquery[1] %in% colnames(infrscp@meta.data)) stop(paste("Missing", metaquery[1], "in meta.data"))
				dalist <- !is.na(match(infrscp@meta.data[[metaquery[1]]] , metaquery[2:length(metaquery)]))
			}
		}else dalist <- !is.na(match(infrscp@meta.data[[metaquery[1]]] , metaquery[2:length(metaquery)]))
	}
        }
        if (negate) dalist <- !dalist
        if (getnames) return(rownames(infrscp@meta.data)[dalist & listquery]) else return(dalist & listquery)
}

.genelistquery <- function(infrscp, metaquery=c(), listquery=c(),default=T,getnames=F,negate=F){
	if (!is.null(listquery)){
		if (length(listquery) != nrow(infrscp@gene.meta.data)) {
			listquery <- !is.null(match(listquery, rownames(infrscp@gene.meta.data)))
		}else{
			listquery[is.na(listquery)] <- F
			listquery <- listquery
		}
		if (is.null(metaquery)) return(listquery)
	}else{
		if (is.null(metaquery)) return(rep(default,nrow(infrscp@gene.meta.data)))
		listquery <- rep(T,nrow(infrscp@gene.meta.data))
	}
	if (class(metaquery) == "list"){
		dalist <- rep(T,nrow(infrscp@gene.meta.data))
		for(i in names(metaquery)){
			if (!i %in% colnames(infrscp@gene.meta.data)) stop(paste("Missing", i, "in gene.meta.data"))
			listquery = listquery & (infrscp@gene.meta.data[[i]] %in% c(metaquery[[i]]))
		}
	}else{
	if (length(metaquery) == 1) {
		if (!metaquery %in% colnames(infrscp@gene.meta.data)) stop(paste("Missing", metaquery, "in gene.meta"))
		dalist <- infrscp@gene.meta.data[[metaquery]]
	}else{
		if (!metaquery[1] %in% colnames(infrscp@gene.meta.data)) {
			if (substr(metaquery[1],1,1) %in% c(">", "<")){
				if (substr(metaquery[1],2,2) == "=") {
					if (substr(metaquery[1],1,1) == "<") dalist <- infrscp@gene.meta.data[[substr(metaquery[1],3,10000)]] <= as.numeric(metaquery[2])
					else dalist <- infrscp@gene.meta.data[[substr(metaquery[1],3,10000)]] >= as.numeric(metaquery[2])
				}else{
					if (substr(metaquery[1],1,1) == "<") dalist <- infrscp@meta.data[[substr(metaquery[1],2,10000)]] < as.numeric(metaquery[2])
					else dalist <- infrscp@gene.meta.data[[substr(metaquery[1],2,10000)]] > as.numeric(metaquery[2])
				}
				dalist[is.na(dalist)] <- F
			}else{
				if (substr(metaquery[1],1,1) == "!"){
					negate = T
					metaquery[1] <- substr(metaquery[1],2,10000)
				}
				if (!metaquery[1] %in% colnames(infrscp@gene.meta.data)) stop(paste("Missing", metaquery[1], "in gene.meta"))
				dalist <- !is.na(match(infrscp@gene.meta.data[[metaquery[1]]] , metaquery[2:length(metaquery)]))
			}
		}else dalist <- !is.na(match(infrscp@gene.meta.data[[metaquery[1]]] , metaquery[2:length(metaquery)]))
	}
        }
        if (negate) dalist <- !dalist
        if (getnames) return(rownames(infrscp@gene.meta.data)[dalist & listquery]) else return(dalist & listquery)
}


.myrainbow <- function(n,nbring=c(), huerange = c(0,360) ){
	if (length(nbring) == 0) nbring =0
	if (length(n) != 1){
		nams <- n
		n <- length(n)
	}else{
		nams <- c()
	}
	if (nbring == 0) nbring <- 1 + floor(n / 12);
	fout = c();	
	phase = pi / (nbring);
	term = pi / (nbring);
	for(i in 0:(n-1)){
		fout <- c(fout, hcl((floor(i / nbring) * (huerange[2]- huerange[1]) / floor(n / nbring)) + huerange[1], 100 + 50 *sin(phase + term* (i%%nbring)) , 50 - 25 *cos(phase + term* (i%%nbring))  ))
	}
	if (is.null(nams)) return(fout)
	fout2 = c()
	for(i in 1:n) fout2[[nams[i]]] <- fout[i]
return(fout2)}

.changeStyle <- function(p, plot.attribs, classprefix=""){
	library(ggplot2)
	if ("flags" %in% names(plot.attribs)) flags.plot <- plot.attribs[["flags"]]
	else flags.plot <- c()
	if ("title" %in% names(plot.attribs)) p <- p + ggtitle(plot.attribs[["title"]])
	if ("xlabel" %in% names(plot.attribs)) p <- p + xlab(plot.attribs[["xlabel"]])
	else if (("no.xlabel" %in% flags.plot)||(paste(classprefix,"no.xlabel",sep=".") %in% flags.plot))  p <- p + xlab(NULL)
	if ("ylabel" %in% names(plot.attribs)) p <- p + ylab(plot.attribs[["ylabel"]])
	else if (("no.ylabel" %in% flags.plot)||(paste(classprefix,"no.ylabel",sep=".") %in% flags.plot))  p <- p + ylab(NULL)

	if ("no.legend" %in% flags.plot) p <- p + theme(legend.position="none")
	else{
		if ("legnames.size" %in% names(plot.attribs))  p <- p + theme(legend.text=element_text(size=plot.attribs[["legnames.size"]]))
	}

	# building xticks style
	themearg <- list()
	if (("no.xnamedticks" %in% flags.plot)||(paste(classprefix,"no.xnamedticks",sep=".") %in% flags.plot )) {
		p <- p + theme(axis.ticks.x=element_blank()) 
		p <- p + theme(axis.text.x=element_blank())
	}else{
		if ("no.xticks" %in% flags.plot) p <- p + theme(axis.ticks.x=element_blank())
		if ("no.xnames" %in% flags.plot) p <- p + theme(axis.text.x=element_blank())
		else if (!(("rot.xnames" %in% flags.plot)||(paste(classprefix,"rot.xnames",sep=".") %in% flags.plot))) themearg <- c(themearg,list(angle = 90, hjust = 1))
		if ("xnames.size" %in% names(plot.attribs)) themearg <- c(themearg,list(size= plot.attribs[["xnames.size"]]))
		if ("xnames.color" %in% names(plot.attribs)) themearg <- c(themearg,list(colour= plot.attribs[["xnames.color"]]))
		if ("xnames.bold" %in% flags.plot) themearg <- c(themearg,list(face="bold"))
	}
	if (length(themearg) > 0) p <- p + theme(axis.text.x= do.call(element_text, themearg))
	themearg <- list()
	if (("no.ynamedticks" %in% flags.plot)||(paste(classprefix,"no.ynamedticks",sep=".") %in% flags.plot )) {
		p <- p + theme(axis.ticks.y=element_blank()) 
		p <- p + theme(axis.text.y=element_blank())
	}else{
		if ("no.yticks" %in% flags.plot) p <- p + theme(axis.ticks.y=element_blank())
		if ("no.ynames" %in% flags.plot) p <- p + theme(axis.text.y=element_blank())
		else if ((("rot.ynames" %in% flags.plot)||(paste(classprefix,"rot.ynames",sep=".") %in% flags.plot))) themearg <- c(themearg,list(angle = 90, hjust = 1))
		if ("ynames.size" %in% names(plot.attribs)) themearg <- c(themearg,list(size= plot.attribs[["ynames.size"]]))
		if ("ynames.color" %in% names(plot.attribs)) themearg <- c(themearg,list(colour= plot.attribs[["ynames.color"]]))
		if ("ynames.bold" %in% flags.plot) themearg <- c(themearg,list(face="bold"))
	}
	if (length(themearg) > 0) p <- p + theme(axis.text.y= do.call(element_text, themearg))

	if ("scale.xrange" %in% names(plot.attribs)) p <- p + scale_x_continuous(limits = plot.attribs$scale.xrange)
	if ("scale.yrange" %in% names(plot.attribs)) p <- p + scale_y_continuous(limits = plot.attribs$scale.yrange)	

	if ("nb.col.legend" %in% names(plot.attribs)) p <- p + guides(color=guide_legend(ncol = plot.attribs$nb.col.legend),fill=guide_legend(ncol = plot.attribs$nb.col.legend))

	return(p)
}


