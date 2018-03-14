## admixture.R
# functions for handling output from ADMIXTURE (Alexander DH 2009 Genome Res)

#' Read ancestry proportions (`Q-matrix`) from ADMIXTURE run
#' 
#' @param infile path to matrix of ancestry proportions (*.Q file) from ADMIXTURE run
#' @param ... extra arguments, ignored
#' @details Assumes that estimation was performed from a PLINK *.bed file,
#' 	and that the corresponding *.fam file resides in same directory as output.
#' 	If standard errors for ancestrty proportions are availbale (*.Q_se file),
#' 	they are read as well and stored in \code{attr(result, "se")}.
#' 
#' @export
read_Q_matrix <- function(infile, ...) {
	
	infile <- infile[1]
	if (grepl("\\.\\d+\\.Q$", infile)) {
		fam.path <- gsub("\\.\\d+\\.Q$",".fam", infile)
	}
	else {
		fam.path <- paste0(infile, ".fam")
		infile <- paste0(infile, ".Q")
	}
	sefile <- paste0(infile, "_se")
	do.se <- file.exists(sefile)
	
	# figure out how many clusters
	K <- as.integer(stringr::str_match(infile, "([0-9]+)\\.Q$")[,2])
	
	message("Reading ancestry proportions from <", infile, ">...")
	Q <- readr::read_delim(infile, delim =" ", col_names = FALSE)
	
	if (do.se) {
		message("Reading standard errors from <", sefile, ">...")
		se <- readr::read_delim(sefile, delim = " ", col_names = FALSE)
	}
	
	message("Reading sample metadata from <", fam.path, ">...")
	fam <- read_fam(fam.path)
	
	#rownames(Q) <- fam$iid
	colnames(Q) <- paste0("pop", seq_len(ncol(Q)))
	Q <- dplyr::bind_cols(iid = fam$iid, fid = as.character(fam$fid), Q)
	
	if (do.se) {
		colnames(se) <- paste0("pop", seq_len(ncol(se)))
		se <- dplyr::bind_cols(iid = fam$iid, fid = as.character(fam$fid), se)
		attr(Q, "se") <- se
	}
	
	class(Q) <- c("admixture_Q", class(Q))
	attr(Q, "K") <- K
	return(Q)
	
}

# deprecated and ill-conceived
read_Q_bootstraps <- function(prefix, k, ...) {
	
	where <- dirname(prefix)
	bname <- basename(prefix)
	pattern <- paste0(bname, "\\..+\\.", k, "\\.Q")
	print(pattern)
	ff <- list.files(path = where, pattern = pattern, full.names = TRUE)
	print(ff)
	#id <- seq_along(ff)
	#df <- tibble::tibble(run = id)
	lapply(ff, read_Q_matrix)
	
}

# deprecated and ill-conceived
bootstrap_Q <- function(ql, ...) {
	
	K <- attr(ql[[1]], "K")
	iids <- ql[[1]][ ,c("iid","fid") ]
	
	Q <- simplify2array(lapply(ql, function(f) as.matrix(f[,-(1:2)])))
	mu <- apply(Q, c(1,2), mean, na.rm = TRUE)
	sigma <- apply(Q, c(1,2), sd, na.rm = TRUE)/sqrt(dim(Q)[3])
	
	mu <- dplyr::bind_cols(iids, tibble::as_tibble(mu))
	sigma <- dplyr::bind_cols(iids, tibble::as_tibble(sigma))
	
	mu <- mutate(tidyr::gather(mu, key = pop, value = mu, -iid, -fid), K = K)
	sigma <- mutate(tidyr::gather(sigma, key = pop, value = se, -iid, -fid), K = K)
	
	dplyr::inner_join(mu, sigma)
	
}

#' Convert admixture proportions to tidy format.
#' 
#' @param Q a matrix of admixture proprtions
#' @param meta additional sample metadata, with column \code{"iid"} as primary key
#' @return a dataframe of admixture proportions (and standard errors if available),
#'   one row per component-individual pair
#' @param ... extra arguments, ignored
#' 
#' @export
tidy.admixture_Q <- function(Q, meta = NULL, ...) {

	if (!inherits(Q, "admixture_Q"))
		warning("Q may not be a properly-formatted ADMIXTURE result.")
	
	K <- attr(Q, "K")
	qq <- mutate(tidyr::gather(Q, key = variable, value = value, -iid, -fid), K = K)
	
	if (!is.null(attr(Q, "se"))) {
		se <- attr(Q, "se")
		se <- mutate(tidyr::gather(se, key = variable, value = se, -iid, -fid), K = K)
		qq <- dplyr::inner_join(qq, se)
	}
	
	if (!is.null(meta)) {
		qq <- dplyr::left_join(qq, meta, by = "iid")
	}
	
	return(qq)
		
}

#' Read log-likelihoods from ADMIXTURE runs
#' @export
read_loglik_files <- function(prefix, ...) {
	
	d <- dirname(prefix)
	f <- basename(prefix)
	ff <- list.files(d, paste0(f, "\\.[0-9]+\\.llik"), full.names = TRUE)
	k <- stringr::str_match(ff, "\\.([0-9]+)\\.llik")[,2]
	llik <- sapply(ff, readLines)
	tibble::tibble(K = as.integer(k), llik = as.numeric(llik))
	
}

#' Make standard stacked-bars plot of ancestry proportions.
#' 
#' @param Q a matrix of admixture proportions, melted or not
#' @param label show individual IDs along axis
#' @param sort.order apply pre-determined sort order to individuals, eg. by cluster membership
#' @param pop.order named vector to map arbitrary population labels ('pop1') onto meaningful ones
#' @param ... extra arguments, ignored
#' 
#' @export
plot_admixture <- function(Q, label = FALSE, sort.order = NULL, pop.order = NULL, ...) {
	
	
	if (inherits(Q, "admixture_Q"))
		qq <- tidy.admixture_Q(Q, meta = meta)
	else
		qq <- tibble::as_tibble(Q)
	
	if (!is.null(pop.order)) {
		qq$variable <- factor(qq$variable, labels = pop.order)
		levels(qq$variable) <- names(pop.order)
	}
	if (!is.null(sort.order)) {
		qq$iid <- factor(qq$iid, sort.order)
	}
	
	qq <- dplyr::arrange(qq, variable)
	p0 <- ggplot2::ggplot(qq) +
		ggplot2::geom_bar(ggplot2::aes(x = iid, y = value, fill = variable), stat = "identity") +
		theme_admixture()
	
	if (!label)
		p0 <- p0 +
		ggplot2::theme(axis.ticks.x = ggplot2::element_blank(),
					   axis.text.x = ggplot2::element_blank())
	
	return(p0)
	
}

#' Sort samples by admixture components using crude heuristics
#' 
#' @param Q an ADMIXTURE Q-matrix from \code{read_Q_matrix()}
#' @param method clustering heuristic; currently only works with "all"
#' @param ... extra arguments, ignored
#' @return character vector of individual IDs in sorted order
#' 
#' @export
sort_by_cluster <- function(Q, method = "all", ...) {
	
	if (!inherits(Q, "admixture_Q"))
		stop("Only works for bona fide un-melted ADMIXTURE result.")
	
	pop.cols <- paste0("pop", seq_len(attr(Q, "K")))
	QQ <- matrix(as.matrix(Q[ ,pop.cols ]), nrow = nrow(Q), dimnames = list(Q$iid))
	rez <- .sort_indiv(QQ, sortind = method)
	so <- rownames(rez$dframe)
	return(so)
	
}

#' A ggplot theme for ADMIXTURE plots
#' @export
theme_admixture <- function(...) {
	
	ggplot2::theme_bw(...) +
		ggplot2::theme(...,
					   axis.title.x = ggplot2::element_blank(),
					   axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 6),
					   panel.background = ggplot2::element_blank(),
					   panel.grid = ggplot2::element_blank(),
					   panel.border = ggplot2::element_blank(),
					   panel.spacing.x = grid::unit(0.15, "lines"),
					   strip.background = ggplot2::element_rect(colour = NA)
		)

}

#' Internal helper function for sorting rows in ADMIXTURE Q-matrix
#' @details Code stolen whole-cloth from <https://github.com/royfrancis/pophelper/blob/master/R/pophelper.R> 
#' @export
.sort_indiv <- function(dframe = NULL, grplab = NA, selgrp = NA, sortind = NA, grplabpos = NA, linepos = NA)
{
	if(is.null(dframe)) stop("sortInd: Argument 'dframe' is empty.")
	if(!all(is.na(grplab))) 
	{
		if(is.na(selgrp)) selgrp <- names(grplab)[1]
		if(!is.character(selgrp)) stop("sortInd: Argument 'selgrp' must be a character datatype.")
		if(length(selgrp)>1) stop("sortInd: Argument 'selgrp' must be a character datatype of length 1.")
		if(!any(selgrp %in% names(grplab))) stop(paste0("sortInd: Argument 'selgrp' contains (",selgrp,") which is not in the 'grplab' titles (",paste0(names(grpla),collapse=", "),")."))
	}
	if(all(!is.na(sortind)))
	{
		if(length(sortind) > 1) stop("sortInd: Argument 'sortind' must be of length 1. Use 'all','label' or a cluster name like 'Cluster1'.")
		#if(sortind != "all" && sortind != "label" && !grepl("Cluster[0-9+]",sortind)) stop("sortInd: Argument 'sortind' must be set to 'all', 'label' or a cluster like 'Cluster1'.")
	}
	
	#sorting without grplab
	if(any(is.na(grplab)))
	{
		if(!is.na(sortind))
		{
			if(sortind=="all")
			{
				maxval <- apply(dframe,1,max)
				matchval <- vector(length=nrow(dframe))
				for(j in 1:nrow(dframe)) matchval[j] <- match(maxval[j],dframe[j,])
				dftemp <- dframe
				dftemp$maxval <- maxval
				dftemp$matchval <- matchval
				
				dframe <- dframe[with(dftemp,order(matchval,-maxval)),,drop=FALSE]
			}
			
			if(sortind=="label") dframe[order(rownames(dframe)),,drop=FALSE]
			
			if(sortind!="all" && sortind!="label")
			{
				if(!(sortind %in% colnames(dframe))) stop(paste0("sortInd: 'sortind' value (",sortind,") not found in file header (",paste0(colnames(dframe),collapse=", "),")."))
				dframe <- dframe[order(dframe[[sortind]]),,drop=FALSE]
			}
		}
		labelposbind <- NA
		markerposbind <- NA
	}
	
	#sorting with grplab
	if(!all(is.na(grplab)))
	{
		gnames <- names(grplab)
		
		if(!is.na(sortind))
		{
			if(sortind=="all")
			{
				maxval <- apply(dframe,1,max)
				matchval <- vector(length=nrow(dframe))
				for(j in 1:nrow(dframe)) matchval[j] <- match(maxval[j],dframe[j,])
				dftemp <- dframe
				dftemp$maxval <- maxval
				dftemp$matchval <- matchval
				
				if(length(intersect(colnames(dftemp),colnames(grplab)))!=0) stop(paste0("sortInd: One or more header labels in the run file are duplicated in grplab header. Change labels to be unique. Following are the duplicate label(s): (",paste0(intersect(colnames(dftemp),colnames(grplab)),collapse=", "),")."))
				dftemp <- cbind(dftemp,grplab)
				
				rle1 <- rle(as.character(unlist(grplab[selgrp])))
				grplabnames <- rle1$values
				tovec <- cumsum(rle1$lengths)
				fromvec <- (tovec - rle1$lengths)+1
				dftemplist <- vector("list",length=length(grplabnames))
				for(k in 1:length(tovec))
				{
					dftemp1 <- dftemp[fromvec[k]:tovec[k],,drop=FALSE]
					dftemp1$grp <- NULL
					dftemplist[[k]] <- dftemp1[with(dftemp1,order(matchval,-maxval)),,drop=FALSE]
				}
				
				dframe <- do.call("rbind",dftemplist)
				grplab <- dframe[,gnames,drop=FALSE]
				dframe[,gnames] <- NULL
				dframe$maxval <- NULL
				dframe$matchval <- NULL
			}
			
			if(sortind=="label")
			{
				dftemp <- dframe
				if(length(intersect(colnames(dftemp),colnames(grplab)))!=0) stop(paste0("sortInd: One or more header labels in the run file are duplicated in grplab header. Change labels to be unique. Following are the duplicate label(s): (",paste0(intersect(colnames(dftemp),colnames(grplab)),collapse=", "),")."))
				dftemp <- cbind(dftemp,grplab)
				
				rle1 <- rle(as.character(unlist(grplab[selgrp])))
				grplabnames <- rle1$values
				tovec <- cumsum(rle1$lengths)
				fromvec <- (tovec - rle1$lengths)+1
				dftemplist <- vector("list",length=length(grplabnames))
				for(k in 1:length(tovec))
				{
					dftemp1 <- dftemp[fromvec[k]:tovec[k],,drop=FALSE]
					dftemp1$grp <- NULL
					dftemplist[[k]] <- dftemp1[order(rownames(dftemp1)),,drop=FALSE]
				}
				
				dframe <- do.call("rbind",dftemplist)
				grplab <- dframe[,gnames,drop=FALSE]
				dframe[,gnames] <- NULL
			}
			
			
			if(sortind!="all" && sortind!="label")
			{
				if(!(sortind %in% colnames(dframe))) stop(paste0("sortInd: 'sortind' value (",sortind,") not found in file header (",paste0(colnames(dframe),collapse=", "),")."))
				clnum <- which(sortind==colnames(dframe))
				dftemp <- dframe
				if(length(intersect(colnames(dftemp),colnames(grplab)))!=0) stop(paste0("sortInd: One or more header labels in the run file are duplicated in grplab header. Change labels to be unique. Following are the duplicate label(s): (",paste0(intersect(colnames(dftemp),colnames(grplab)),collapse=", "),")."))
				dftemp <- cbind(dftemp,grplab)
				
				rle1 <- rle(as.character(unlist(grplab[selgrp])))
				grplabnames <- rle1$values
				tovec <- cumsum(rle1$lengths)
				fromvec <- (tovec - rle1$lengths)+1
				dftemplist <- vector("list",length=length(grplabnames))
				for(k in 1:length(tovec))
				{
					dftemp1 <- dftemp[fromvec[k]:tovec[k],,drop=FALSE]
					dftemp1$grp <- NULL
					dftemplist[[k]] <- dftemp1[order(dftemp1[[sortind]]),,drop=FALSE]
				}
				
				dframe <- do.call("rbind",dftemplist)
				grplab <- dframe[,gnames,drop=FALSE]
				dframe[,gnames] <- NULL
			}
		}
		
		#create labelpos and markerpos
		gnames <- names(grplab)
		markerlist <- vector("list",length=length(gnames))
		labellist <- vector("list",length=length(gnames))
		intlablist <- vector("list",length=length(gnames))
		for(k in 1:length(gnames))
		{
			rlegrp <- rle(unlist(grplab[gnames[k]]))
			labelpos <- data.frame(label=rlegrp$values,freq=rlegrp$lengths,stringsAsFactors=F)
			markerpos <- data.frame(markerxpos=c(0,cumsum(labelpos$freq)),stringsAsFactors=F)
			labelpos$labxpos <- round((diff(markerpos$markerxpos)/2)+markerpos$markerxpos[1:length(markerpos$markerxpos)-1],1)
			labelpos$labypos <- rep(grplabpos,nrow(labelpos))
			markerpos$temp <- factor(rep(1,nrow(markerpos)))
			markerpos$markerypos <- rep(linepos,nrow(markerpos))
			
			markerpos$count <- gnames[k]
			markerlist[[k]] <- markerpos
			labelpos$count <- gnames[k]
			labellist[[k]] <- labelpos
		}
		
		markerposbind <- do.call("rbind",markerlist)
		markerposbind$count <- factor(markerposbind$count,levels=gnames)
		labelposbind <- do.call("rbind",labellist)
		labelposbind$count <- factor(labelposbind$count,levels=gnames)
		
		rownames(markerposbind) <- 1:nrow(markerposbind)
		rownames(labelposbind) <- 1:nrow(labelposbind)
		
		#adjust divider position
		markerposbind$divxpos <- markerposbind$markerxpos+0.5
	}
	
	
	
	return(list(dframe=dframe,grplab=grplab,labelpos=labelposbind,markerpos=markerposbind))
}