## admixture.R
# functions for handling output from ADMIXTURE (Alexander DH 2009 Genome Res)

#' Read ancestry proportions (`Q-matrix`) from ADMIXTURE run
#' 
#' @param infile path to matrix of ancestry proportions (*.Q file) from ADMIXTURE run
#' @param K number of clusters; if \code{NULL} (default), attempt to infer from filename
#' @param iids sample IDs (in same order as rows in matrix); if \code{NULL} (default),
#'  look for corresponding *.fam file
#' @param ... extra arguments, ignored
#' @details If \code{K} and \code{iids} are both \code{NULL} (default), assumes that
#'  estimation was performed from a PLINK *.bed file, and that the corresponding *.fam
#'  file resides in same directory as output. If standard errors for ancestrty proportions
#'  are availbale (*.Q_se file), they are read as well and stored in \code{attr(result, "se")}.
#' 
#' @export
read_Q_matrix <- function(infile, K = NULL, iids = NULL, ...) {
	
	infile <- infile[1]
	fam.path <- NULL
	if (is.null(iids)) {
		if (grepl("\\.\\d+\\.Q$", infile)) {
			fam.path <- gsub("\\.\\d+\\.Q$",".fam", infile)
		}
		else {
			fam.path <- paste0(infile, ".fam")
			infile <- paste0(infile, ".Q")
		}
	}
	
	sefile <- paste0(infile, "_se")
	do.se <- file.exists(sefile)
	
	# figure out how many clusters
	if (is.null(K))
		K <- as.integer(stringr::str_match(infile, "([0-9]+)\\.Q$")[,2])
	
	message("Reading ancestry proportions from <", infile, ">...")
	Q <- readr::read_delim(infile, delim =" ", col_names = FALSE)
	
	if (do.se) {
		message("Reading standard errors from <", sefile, ">...")
		se <- readr::read_delim(sefile, delim = " ", col_names = FALSE)
	}
	
	if (!is.null(fam.path)){
		message("Reading sample metadata from <", fam.path, ">...")
		fam <- read_fam(fam.path)
	}
	else {
		message("Using input sample IDs ...")
		fam <- tibble::tibble(fid = iids, iid = iids, mom = 0, dad = 0, sex = 0, pheno = -9)
	}
		
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

#' @export
read_Fst_matrix <- function(infile, melted = TRUE, lower.tri = TRUE, ...) {
	
	infile <- infile[1]
	
	message("Reading Fst matrix from <", infile, "> ...")
	lines <- readLines(infile)
	popfile <- gsub("\\.\\d+\\.out$", ".pop", infile)
	
	keep <- grepl("^Pop\\d+", lines, perl = TRUE)
	lines <- lines[keep]
	ll <- stringr::str_split(lines, "\\t")
	
	npop <- length(ll[[1]])+1
	if (file.exists(popfile)) {
		message("Reading population labels from <", popfile, "> ...")
		pops <- readr::read_tsv(popfile, col_names = "pop", na = "-", progress = FALSE)$pop
		pops <- pops[ !is.na(pops) ]
		pops <- pops[ !duplicated(pops) ]
		po <- sort(pops)
	}
	else {
		message("Using default population labels from ADMIXTURE ...")
		pops <- ll[[1]]
		po <- pops
	}
	
	Fst <- matrix(0.0, nrow = npop, ncol = npop)
	for (ii in seq_along(ll[-(1:2)])) {
		values <- as.numeric( ll[[ ii+2 ]][ -1 ] )
		Fst[ ii+1,1:length(values) ] <- values
	}
	
	Fst <- Fst + t(Fst)
	rownames(Fst) <- colnames(Fst) <- pops
	cl <- stats::hclust( stats::as.dist(Fst+0.001), method = "centroid" )
	o <- cl$order
	po <- po[o]
	Fst <- Fst[ po,po ]
	
	if (melted) {
		
		if (lower.tri) {
			Fst[ upper.tri(Fst) ] <- NA
		}
		Fst <- dplyr::bind_cols(pop = rownames(Fst), tibble::as_tibble(Fst))
		Fst <- tidyr::gather(Fst, key = other, value = fst, -pop)
		Fst <- subset(Fst, !is.na(fst))
		Fst$pop <- factor(Fst$pop, levels = po)
		Fst$other <- factor(Fst$other, levels = po)
		
	}
	
	attr(Fst, "source") <- infile
	attr(Fst, "clustering") <- cl
	return(Fst)
	
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
#' @param ... extra arguments; if first one is a data frame, treat it as additional sample metadata,
#'   with column \code{"iid"} as primary key
#' @return a dataframe of admixture proportions (and standard errors if available),
#'   one row per component-individual pair
#' 
#' @export
tidy <- function(Q, ...){
	UseMethod("tidy")
}

#' @export
tidy.admixture_Q <- function(Q, ...) {

	if (!inherits(Q, "admixture_Q"))
		warning("Q may not be a properly-formatted ADMIXTURE result.")
	
	K <- attr(Q, "K")
	qq <- dplyr:: mutate(tidyr::gather(Q, key = variable, value = value, -iid, -fid), K = K)
	
	if (!is.null(attr(Q, "se"))) {
		se <- attr(Q, "se")
		se <- dplyr::mutate(tidyr::gather(se, key = variable, value = se, -iid, -fid), K = K)
		qq <- dplyr::inner_join(qq, se)
	}
	
	extras <- list(...)
	if (length(extras)) {
		if (is.data.frame(extras[[1]])) {
			meta <- tibble::as_tibble(extras[[1]])
			qq <- dplyr::left_join(qq, meta, by = "iid")
		}
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
#' @param sort_order apply pre-determined sort order to individuals, eg. by cluster membership
#' @param pop_order named vector to map arbitrary population labels ('pop1') onto meaningful ones
#' @param ... extra arguments, ignored
#' 
#' @export
plot_admixture <- function(Q, label = FALSE, sort_order = NULL, pop_order = NULL, ...) {
	
	
	if (inherits(Q, "admixture_Q"))
		qq <- tidy.admixture_Q(Q, ...)
	else
		qq <- tibble::as_tibble(Q)
	
	if (!is.null(pop_order)) {
		qq$variable <- factor(qq$variable, labels = pop_order)
		levels(qq$variable) <- names(pop_order)
	}
	if (!is.null(sort_order)) {
		qq$iid <- factor(qq$iid, sort_order)
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

#' Make stacked-bar plot of ancestry proportions, with population labels
#' 
#' @param Q a matrix of admixture proportions, melted or not; must have a column `pop` with population labels
#' @param label show individual IDs along axis
#' @param sort_order pre-determined sort order for individuals, eg. by cluster membership
#' @param pop_order pre-determined sort order for populations
#' @param cluster_order named vector to map arbitrary ancestry-component labels ('pop1') onto meaningful ones
#' @param label_pops should population labels be drawn?
#' @param label_ind should individuals be labelled along x-axis (probaly not, if more than a few samples)
#' @param nudge extra buffer around population labels
#' @param ymax how much extra space to allow for population labels; about 1.5 is probably enough
#' @param ... extra arguments, ignored
#' 
#' @export
plot_admixture_labelled <- function(Q, sort_order = NULL, pop_order = NULL, cluster_order = NULL, label_pops = TRUE, label_ind = FALSE, nudge = 0.05, ymax = 1.5, ...) {
	
	if (inherits(Q, "admixture_Q"))
		qq <- tidy.admixture_Q(Q, ...)
	else
		qq <- tibble::as_tibble(Q)
	
	if (!is.null(cluster_order)) {
		qq$variable <- factor(qq$variable, labels = cluster_order)
		if (!is.null(names(cluster_order))) {
			levels(qq$variable) <- names(cluster_order)
		} else {
			levels(qq$variable) <- cluster_order
		}
	} else {
		qq$variable <- factor(qq$variable)
	}
	if (!is.null(sort_order)) {
		qq$iid <- factor(qq$iid, sort_order)
	} else {
		qq$iid <- factor(qq$iid)
		sort_order <- levels(qq$iid)
	}
	if (!is.null(pop_order)) {
		qq$pop <- factor(qq$pop, levels = pop_order)
	} else {
		qq$pop <- factor(qq$pop)
		pop_order <- levels(qq$pop)
	}
	
	qq <- dplyr::arrange(qq, variable)
	qq1 <- subset(qq, as.numeric(variable) == 1)
	so <- with(qq1, order(pop, iid))
	new_so <- as.character(qq1$iid[so])
	qq$iid <- factor(qq$iid, new_so)
	qq1$iid <- factor(qq1$iid, new_so)
	qq1 <- dplyr::arrange(qq1, iid)
	
	pop_labs <- pop_order
	nudge_by <- length(new_so)*nudge
	pop_spacing <- (length(new_so)-nudge_by)/length(pop_labs)
	pop_pos <- (seq_along(pop_labs)-1)*pop_spacing+nudge_by
	#print(setNames(pop_pos, pop_labs))
	first_ind <- sapply(pop_labs, function(f) min(which(qq1$pop == f)))
	
	#print(first_ind)
	pop_df <- tibble::tibble(pop = pop_labs,
							 pos = pop_pos,
							 first = first_ind)
	
	p0 <- ggplot2::ggplot(qq) +
		ggplot2::geom_bar(ggplot2::aes(x = iid, y = value, fill = variable), stat = "identity") +
		ggplot2::geom_segment(data = pop_df, 
							  ggplot2::aes(x = first-1, xend = first-1, y = 0, yend = 1)) +
		ggplot2::scale_y_continuous(limits = c(0, ymax),
									breaks = c(0, 0.5, 1.0)) +
		theme_admixture()
	
	if (label_pops) {
		p0 <- p0 +
			ggplot2::geom_text(data = pop_df, 
							   ggplot2::aes(x = pos, y = 1.2, label = pop),
							   angle = 90, hjust = 0, nudge_x = 0.05, nudge_y = 0.02) +
			ggplot2::geom_segment(data = pop_df, 
								  ggplot2::aes(x = pos, xend = first-1, y = 1.2, yend = 1.0))
	}
	
	if (!label_ind)
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