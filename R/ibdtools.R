# ibdtools.R
# functions for dealing with IBD segments

#' @export
make_ibd_matrix <- function(df, iids = NULL, as_matrix = FALSE, keep_diagonal = FALSE, dedup = FALSE, what = c("cM","kb"), ...) {
	
	if (is.null(iids))
		iids <- union(df$iid1, df$iid2)

	## should we tally up genetic (cM) or physical (kb) segment lengths?
	what <- match.arg(what)
	if (what == "cM") {
		df$score <- df$cM
	}
	else if (what == "kb") {
		df$score <- with(df, end/1e3-start/1e3)
	}
	else {
		stop("Specify IBD to be measured in either 'kb' or 'cM'.")
	}
	
	## attempt to merge overlapping redundant segments, only useful in cases of strange ploidy (eg. Plasmodium)
	if (dedup) {
		
		if (what != "kb")
			stop("If deduplicating, can only score on 'kb' to avoid ambiguity on genetic <-> physical map")
		
		.segment_to_ranges <- function(x) {
			oiid1 <- x$iid1[1]
			oiid2 <- x$iid2[1]
			gr <- GenomicRanges::makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
			gr <- GenomicRanges::reduce(gr)
			#print(c(length(gr), nrow(x)))
			rez <- GenomicRanges::as.data.frame(gr) %>%
				tibble::as_tibble() %>%
				dplyr::rename(chr = seqnames) %>%
				dplyr::select(chr, start, end) %>%
				dplyr::mutate(score = (end/1e3-start/1e3)) %>%
				dplyr::mutate(iid = oiid1, iid2 = oiid2)
			#print(rez)
			return(rez)
		}
		
		df <- dplyr::group_by(df, iid1, iid2) %>%
			dplyr::do(.segment_to_ranges(.))

	}
		
	df$iid1 <- factor(df$iid1, iids)
	df$iid2 <- factor(df$iid2, iids)
	shared <- with(df, tapply(score, list(iid1,iid2), sum))
	shared[ is.na(shared) ] <- 0
	
	if (as_matrix) {
		if (!keep_diagonal)
			diag(shared) <- NA
		return(shared)
	}
	else {
		
		longest <- with(df, tapply(score, list(iid1,iid2), max))
		longest[ is.na(longest) ] <- 0
		howmany <- with(df, tapply(score, list(iid1,iid2), length))
		howmany[ is.na(howmany) ] <- 0
		
		rez <- tibble::as_tibble(shared, rownames = "iid1") %>%
			tidyr::gather(iid2, score, -iid1) %>%
			dplyr::mutate(iid1 = factor(iid1, iids),
						  iid2 = factor(iid2, iids))
		
		rez2 <- tibble::as_tibble(longest, rownames = "iid1") %>%
			tidyr::gather(iid2, longest, -iid1) %>%
			dplyr::mutate(iid1 = factor(iid1, iids),
						  iid2 = factor(iid2, iids))
		
		rez3 <- tibble::as_tibble(howmany, rownames = "iid1") %>%
			tidyr::gather(iid2, nseg, -iid1) %>%
			dplyr::mutate(iid1 = factor(iid1, iids),
						  iid2 = factor(iid2, iids))
		
		rez <- dplyr::left_join(rez, rez2) %>%
			dplyr::left_join(rez3)
		
		if (!keep_diagonal)
			rez <- subset(rez, iid1 != iid2)
		
		if (what == "cM")
			rez <- dplyr::rename(rez, cM = score)
		
		return(rez)
		
	}
	
}