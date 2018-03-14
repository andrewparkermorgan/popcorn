# io.R

#' Encode sex as integer
#' 
#' @param x vector of sex labels, encoded in any sensible way
#' @return integer vector of sex codes: 0=uknown, 1=male, 2=female
#' @details Think of the result as the count of X chromosomes.
#' 
#' @export
.fix_sex <- function(x) {
	
	if (is.character(x) || is.factor(x))
		ifelse(grepl("^[fF]", x), 2L, ifelse(grepl("^[mM]", x), 1L, 0L))
	else
		return(as.integer(x))
	
}

#' Read a PLINK family (*.fam) file
#' 
#' @param ff path to a *.fam file
#' @return dataframe of sample info
#' @details A *.fam file should have columns family ID, individual ID, mom ID, dad ID, sex, phenotype.  PLINK uses
#'   special value -9 to denote missing phenotype, for some reason, while 0=missing for mom and dad ID.
#' 
#' @export
read_fam <- function(ff, ...) {
	
	fam <- tibble::as_tibble(read.table(ff, header = FALSE, stringsAsFactors = FALSE))
	cols <- c("fid","iid","mom","dad","sex","pheno")
	nc <- min(length(cols), ncol(fam))
	colnames(fam)[ 1:nc ] <- cols[1:nc]
	if ("mom" %in% colnames(fam)) {
		fam$mom <- as.character(fam$mom)
	}
	if ("dad" %in% colnames(fam)) {
		fam$dad <- as.character(fam$dad)
	}
	if ("sex" %in% colnames(fam)) {
		fam$sex <- .fix_sex(fam$sex)
	}
	return(fam)
	
}

#' Read a PLINK marker map (*.bim) file
#' 
#' @param ff path to a *.bim file
#' @return dataframe of marker info
#' @details A *.bim file should have columns chromosome, marker ID, genetic position (cM), physical position,
#'   allele 1 (A1=alternate allele), allele 2 (A2=reference allele).
#' 
#' @export
read_map <- function(ff, ...) {
	
	df <- tibble::as_tibble(read.table(ff, header = FALSE, stringsAsFactors = FALSE))
	cols <- c("chr","marker","cM","pos","A1","A2")
	nc <- min(length(cols), ncol(df))
	colnames(df)[ 1:nc ] <- cols[1:nc]
	df$chr <- mouser::factor_chrom(df$chr)
	return(df)
	
}

#' @export
write_fam <- function(ff, fname, ...) {
	
	fam <- tibble::as_tibble(ff)
	cols <- c("fid","iid","mom","dad","sex","pheno")
	nc <- min(length(cols), ncol(fam))
	colnames(fam)[ 1:nc ] <- cols[1:nc]
	if ("mom" %in% colnames(fam)) {
		fam$mom <- as.character(fam$mom)
	}
	if ("dad" %in% colnames(fam)) {
		fam$dad <- as.character(fam$dad)
	}
	if ("sex" %in% colnames(fam)) {
		fam$sex <- .fix_sex(fam$sex)
	}
	readr::write_tsv(ff, fname, col_names = FALSE)
	
}

#' @export
update_fam_file <- function(ff, fam, ...) {
	
	if (!inherits(ff, "data.frame"))
		oldfam <- read_fam(ff)
	else
		oldfam <- ff
	ii <- match(oldfam$iid, fam$iid)
	nona <- !is.na(ii)
	newfam <- oldfam
	newfam$fid[nona] <- as.character(fam$fid[ii[nona]])
	newfam$sex[nona] <- .fix_sex(fam$sex[ii[nona]])
	return(newfam)
	
}

#' @export
read_missingness_plink <- function(ff, ...) {
	
	df <- read.table(ff, header = TRUE, stringsAsFactors = FALSE)
	colnames(df) <- c("chr","marker","nmiss","nind","fmiss")
	df <- tibble::as_tibble(df)
	#print(head(df))
	df$ntyped <- with(df, nind-nmiss)
	df$chr <- factor_chrom(df$chr)
	return(df)
	
}

#' @export
read_freq_plink <- function(ff, ...) {
	
	df <- read.table(ff, header = TRUE, stringsAsFactors = FALSE)
	colnames(df) <- c("chr","marker","A1","A2","freq","nobs")
	df <- tibble::as_tibble(df)
	#print(head(df))
	df$chr <- factor_chrom(df$chr)
	return(df)
	
}

#' Read BEAGLE-format IBD segments
#' 
#' @export
read_beagle_ibd <- function(ff, map = NULL, expand = FALSE, as.ranges = FALSE, ...) {
	
	ibd <- readr::read_tsv(ff, col_names = FALSE)
	colnames(ibd) <- c("iid1","p1","iid2","p2","chr","start","end","lod")
	ibd$chr <- mouser::factor_chrom(ibd$chr)
	
	if (!is.null(map)) {
		start.cM <- dplyr::left_join(ibd, map[ ,c("chr","cM","pos") ], by = c("chr" = "chr", "start" = "pos"))
		end.cM <- dplyr::left_join(ibd, map[ ,c("chr","cM","pos") ], by = c("chr" = "chr", "end" = "pos"))
		ibd$start.cM <- start.cM$cM
		ibd$end.cM <- end.cM$cM
		ibd$width.cM <- with(ibd, end.cM - start.cM)
	}
	
	if (expand) {
		ibd2 <- ibd
		tmp <- ibd2$iid1
		ibd2$iid1 <- ibd2$iid2
		ibd2$iid2 <- tmp
		ibd <- dplyr::bind_rows(ibd, ibd2)
	}
	
	if (as.ranges) {
		ibd.gr <- GenomicRanges::makeGRangesFromDataFrame(ibd, keep.extra.columns = TRUE)
		return(ibd.gr)
	} else {
		return(ibd)
	}
	
}

make_ibd_matrix <- function(df, iids = NULL, ...) {
	
	if (is.null(iids))
		iids <- union(df$iid1, df$iid2)
	
	df$iid1 <- factor(df$iid1, iids)
	df$iid2 <- factor(df$iid2, iids)
	X <- with(df, tapply(width.cM, list(iid1, iid2), sum, na.rm = TRUE))
	return(X)
	
}

#' Read `bcftools roh` output (mode `-Or``)
#' 
#' @export
read_roh <- function(ff, force.chroms = TRUE, ...) {
	
	roh <- readr::read_tsv(ff, col_names = c("flag","iid","chr","start","end","width","nsites","qual"), comment = "#")
	roh$flag <- NULL
	if (force.chroms)
		roh$chr <- mouser::factor_chrom(roh$chr)
	return(roh)
	
}

# Read BLAT *.psl file
#'
#' @export
read_blat_psl <- function(ff, with.header = FALSE, force.chroms = TRUE, ...) {
	
	nskip <- if (with.header) 3 else 0
	df <- readr::read_tsv(ff, skip = nskip,
						  col_names = c("matches","mismatches","rep.matches","ns","qinserts","qbpins",
						  			  "tinserts","tbpins","strand","qname","qsize","qstart","qend",
						  			  "tname","tsize","tstart","tend","nblocks","blen","qstarts","tstarts"))
	if (force.chroms) {
		df$chr <- mouser::factor_chrom(df$tname)
		df$start <- df$tstart
		df$end <- df$tend
	}
	
	return(df)
		
}