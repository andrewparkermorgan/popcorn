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

#' Read ROH report from PLINK
#' @export
read_homozyg_plink <- function(ff, ...) {
	
	df <- read.table(ff, header = TRUE, stringsAsFactors = FALSE)
	colnames(df) <- c("fid","iid","phe","chr","left","right","start","end","width.kb","nsites","dens","phom","phet")
	df <- tibble::as_tibble(df)
	df <- df[ ,c("iid","chr","left","start","right","end","nsites","dens","phom","phet") ]
	df$width <- with(df, end-start)
	df$chr <- factor_chrom(df$chr)
	return(df)
	
}

#' Read estimated inbreeding report from PLINK
#' @export
read_het_plink <- function(ff, ...) {
	
	df <- read.table(ff, header = TRUE, stringsAsFactors = FALSE)
	colnames(df) <- c("fid","iid","obs","expect","nhet","fhat")
	df <- tibble::as_tibble(df)
	return(df)
	
}

#' Read BEAGLE-format IBD segments
#' 
#' @export
read_beagle_ibd <- function(ff, map = NULL, expand = FALSE, as.ranges = FALSE, ...) {
	
	ibd <- readr::read_tsv(ff, col_names = FALSE)[ ,1:8 ]
	colnames(ibd) <- c("iid1","p1","iid2","p2","chr","start","end","lod")
	ibd$chr <- mouser::factor_chrom(ibd$chr)
	map <- subset(map, !is.na(cM) & cM > 0)
	dups <- duplicated(map[ ,c("chr","pos") ])
	map <- map[ !dups, ]
	
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

#' @export
read_refinedIBD <- function(ff, expand = FALSE, row_ids = TRUE, ...) {
	
	ibd <- readr::read_tsv(ff, col_names = c("iid1","p1","iid2","p2","chr","start","end","lod","cM"))
	ibd$.row <- seq_len(nrow(ibd))
	if (expand) {
		ibd2 <- ibd
		tmp <- ibd2$iid1
		ibd2$iid1 <- ibd2$iid2
		ibd2$iid2 <- tmp
		ibd <- dplyr::bind_rows(ibd, ibd2)
	}
	
	return(ibd)
	
}

#' @export
read_ibdne_result <- function(ff, ...) {

	ne <- readr::read_tsv(ff, ...)
	colnames(ne) <- c("gen","Ne","lo","hi")[ seq_len(ncol(ne)) ]
	ne$gen <- as.integer(gsub("^G","", ne$gen))
	return(ne)
	
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
						  			  "tname","tsize","tstart","tend","nblocks","blen","qstarts","tstarts"),
						  col_types = "iiiiiiiicciiiciiiiccc")
	
	if (force.chroms) {
		df$chr <- mouser::factor_chrom(df$tname)
		df$start <- df$tstart
		df$end <- df$tend
	}
	
	return(df)
		
}

# Write genomic intervals in BED format
#' @export
write_bed <- function(x, ff = "/dev/stdout", extra = TRUE, is.zero.based = FALSE, ...) {
	
	if (inherits(x, "GRanges")) {
		x <- as.data.frame(x)
		colnames(x)[1] <- "chr"
		x$start <- pmax(x$start - 1L, 1L)
	}
	cn <- colnames(x)
	if ("seqnames" %in% cn) {
		cn[ cn == "seqnames" ] <- "chr"
		colnames(x) <- cn
	}
	
	if (!is.zero.based) {
		x$start <- pmax(x$start - 1L, 1L)
	}
	
	cols <- c("chr","start","end","name","score","strand")
	others <- cn[ !(cn %in% cols) ]
	keep.cols <- cols[ sort(which(cols %in% cn)) ]
	if (extra)
		write.cols <- c(keep.cols, others)
	else
		write.cols <- keep.cols
	write.table(x[ ,write.cols ], ff, 
				col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
	
}

#' Read genomic intervals in BED format
#' @export
read_bed <- function(ff, ..., force.chroms = TRUE) {
	
	df <- readr::read_tsv(ff, col_names = FALSE, comment = "#", ...)
	idx <- 1:min(ncol(df), 6)
	colnames(df)[ idx ] <- c("chr","start","end","name","score","strand")[ idx ]
	if (force.chroms)
		df$chr <- mouser::factor_chrom(df$chr)
	return(df)
	
}

#' Read output from my \code{coverage.py} script
#' 
#' @param path name of file to read
#' @param scrub if \code{TRUE}, remove trailing suffixes (\code{*.sorted}, \code{*.merged} etc) from individual IDs
#' @param as.ranges if \code{TRUE}, convert to a \code{GRanges} object
#' 
#' @details Expected input is a BED-like file with columns *chr*, *start*, *end*, *iid*, *raw read count*,
#' 	*normalized read count* and optionally *count of MQ0 reads*, *proportion of total reads which are MQ0*.
#' 
#' @return a dataframe
read_coverage <- function(path, scrub = TRUE, as.ranges = FALSE, has.MQ0 = FALSE, ...) {
	
	counts <- tryCatch({
		readr::read_tsv(path, col_names = FALSE, comment = "#")
	},
	error = function(e) {
		read.table(path, header = FALSE, comment.char = "#")
	})
	colnames(counts)[1:6] <- c("chr","start","end","iid","nreads","depth")
	if (has.MQ0)
		colnames(counts)[7:8] <- c("MQ0","MQ0.prop")
	
	if (scrub)
		counts$iid <- gsub("\\.\\w+$", "", counts$iid)
	
	if (as.ranges)
		counts <- GenomicRanges::makeGRangesFromDataFrame(counts, starts.in.df.are.0based = TRUE,
														  keep.extra.columns = TRUE)
	
	attr(counts, "filename") <- normalizePath(path)
	return(counts)
	
}

#' @export
read_psmc_result <- function(ff, niter = 25, mu = 1e-8, s = 100, ...) {
	
	cmd <- paste0("awk '/RD.", niter, "/{f=1}/\\/\\//{f=0}f&&/RS/{print $0;}' ", ff)
	rez <- readr::read_tsv(pipe(cmd), col_names = c("what","k","t_k","lambda_k","pi_k","sum_A_kl","A_kk"))
	maxk <- max(rez$k)
	nboot <- sum(rez$k == maxk)
	rez$repno <- floor(seq_len(nrow(rez))/maxk)
	
	cmd <- paste0("awk 'BEGIN{OFS=\"\\t\";}/RD.", niter, "/{f=1}/\\/\\//{f=0}f&&/TR/{print $2,$3}' ", ff)
	pis <- readr::read_tsv(pipe(cmd), col_names = c("theta_0","rho_0"))
	pis$repno <- seq_len(nrow(pis))-1

	rez <- dplyr::left_join(rez, pis) %>%
		dplyr::mutate(Nk = (lambda_k*theta_0/s)/(4*mu))
	return(rez)
	
}