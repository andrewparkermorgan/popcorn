# vcfutils.R
# Some utilities for interacting with genetic variation data stored in VCF format, mostly via command-line tools.

#' Perform PCA on a VCF file using `akt pca` and consume result
#' 
#' @param vcf path to VCF file; must be indexed
#' @param samples optional character vector of sample IDs to which the analysis should be restricted
#' @param K number of PCs to compute
#' @param region genomic region on which to perform PCA, as \code{samtools}-style region expression (eg. "chr1:1-10")
#' @param ... ignored
#' 
#' @return An object of classes \code{data.frame} and \code{pca_result} whose columns are the projections of
#' 	samples onto PCs.  Eigenvalues, expressed as proportion of variance explained, are returned as an attribute.
#' 	
#' @details This function generates a command-line call to \code{akt pca}, whose output is directed to a set of
#' 	temporary files and then read into the \code{R} session.
#' 	 	
#' @export
pca_vcf <- function(vcf, samples = NULL, skeleton = NULL, K = 20, meta = NULL, region = NULL, useX = FALSE, ...) {
	
	ff <- tempfile()
	ffev <- paste0(ff, ".ev")
	cmd <- paste0("akt pca -N ", K[1], " -F ", ffev)
	
	if (!is.null(region)) {
		cmd <- paste0(cmd, " -r ", as.character(region)[1])
	}
	else if (!useX) {
		region <- paste0(mouser::chromnames()[1:19], collapse = ",")
		cmd <- paste0(cmd, " -r ", as.character(region)[1])
	}
	
	if (!is.null(skeleton)) {
		rfile <- paste0(ff, ".skeleton")
		wfile <- paste0(ff, ".vcf")
		write.table(as.character(skeleton), rfile, row.names = FALSE, col.names = FALSE, quote = FALSE)
		cmd1 <- paste0(cmd, " -S ", rfile, " -o ", wfile, " --force ", vcf[1], " >/dev/null")
		rez <- system(cmd1)
		rez <- system(paste0("bgzip ", wfile, " && tabix ", wfile, ".gz"))
		cmd <- paste0(cmd, " -W ", wfile, ".gz")
	}
	
	if (!is.null(samples)) {
		sfile <- paste0(ff, ".samples")
		write.table(as.character(samples), sfile, row.names = FALSE, col.names = FALSE, quote = FALSE)
		cmd <- paste0(cmd, " -S ", sfile)
	}
	
	cmd <- paste0(cmd, " --force ", vcf[1], " >", ff)
	message("The command will be:\n", cmd)
	
	rez <- system(cmd)
	
	pc <- tibble::as_tibble(read.table(ff, col.names = c("iid", paste0("PC", seq_len(K)))))
	eigvals <- as.numeric(readr::read_lines(ffev))^2
	if (!is.null(meta)) {
		pc <- dplyr::left_join(pc, meta)
	}
	attr(pc, "eigvals") <- eigvals
	attr(pc, "explained") <- eigvals/sum(eigvals)
	attr(pc, "result") <- ff
	class(pc) <- c("pca_result", class(pc))
	return(pc)
	
}

#' @export
write_pca_result <- function(rez, ff, ...) {
	readr::write_tsv(rez, paste0(ff, ".pca"))
	readr::write_lines(attr(rez, "explained"), paste0(ff, ".pca.ev"))
}

#' @export
read_pca_result <- function(ff, ...) {
	ff <- gsub("\\.pca$", "", ff)
	pcs <- readr::read_tsv(paste0(ff, ".pca"))
	evs <- as.numeric(readr::read_lines(paste0(ff, ".pca.ev")))
	attr(pcs, "explained") <- evs
	class(pcs) <- c("pca.result", class(pcs))
	return(pcs)
}

#' @export
merge_pca_result <- function(x, y, ...) {
	
	if (!all(inherits(x, "pca_result"), !inherits(y, "pca_result"), inherits(y, "data.frame"))) {
		stop("Can only merge a `pca.result` (x) with a `data.frame` (y).")
	}
	
	evs <- attr(x, "eigvals")
	explained <- attr(x, "explained")
	rez <- dplyr::left_join(x, y)
	class(rez) <- c("pca_result", class(rez))
	
	attr(rez, "explained") <- explained
	attr(rez, "eigvals") <- evs
	
	return(rez)
	
}

#' @export
melt_pc_pairs <- function(pc, K = 5, ...) {
	
	is.pc <- grepl("^PC\\d+", colnames(pc))
	pc.cols <- which(is.pc)
	v <- as.integer(gsub("^PC", "", colnames(pc)[ pc.cols ]))
	pc.cols <- pc.cols[ v <= K ]
	o <- order(v[ v <= K ])
	pc <- as.data.frame(pc)
	
	other.cols <- which(!is.pc)
	rez <- myldply(as.list(as.data.frame(combn(pc.cols, 2))), function(f) {
		this.pair <- pc[ ,other.cols ]
		this.pair$x <- pc[ ,f[1] ]
		this.pair$xvar <- colnames(pc)[ f[1] ]
		this.pair$y <- pc[ ,f[2] ]
		this.pair$yvar <- colnames(pc)[ f[2] ]
		return(tibble::as_tibble(this.pair))
	})
	
	rez$xvar <- factor(rez$xvar, colnames(pc)[ pc.cols[o] ])
	rez$yvar <- factor(rez$yvar, colnames(pc)[ pc.cols[o] ])
	return(rez)
	
}

#' @export
pc_label <- function(x, k = 1, ...) {
	pvar <- 100*attr(x, "explained")
	paste0("PC", k, " (", sprintf("%0.01f", pvar[k]), "%)")
}

#' @export
get_pca_palette <- function(f, palette = "Paired", ...) {
	
	if (is.factor(f))
		f <- factor(f)
	
	n <- nlevels(f)
	maxcol <- min(n, RColorBrewer::brewer.pal.info[ palette,"maxcolors" ])
	cluster.palette <- colorRampPalette(RColorBrewer::brewer.pal(maxcol, palette))
	return(cluster.palette)
	
}