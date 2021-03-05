# functions for PCA and projection with `plink`

read_pca_plink <- function(ff, ...) {
	
	ev <- as.double(readLines(paste0(ff, ".eigenval")))
	npc <- length(ev)
	ct <- paste0(rep("d", npc), collapse = "")
	
	df <- readr::read_tsv(paste0(ff, ".eigenvec"), col_types = paste0("cc", ct))
	colnames(df) <- c("fid","iid", paste0("PC", seq_len(npc)))
	class(df) <- c("pca_result", class(df))
	
	attr(df, "eigvals") <- ev
	attr(df, "explained") <- ev/sum(ev)
	
	return(df)
	
}

#' Perform PCA using `plink2`
#' 
#' @param vcf path to `plink` binary fileset (bed/bim/fam)
#' @param skeleton optional character vector fo sample IDs to use for constructing PCs
#' @param samples optional character vector of sample IDs to project (if using `skeleton`, should include all reference samples)
#' @param max_missing input filter for `plink` corresponding to `--geno`
#' @param min_maf input filter for `plink` corresponding to `--maf`
#' @param ... ignored
#' 
#' @return An object of classes \code{tbl_df} and \code{pca_result} whose columns are the projections of
#' 	samples onto PCs.  Eigenvalues, expressed as proportion of variance explained, are returned as an attribute.
#' 	
#' @details This function generates command-line calls to \code{plink2}, whose output is directed to a set of
#' 	temporary files and then read into the \code{R} session. This requires the hottest version (alpha 3 or newer)
#' 	of \code{plink2}.
#' 	 	
#' @export
pca_plink <- function(bfile, skeleton = NULL, samples = NULL, max_missing = 0, min_maf = 0, ...) {
	
	ff <- tempfile()
	keeps <- paste0(ff, ".ref")
	targets <- paste0(ff, ".target")
	
	if (!is.null(samples)) {
		samples <- as.character(samples)
		writeLines(samples, targets)
	}
	
	if (!is.null(skeleton)) {
		
		skeleton <- as.character(skeleton)
		writeLines(skeleton, keeps)
		
		cmd1 <- paste0("plink2",
					   " --bfile ", bfile,
					   " --keep ", keeps,
					   " --geno ", max_missing,
					   " --maf ", min_maf,
					   " --freq counts",
					   " --pca allele-wts",
					   " --out ", ff, ".ref_pcs")
		message("Constructing PC axes with reference samples ...")
		message("The command will be: ", cmd1)
		rez <- system(cmd1)
		
		cmd2 <- paste0("plink2 --bfile ", bfile)
		if (!is.null(samples)) {
			cmd2 <- paste0(cmd2, " --keep ", targets)
		}
		cmd2 <- paste0(cmd2,
					   " --read-freq ", ff, ".ref_pcs.acount",
					   " --score ", ff, ".ref_pcs.eigenvec.allele 2 5 header-read",
					   " no-mean-imputation variance-standardize",
					   " --score-col-nums 6-15",
					   " --out ", ff, ".proj")
		
		message("Projecting onto reference PCs ...")
		message("The command will be: ", cmd2)
		rez <- system(cmd2)
		
		df <- readr::read_tsv(paste0(ff, ".proj.sscore"), col_types = c("ccdddddddddddd"))
		colnames(df) <- c("fid","iid","allele_count","sum_dosage", paste0("PC", 1:10))
		class(df) <- c("pca_result", class(df))
		
		ev <- as.double(readLines(paste0(ff, ".ref_pcs.eigenval")))
		attr(df, "eigvals") <- ev
		attr(df, "explained") <- ev/sum(ev)
		
		return(df)
		
	}
	
	else {
		
		cmd2 <- paste0("plink2 --bfile ", bfile)
		if (!is.null(samples)) {
			cmd2 <- paste0(cmd2, " --keep ", targets)
		}
		cmd2 <- paste0(cmd2,
					   " --geno ", max_missing,
					   " --maf ", min_maf,
					   " --pca ",
					   " --out ", ff)
		
		message("Constructing PCs and projecting in same sample set ...")
		message("The command will be: ", cmd2)
		rez <- system(cmd2)
		
		return( read_pca_plink(ff) )
		
	}
	
}

