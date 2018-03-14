
#' @export
.layout_agraph <- function(x,
						show_leaf_labels = TRUE,
						draw_leaves = TRUE,
						color = "yellowgreen",
						show_inner_node_labels = FALSE,
						draw_inner_nodes = FALSE,
						inner_node_color = color,
						show_admixture_labels = FALSE,
						parent_order = list(),
						child_order = list(),
						leaf_order = NULL,
						fix = list(),
						platform = 1,
						title = NULL,
						bounded_edges = NULL,
						...) {
	# Combine the user instructions and automated heuristics about the graph orderings.
	graph <- x
	arranged <- admixturegraph:::arrange_graph(graph)
	parent_order <- utils::modifyList(arranged$parent_order, parent_order)
	child_order <- utils::modifyList(arranged$child_order, child_order)
	if (is.null(leaf_order) == TRUE) {
		leaf_order <- admixturegraph:::leaf_order(graph, parent_order, child_order)
	} else if (typeof(leaf_order) != "character") {
		machine_order <- admixturegraph:::leaf_order(graph, parent_order, child_order)
		human_order <- leaf_order
		leaf_order <- rep("", length(leaf_order))
		for (i in seq(1, length(human_order))) {
			leaf_order[i] <- machine_order[human_order[i]]
		}
	}
	# Assign initial coordinates for all nodes.
	leaves <- list()
	if (length(leaf_order) > 1) {
		for (i in seq(1, length(leaf_order))) {
			leaves[[i]] <- c(100*(i - 1)/(length(leaf_order) - 1), 0)
		}
	} else {leaves[[1]] <- c(50, 0)}
	names(leaves) <- leaf_order
	parents <- graph$parents
	for (i in seq(1, length(graph$inner_nodes))) {
		candidate <- graph$inner_nodes[i]
		abandon <- FALSE
		for (j in seq(1, NCOL(parents))) {
			if (parents[candidate, j] == TRUE) {
				abandon <- TRUE
			}
		}
		if (abandon == FALSE) {
			root <- list(c(50, 100))
			names(root) <- c(candidate)
			delete <- i
		}
	}
	root_removed <- graph$inner_nodes[-delete]
	inner <- list()
	if (length(root_removed) > 0) {
		for (i in seq(1, length(root_removed))) {
			inner[[i]] <- c(0, 0)
		}
	}
	names(inner) <- root_removed
	# Assign the y-coordinates according to my arbitrary principles.
	refined_graph <- admixturegraph:::refined_graph(graph)
	heights <- rep(0, length(inner))
	names(heights) <- names(inner)
	global_longest <- 0
	for (inner_node in names(inner)) {
		paths <- admixturegraph:::all_paths_to_leaves(refined_graph, inner_node)
		longest <- 0
		for (path in paths) {
			if (length(path) > longest) {
				longest <- length(path)
			}
		}
		heights[inner_node] <- longest - 1
		if (longest > global_longest) {
			global_longest <- longest
		}
	}
	for (inner_node in names(inner)) {
		heights[inner_node] <- global_longest - heights[inner_node]
		inner[[inner_node]][2] <- 100*(1 - heights[inner_node]/global_longest)
	}
	# Perform Nelder-Mead to optimize the x-coordinates of the non-root inner nodes.
	if (length(inner) > 0) {
		x0 <- rep(50, length(inner))
		min <- rep(0, length(inner))
		max <- rep(100, length(inner))
		cfunc <- admixturegraph:::drawing_cost(graph, leaves, root, inner, child_order, parent_order, platform)
		opti <- suppressWarnings(neldermead::fminbnd(cfunc, x0 = x0, xmin = min, xmax = max))
		x <- neldermead::neldermead.get(opti, "xopt")
		for (i in seq(1, length(inner))) {
			inner[[i]][1] <- x[i]
		}
	}
	# Plot everything asked for.
	#xpd <- graphics::par()$xpd
	#graphics::par(xpd = NA)
	level <- platform*25/(max(2, length(leaves)) - 1)
	for (inner_node in names(inner)) {
		inner[[inner_node]][2] <- 100*(1 - heights[inner_node]/global_longest)
	}
	for (fixed in names(fix)) {
		inner[[fixed]] <- inner[[fixed]] + fix[[fixed]]
	}
	nodes <- graph$nodes
	coordinates <- c(leaves, root, inner)
	#graphics::plot(c(-level, 100 + level), c(0, 100), type = "n", axes = FALSE,
	#			   frame.plot = FALSE, xlab = "", ylab = "", main = title, ...)
	ee <- data.frame()
	labs <- data.frame()
	for (i in nodes) {
		for (j in nodes) {
			if (parents[i, j] == TRUE) {
				i_thing <- 0
				if (length(parent_order[[i]]) == 2) {
					if (parent_order[[i]][1] == j) {
						i_thing <- -level
						elab <- paste0("edge_", j, "_", i)
						dd <- data.frame(x0 = coordinates[[i]][1] + i_thing,
										 y0 = coordinates[[i]][2],
										 x1 = coordinates[[i]][1] - i_thing,
										 y1 = coordinates[[i]][2],
										 label = elab, stringsAsFactors = FALSE)
						#graphics::segments(coordinates[[i]][1] + i_thing, coordinates[[i]][2], coordinates[[i]][1] - i_thing,
						#				   coordinates[[i]][2], col = "black", lwd = 2)
						ee <- rbind(ee, dd)
					}
					if (parent_order[[i]][2] == j) {
						i_thing <- level
					}
					if (show_admixture_labels == TRUE) {
						label <- graph$probs[i, j]
						if (substr(label, 1, 1) == "(") {
							label <- substr(label, 2, nchar(label) - 1)
						}
						labs <- rbind(labs,
									  data.frame(x = coordinates[[i]][1] + 0.75*i_thing,
									  		   y = coordinates[[i]][2],
									  		   label = label, stringsAsFactors = FALSE)
						)
						#graphics::text(coordinates[[i]][1] + 0.75*i_thing, coordinates[[i]][2], label,
						#			   adj = c(0.5, 1.6), cex = 0.8)
					}
				}
				elab <- paste0("edge_", j, "_", i)
				dd <- data.frame(x0 = coordinates[[i]][1] + i_thing,
								 y0 = coordinates[[i]][2],
								 x1 = coordinates[[j]][1],
								 y1 =  coordinates[[j]][2],
								 label = elab, stringsAsFactors = FALSE)
				ee <- rbind(dd, ee)
				#graphics::segments(coordinates[[i]][1] + i_thing, coordinates[[i]][2], coordinates[[j]][1],
				#				   coordinates[[j]][2], col = "black", lwd = 2)
			}
		}
	}
	nn <- data.frame()
	for (i in seq(1, length(leaves))) {
		leaf <- leaves[[i]]
		if (draw_leaves == TRUE) {
			nn <- rbind(nn,
						data.frame(x = leaf[1], y = leaf[2], label = names(leaves)[i], type = "leaf", stringsAsFactors = FALSE)
						)
			#graphics::points(leaf[1], leaf[2], lwd = 2, pch = 21, col = "black", bg = color, cex = 2)
		}
		if (show_leaf_labels == TRUE) {
			#graphics::text(leaf[1], leaf[2], names(leaves)[i], adj = c(0.5, 2.6), cex = 0.8)
		}
	}
	if (length(inner) > 0) {
		for (i in seq(1, length(inner))) {
			vertex <- inner[[i]]
			if (draw_inner_nodes == TRUE) {
				nn <- rbind(nn,
							data.frame(x = vertex[1], y = vertex[2], label = names(inner)[i], type = "inner", stringsAsFactors = FALSE)
							)
				#graphics::points(vertex[1], vertex[2], lwd = 2, pch = 21, col = "black", bg = inner_node_color, cex = 2)
			}
			if (show_inner_node_labels == TRUE) {
				#graphics::text(vertex[1], vertex[2], names(inner)[i], adj = c(0.5, -1.6), cex = 0.8)
			}
		}
	}
	juuri <- root[[1]]
	if (draw_inner_nodes == TRUE) {
		nn <- rbind(nn,
					data.frame(x = juuri[1], y = juuri[2], label = names(root)[1], type = "inner", stringsAsFactors = FALSE)
					)
		#graphics::points(juuri[1], juuri[2], lwd = 2, pch = 21, col = "black", bg = inner_node_color, cex = 2)
	}
	if (show_inner_node_labels == TRUE) {
		#graphics::text(juuri[1], juuri[2], names(root)[1], adj = c(0.5, -1.6), cex = 0.8)
	}
	#graphics::par(xpd = xpd)
	
	
	ee$is.bound <- FALSE
	if (!is.null(bounded_edges)) {
		bounded_edges <- gsub(" \\= 0", "", bounded_edges)
		ee$is.bound <- ee$label %in% bounded_edges
	}
	
	return( structure(list(edges = tibble::as_tibble(ee),
						   nodes = tibble::as_tibble(nn)),
					  class = "agraph.layout") )
	
}

#' @export
.multilayout_agraph <- function(ll, fits = NULL, ...) {
	
	layouts <- lapply(seq_along(ll), function(f) .layout_agraph(ll[[f]], bounded_edges = fits[[f]]$bounded_edges), ...)
	names(layouts) <- names(ll)
	if (is.null(names(layouts)))
		names(layouts) <- paste0("m", seq_along(ll))
	nodes <- lapply(names(layouts), function(f) dplyr::mutate(layouts[[f]]$nodes, graph = f))
	edges <- lapply(names(layouts), function(f) dplyr::mutate(layouts[[f]]$edges, graph = f))
	nodes <- do.call(dplyr::bind_rows, nodes)
	edges <- do.call(dplyr::bind_rows, edges)
	return(list(nodes = nodes, edges = edges))
	
}

#' @export
plot_agraph <- function(g, bound.edges = NULL, ...) {
	
	if (inherits(g, "agraph.layout") || (is.list(g) & all(c("edges","nodes") %in% names(g))))
		layout <- g
	else if (inherits(g, "agraph"))
		layout <- .layout_agraph(g, ...)
	else
		stop("Unsupported input.")
		
	p <- ggplot2::ggplot(NULL) +
		ggplot2::geom_text(data = subset(layout$nodes, type == "leaf"),
				  ggplot2::aes(x = x, y = y, label = label),
				  angle = 90, hjust = 1) +
		ggplot2::geom_segment(data = layout$edges,
							  ggplot2::aes(x = x0, y = y0, xend = x1, yend = y1, colour = is.bound)) +
		ggplot2::scale_colour_manual(values = c("black","red"), guide = FALSE) +
		ggplot2::theme_void()
	
	return(p)
	
}

#' @export
.prep_fits <- function(fit, sigma = 6, ...) {
	
	D <- fit$D
	fit$stderr <- with(fit, D/Z.value)
	fit$error_bar_start <- with(fit, D - sigma/2 * stderr)
	fit$error_bar_end <- with(fit, D + sigma/2 * stderr)
	fit$test <- with(fit, paste("D(", W, ",", X, ";", Y, ",", Z, ")"))
	fit$test <- with(fit, reorder(factor(test), D))
	fit$hit <- with(fit, Vectorize(dplyr::between)(graph_f4, error_bar_start, error_bar_end))

	attr(fit, "prepped") <- TRUE
	return(fit)
	
}

#' @export
glance.agraph_fit <- function(fit, ...) {
	tibble::tibble(bounded_edges = length(fit$bounded_edges),
				   free_edges = length(fit$free_edges),
				   best_error = fit$best_error,
				   sse = fit$best_error,
				   complaints = length(fit$complaint),
				   nadmix = length(fit$best_fit))
}


#' @export
#' @method tidy agraph_fit
tidy.agraph_fit <- function(fit, ...) {
	.prep_fits(admixturegraph:::fitted.agraph_fit(fit))
}

#' @export
sort_f4stats <- function(fit, fitted = FALSE, ...) {
	if (fitted && all("graph_f4" %in% colnames(fit))) {
		fit$test <- with(fit, reorder(factor(test), graph_f4))
	}
	else {
		fit$test <- with(fit, reorder(factor(test), D))
	}
	return(fit)
		
}

#' @export
plot_fitted <- function(fit, ...) {
	
	if (inherits(fit, "agraph_fit")) {
		fit2 <- tidy.agraph_fit(.prep_fits(admixturegraph:::fitted.agraph_fit(fit)))
	}
	else if (!is.null(attr(fit, "prepped"))) {
		fit2 <- fit
	}
	else {
		fit2 <- fit
		warning("Don't understand input; going ahead blindly ...")
	}
	
	ggplot2::ggplot(fit2) + 
		ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
		ggplot2::geom_segment(ggplot2::aes(x = test, xend = test, y = D, yend = graph_f4), colour = "red", lty = "dashed") +
		ggplot2::geom_errorbar(ggplot2::aes(x = test, ymin = error_bar_start, ymax = error_bar_end)) + 
		ggplot2::geom_point(ggplot2::aes(x = test, y = D), shape = 3) +
		ggplot2::geom_point(ggplot2::aes(x = test, y = graph_f4, colour = hit), shape = 16) +
		ggplot2::scale_colour_manual(values = c("red","green"), guide = FALSE) +
		ggplot2::coord_flip()

}
