avg_matrix = function(mat, lengths){
	require(Matrix)
  	
	if(sum(lengths) != nrow(mat)){
		stop("Lengths does not sum up to number of rows in mat")
	}
	
	diag = spMatrix(nrow = length(lengths), ncol = sum(lengths), i = rep(1:length(lengths), lengths), j = 1:sum(lengths), x = 1 / rep(lengths, lengths))
	
	return(diag %*% mat)
}

## Functions for calculating the description length
row_subsets = function(mat, k){
	n = ncol(mat)
	m = nrow(mat)
	
	res = matrix(NA, n * k, m - k + 1)
	for(i in 1:(m - k + 1)){
		res[, i] = as.vector(t(mat[i:(i + k - 1), ]))
	}
	return(res)
}

row_model_spec = function(model, annotation, k, individual = NA, lme = FALSE){
	rownames(annotation) = NULL
	annotation = do.call("rbind", replicate(k, annotation, simplify = FALSE))
	
	if(lme){
		# annotation = cbind(annotation, "individual" = factor(rep(individual, k)))
		res = list(
			model = paste("y ", model, sep = ""),
			annotation = annotation
		)
	}
	else{
		mmformula = formula(model)
		mm = model.matrix(mmformula, annotation)
		res = list(
			mm = mm
		)
	}
	
	return(res)
}

vars = function(rss, n, df, n0, sig0){
	s = ((n0 * sig0 + n * rss / df) / (n0 + n))
	v = n0 + n
	
	var = (v * s) / (v - 2)
	
	return(var)
}

calculate_rss = function(fit){
	if (is.matrix(fit$residuals)){
		rss = apply(fit$residuals, 2, function(r) sum(r * r))
	}
	else{
		rss = sum(fit$residuals * fit$residuals)
	}
	return(rss)
}

number_of_coefs = function(fit){
	if (is.matrix(fit$residuals)){
		n = nrow(fit$residuals)
	}
	else{
		n = length(fit$residuals)
	}
	return(n)
}

row_bics = function(fit, n0,  m0, sig0, alpha){
	k = fit$rank
	rss = calculate_rss(fit)
	n = number_of_coefs(fit)
	vars = vars(rss, n, fit$df.residual, n0, sig0)

	mdl = ifelse(vars>0, 0.5 * (n*(log2(2*pi) + log2(vars)) + 1/vars/log(2)*rss + alpha*k*log2(m0)), 0.5*alpha*k*log2(m0))
	return(mdl)
}

row_fit = function(mat, fit_par, k, return_model = FALSE){
	y = row_subsets(mat, k)
	
	if(return_model){
		# temp = cbind(y = y, fit_par$annotation)
		# fit = tryCatch({
			# lme(formula(fit_par$model), random = ~ 1|individual, data = temp)
		# }, error = function(x){return(NULL)})
		fit = lm(formula(fit_par$model))
		return(fit)
	}
	else{
		fit = lm.fit(fit_par$mm, y)
	}
	return(fit)
}

prior_var = function(values, model, annotation){
	fp = row_model_spec(model, annotation, k = 1, lme = F)
	f = row_fit(values, fp, k = 1)
	sig0 = mean(calculate_rss(f) / f$df.residual)
	
	return(sig0)
}

description_length_lm = function(values, maxlength, model, annotation, n0, m0, sig0, alpha, max_block_length){
	n = ncol(values)
	m = nrow(values)
	
	# Estimate variance
	if(is.na(sig0)){
		if(m > 100){
			sig0 = prior_var(values[sample(1:m, 100), ], model, annotation)
		}
		else{
			sig0 = prior_var(values, model, annotation)
		}
	}
	
	n0 = n * n0
	m0 = n * m0
	
	# Calculate scores
	bics = matrix(NA, m, m)
	m2 = min(m, max_block_length)
	for(i in 1:m2){
		a = rep(NA, (m - i + 1))

		fp = row_model_spec(model, annotation, k = i, lme = FALSE)
		a = row_bics(row_fit(values, fp, i, return_model = FALSE), n0, m0, sig0, alpha)
		
		a = a + alpha*log2(maxlength)
		
		bics[cbind(1 : (m - i + 1), (1 : (m - i + 1)) + i - 1)] = a
	}
	return (list("bics" = bics, "max_block_length" = max_block_length))
}
##

## Functions for annotating the segments
# Calculate pvalues and lm coefficients for all regions
 
#' Fits contrasts on segments
#' 
#' Fits contrasts on segments
#' 
#' Contrasts have to be specified as a list where the element name shows the factor and the value is the character 
#' representation of the contrasts, using the levels of the variable check \code{\link{glht}} on how to specify contrasts
#'
#' @param seqlmresults a seqlmresults object
#' @param model_nr model to consider for the contrasts
#' @param contr list with contrasts
#' @return  A list with one segmentation table, with coefficients and p-values, for every contrast
#' 
#' @author  Kaspar Martens <kmartens@@ut.ee> Raivo Kolde <rkolde@@gmail.com>
#' 
#' 
#' @export

seqlm.contrasts = function(seqlmresults){
	coord = seqlmresults$data$genome_information$pos
	chr = seqlmresults$data$genome_information$chr

	segments = seqlmresults$segments

	if(nrow(segments) == 0){
		return(NULL)
	}

	model = seqlmresults$description_length_par$model
	annotation = seqlmresults$description_length_par$annotation
	values = seqlmresults$data$values
	individual = colnames(values)
	
	# contr.names = rep(names(contr), sapply(contr, length))
	
	avg_mat = avg_matrix(values, segments$length)
	fit_par = row_model_spec(model, annotation, k = 1, individual = individual, lme = FALSE) 
	m = row_fit(avg_mat, fit_par, 1, return_model = FALSE)
	
	calc.lm <- function(x){
		Qr <- x$qr
		p <- x$rank
		p1 <- 1L:p
		rss <- colSums(x$residuals * x$residuals)

		n <- NROW(Qr$qr)
		rdf <- n - p

		resvar <- rss/rdf
		R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
		se <- sqrt(diag(R) * resvar)

		est <- x$coefficients[2, ]
		tval <- est/se

		res <- cbind(coef = est, se = se, tstat = tval, p.value = 2 * pt(abs(tval), df=rdf, lower.tail=FALSE))
		res
	}
	
	output = calc.lm(m)

	return(output)
	
	# pvalues = foreach(i = 1:nrow(segments), .export=c("glht", "row_fit", "row_model_spec", "row_subsets", "mcp", "pmvt")) %dopar% {
		# start = segments$startIndex[i]
		# end = segments$endIndex[i]
		# mat = apply(values[(start):(end), , drop=FALSE], 2, mean)
		# fit_par = row_model_spec(model, annotation, k = 1, individual = individual, lme = TRUE) 
		
		# m = row_fit(mat, fit_par, nrow(mat), return_model = TRUE)
		# if(is.null(m)){
			# res = matrix(rep(c(NA, NA, NA), length(contr.names)), ncol=3)
		# }
		# else{
			# glht.m = glht(m, linfct = do.call("mcp", contr))
			# glht.m$vcov[glht.m$vcov == 0] = .Machine$double.eps
			# testsummary = summary(glht.m)
			# test = testsummary$test
			# res = matrix(c(test$coefficients, test$tstat, test$pvalues), ncol=3)
		# }
		# res
	# }
	
	# output = foreach(i = 1:nrow(pvalues[[1]])) %do% {
		# pieces = lapply(pvalues, function(x){c(x[i, 1], tstat = x[i, 2], p.value = x[i, 3])})
		# res = do.call("rbind", pieces)
		# colnames(res) = c("coef", "tstat", "p.value")
		# contrast = list(unlist(contr)[i])
		# names(contrast) = contr.names[i]
		# attr(res, "contrast") = contrast
		# res
	# }
	
	
	# additionalAnnotation = elementMetadata(seqlmresults$data$genome_information)
	# if(ncol(additionalAnnotation) == 0){
		# output = foreach(x = output) %do% {
			# gr = GRanges(seqnames = segments$chr, ranges = IRanges(start = segments$startPos, end = segments$endPos))
			# df = cbind(segments[, c(-1, -2, -3)], x)
			# elementMetadata(gr) = as(df, "DataFrame")
			# attr(gr, "contrast") = attr(x, "contrast")
			# attr(gr, "model") = model
			# gr[order(abs(df[, "tstat"]), decreasing=TRUE), ]
		# }
	# }
	# else{
		# # Add annotation
		# contr.ann = additional_annotation(segments$startIndex, segments$endIndex, additionalAnnotation)
			
		# output = foreach(x = output) %do% {
			# gr = GRanges(seqnames = segments$chr, ranges = IRanges(start = segments$startPos, end = segments$endPos))
			# df = cbind(segments[, c(-1, -2, -3)], x, contr.ann)
			# elementMetadata(gr) = as(df, "DataFrame")
			# attr(gr, "contrast") = attr(x, "contrast")
			# attr(gr, "model") = model
			# gr[order(abs(df[, "tstat"]), decreasing=TRUE), ]
		# }
	# }
	
	# # names(output) = pvalues[[1]]$names
	# return (output)
}



additional_annotation = function(startIndexes, endIndexes, df, variables = names(df)){
	output = list()
	cat("Adding annotation to the segments.\n")
	for (j in 1:length(variables)){
		variable = variables[j]
		currentData = df[, variable]
		cat(sprintf("\tVariable: %s\n", variable))
		if(is.numeric(currentData)){
			output[[j]] = foreach(i = 1:length(startIndexes), .combine='c') %dopar% {
				start = startIndexes[i]
				end = endIndexes[i]
				sum(currentData[(start):(end)])
			}
		}
		else{
			output[[j]] = foreach(i = 1:length(startIndexes), .combine='c') %dopar% {
				start = startIndexes[i]
				end = endIndexes[i]
				temp = unique(unlist(strsplit(paste(currentData[(start):(end)], collapse=";"), ";")))
				paste(temp[temp!=""], collapse=";")
			}
		}
		names(output)[j] = variables[j]
	}
	return(as.data.frame(output))
}

##

## Function to call seqfit for the linear model
 
#' Sequential lm
#' 
#' Segments genome based on given linear models
#' 
#' Implementation details
#'
#' @param values a matrix where columns are samples and rows correspond to the sites
#' @param genome_information \code{data.frame} giving the genomic coordinates and optionally additional 
#' description for the sites.
#' @param models vector of models given as character representation of formulas. The formulas should use 
#' variables present in the annotation table.
#' @param annotation table with sample annotations that are used to specify models
#' @param n0 prior number of observations to stablilize the variation estimate when calculating likelihood
#' @param sig0 prior standard deviation
#' @param max_block_length maximal length of the block we are searching. This is used to speed up computation
#' @param max_block_length_second_stage maximal block length for the second stage of search if two-stage search is used 
#' for speeding up the analysis. Second stage is initiated only if \code{max_block_length_second_stage} > 
#' \code{max_block_length}
#' @param max_dist maximal genomic distance between the sites to be considered the same region
#' @return  A list containing the input data, parameters and the segmentation.
#' 
#' @author  Kaspar Martens <kmartens@@ut.ee> Raivo Kolde <rkolde@@gmail.com>
#' 
#' @examples
#' # library(pheatmap)
#' # Generate data
#' rmat = function(n, m, mean=0, sd=1){ return(matrix(rnorm(n * m, mean, sd), n, m))}
#' sd = 0.1
#' mat = cbind(rmat(10, 3, 0.3, sd),
#' 	rbind(rmat(3, 5, 0.8, sd), rmat(7, 5, 0.2, sd)),
#' 	rbind(rmat(5, 8, 0.2, sd), rmat(5, 8, 0.9, sd))
#' )
#' colnames(mat) = paste("X", 1:ncol(mat), sep = "")
#' 
#' annotation = data.frame(Factor1 = rep(c("A", "B"), c(3, 7)), Factor2 = rep(c("A", "B"), c(5, 5)))
#' 
#' # pheatmap(mat, cluster_rows=F, cluster_cols=F)
#' 
#' # Generate necessary files for running the algorithm
#' genome_information = as.data.frame(list("chr" = rep(1, 16), "pos" = 1:16, "name" = sample(letters, 16)))
#' rownames(genome_information) = colnames(mat)
#' 
#' # Specify models
#' models = c("~ 1", "~ Factor1 + Factor2")
#' 
#' # Calculate segmentation
#' seqlmresults = seqlm(t(mat), genome_information, models, annotation, max_dist=2)
#' 
#' # Use contrasts to calculate p-values and effect sizes
#' contr.res = seqlm.contrasts(seqlmresults, model_nr = 2, contr = list(Factor1 = "A - B = 0", Factor2 = "A - B = 0"))
#' 
#' # Draw figures of the segments
#' seqlmplots(contr.res[[1]], seqlmresults, annotation_column = "Factor1", coef = 0, length = 0, adjusted.pvalue = 0.05)
#' 
#' # Generate a report about these figures
#' dir = tempdir()
#' seqlmreport(contr.res[[1]], seqlmresults, annotation_column = "Factor1", coef = 0, length = 0, adjusted.pvalue = 0.05, 
#' dir = dir)
#' 
#' @export
seqlm = function(values, genome_information, model, annotation, n0 = 1, m0 = 10, sig0 = NA, alpha = 2, max_block_length = 50, max_dist = 4000){
	
	return (seqfit(values, genome_information, max_dist, description_length_fun = description_length_lm, max_block_length = max_block_length, description_length_par = list(model = model, annotation = annotation,  n0 = n0, m0 = m0, sig0 = sig0, alpha = alpha)))
}
##

## Functions to visualize the regions
fortify_seqlmplot = function(segment, values, sample_annotation, genome_information, genomic_coordinates, expand){
	# Get rows from the matrix
	which_rows = (segment$start_which):(segment$end_which)
	which_rows_expand = max(0, segment$start_which - expand):min(nrow(values), segment$end_which + expand)
	values0 = as.data.frame(t(values[which_rows_expand, ]))
	
	# Compile with annotations
	df0 = cbind(values0, group = sample_annotation)
	
	# Bring into long format
	df = melt(df0, id.vars=c("group"))
	
	# Add coordinates
	if(genomic_coordinates){
		df$variable = rep(genome_information$pos[which_rows_expand], each=ncol(values))
	}
	
	# Calculate group means to show lines on the figure
	group.means = ddply(df, .(group, variable), summarise, value.mean = mean(value))
	
	return(list(df = df, group.means = group.means))
}

draw_seqlmplot = function(df, group.means, ylim){
	 plot = ggplot(data=df, aes(x=variable, y=value, color=group)) + geom_jitter(position = position_jitter(width = .1)) + scale_y_continuous(limits = ylim) +
	theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
	geom_path(data=group.means, aes(x=variable, y=value.mean, group=group, colour=group), size=1) +
	geom_point(data=group.means, aes(x=variable, y=value.mean, colour=group), size=4, shape=21, fill="white")
}
	
seqlmplot = function(segment, values, sample_annotation, genome_information, genomic_coordinates = F, expand = 0, ylim = extendrange(values), filename = NA, ...){
	data = fortify_seqlmplot(segment = segment, values = values, sample_annotation = sample_annotation, genome_information = genome_information, genomic_coordinates = genomic_coordinates, expand = expand)

	plot = draw_seqlmplot(df = data$df, group.means = data$group.means, ylim)
	
	if(is.na(filename)){
		print(plot)
	}
	else{
		ggsave(filename, plot, ...)
	}
}
	
 
#' Plot the significant segments that 
#' 
#' Function that draws the plots of significant segments. 
#' 
#' See examples in \code{\link{seqlm}}
#'
#' @param contr.res a table corresponding one contrast from the results of \code{\link{seqlm.contrasts}}
#' @param seqlmresults original seqlm results file
#' @param length minimal length of region to draw
#' @param coefficient minimal absolute coefficient value
#' @param adjusted.pvalue threshold for adjusted pvalues
#' @param annotation_column the figure is the column in the sample annotation file that we use to group the data
#' @param genomic_coordinates boolean determining if the sites are drawn equidistantly or based on genomic coordinates
#' @param expand number of extra probes to show on both sides of the region
#' @param ylim two element vector giving the lower and higher limit of the y axis
#' @param dir directory where to put the images, of NA then plots are drawn into the plotting window
#' @param filetype picture filetype 
#' @param ... extra parameters to \code{\link{ggsave}}
#' 
#' @author  Kaspar Martens <kmartens@@ut.ee> Raivo Kolde <rkolde@@gmail.com>
#' 
#' 
#' 
#' @export
seqlmplots = function(contr.res, seqlmresults, length = 10, coefficient = 0.2, adjusted.pvalue = 0.05, annotation_column, genomic_coordinates = F, expand = 0, ylim = extendrange(seqlmresults$data$values), dir = NA, filetype = "png", ...){
	# Select only relevant regions
	contr.res = subset(contr.res, lengths >= length & abs(coef) >= coefficient & adjusted.pvalues <= adjusted.pvalue)
	
	for(i in 1:nrow(contr.res)){
		if(is.na(dir)){
			filename = NA
		}
		else{
			filename = file.path(dir, sprintf("%d.%s", i, filetype))
		}
	
		seqlmplot(segment = contr.res[i, ], values = seqlmresults$data$values, sample_annotation = seqlmresults$description_length_par$annotation[, annotation_column], genome_information = seqlmresults$data$genome_information, genomic_coordinates = genomic_coordinates, expand = expand, ylim = ylim, filename = filename, ...)
	}
}
	
# seqlmplots(pval[[1]], seqlmresults, annotation_column = "Factor2", coef = 0, length = 0, dir = "~/Desktop/")
##

## seqlm raport
raport_template = '
<!DOCTYPE html>
<html>
<head>
<style type="text/css">.knitr.inline {
	background-color: #f7f7f7;
	border: solid 0px #b0b0b0
}
.message {
	font-style: italic
}
.source,.output,.warning,.error,.message {
	padding: 0em 1em;
	border: solid 1px #f7f7f7
}
.source {
	background-color: #f7f7f7
}
.rimage.left {
	text-align: left
}
.rimage.right {
	text-align: right
}
.rimage.center {
	text-align: center
}
.source {
	color: #333
}
.background {
	color: #f7f7f7
}
</style>
<title>%s</title>
</head>
<body>

<code class="knitr inline">
<h1> %s </h1>

%s
</code>
</body>
</html>
'
chunk_template = '
<h2> Segment %d </h2> 

<table>
	<tr>
		<td><b>Location</b></td>
		<td>%s</td>
	</tr>
	%s
</table>

<div class="rimage default"><img src="%s" class="plot"/></div>
'

annotation_template = '
<tr>
	<td><b>%s</b></td>
	<td>%s</td>
</tr>
'

location_template = 'chr%s:%d-%d'

annotation_table = function(x){
	n = which(colnames(x) == "adjusted.pvalues")
	
	res = paste(sprintf(annotation_template, "Coefficient", round(x[1, "coef"], 3)), sprintf(annotation_template, "Adjusted p-value", x[1, "adjusted.pvalues"]), sprintf(annotation_template, "No. probes", x[1, "lengths"]), sprintf(annotation_template, "Length in bp", x[1, "endPos"] - x[1, "startPos"]))
	
	if(!(n == ncol(x))){
		for(i in (n + 1):ncol(x)){
			res = paste(res, sprintf(annotation_template, colnames(x)[i], x[1, i]), sep = "\n")
		}
	}
	
	return(res)
}

 
#' Generate the HTML report for the seqlm results
#' 
#' Generate the HTML report for the seqlm results
#'
#' @param contr.res a table corresponding one contrast from the results of \code{\link{seqlm.contrasts}}
#' @param seqlmresults original seqlm results file
#' @param length minimal length of region to draw
#' @param coefficient minimal absolute coefficient value
#' @param adjusted.pvalue threshold for adjusted pvalues
#' @param orderby name of the contrast results table column to order results by
#' @param annotation_column the figure is the column in the sample annotation file that we use to group the data
#' @param genomic_coordinates boolean determining if the sites are drawn equidistantly or based on genomic coordinates
#' @param expand number of extra probes to show on both sides of the region
#' @param ylim two element vector giving the lower and higher limit of the y axis
#' @param dir directory where to put the page, if the directory does not exist it will be created
#' @param width picture width in inches
#' @param height picture height in inches
#' @param dpi dots per inch, to calibrate the picture size in pixels
#' @param main title for the raport
#' 
#' @author  Kaspar Martens <kmartens@@ut.ee> Raivo Kolde <rkolde@@gmail.com>
#' 
#' @export
seqlmreport = function(contr.res, seqlmresults, length = 10, coefficient = 0.2, adjusted.pvalue = 0.05, orderby = "coef", annotation_column, genomic_coordinates = F, expand = 0, ylim = extendrange(seqlmresults$data$values), dir = NA, width = 8, height = 5, dpi = 100, main = "seqlm results"){
	# Select only relevant regions
	contr.res = subset(contr.res, lengths >= length & abs(coef) >= coefficient & adjusted.pvalues <= adjusted.pvalue)
	contr.res = contr.res[order(contr.res[, orderby]), ]
	
	# Create main directory 
	if(!file.exists(dir)){
		dir.create(dir)
	}
	
	# Create image directory
	img_dir = file.path(dir, "img")
	if(!file.exists(img_dir)){
		dir.create(img_dir)
	}
	
	# Create images
	seqlmplots(contr.res, seqlmresults, length = 0, coefficient = 0, adjusted.pvalue = 1, annotation_column = annotation_column, genomic_coordinates = genomic_coordinates, expand = expand, ylim = ylim, dir = img_dir, filetype = "png", width = width, height = height, dpi = dpi)
	
	# Create HTML file
	chunks = ''
	
	for(i in 1:nrow(contr.res)){
		location = sprintf(location_template, contr.res[i, "chr"], contr.res[i, "startPos"], contr.res[i, "endPos"])
		
		chunk = sprintf(chunk_template, i, location, annotation_table(contr.res[i, ]), sprintf("img/%d.png", i))
		
		chunks = paste(chunks, chunk, sep = "\n\n")
	}
	
	page = sprintf(raport_template, main, main, chunks)
	
	cat(page, file = file.path(dir, "index.html"))
}


##
