## Functions to visualize the regions
fortify_seqlmplot = function(segment, values, annotation, genome_information, expand){
	# Get rows from the matrix
	segment_expanded = segment
	start(segment_expanded) = start(segment_expanded) - expand
	end(segment_expanded) = end(segment_expanded) + expand

	gi = subsetByOverlaps(genome_information, segment)
	gi_expanded = subsetByOverlaps(genome_information, segment_expanded)
	
	values0 = as.matrix(values[names(gi_expanded), ])
	
	# Bring into long format
	df = melt(values0, varnames = c("Probe", "Sample"))
	
	# Add annotations
	a = data.frame(Sample = colnames(values), Annotation = annotation)
	df = merge(df, a)
	
	# Add region information
	probe_det = data.frame(Probe = names(gi_expanded), Region = gi_expanded %over% gi, Position = start(gi_expanded)) 
	df = merge(df, probe_det)
	
	# Calculate the box that shows region
	reg = probe_det[probe_det$Region,]
	nonreg = probe_det[!probe_det$Region,]
	
	smaller = nonreg[nonreg$Position < min(reg$Position),]
	bigger = nonreg[nonreg$Position > max(reg$Position),]
	
	if(nrow(smaller) > 0){
		diff = min(reg$Position) - max(smaller$Position)
		start = min(reg$Position) - min(100, diff / 2)
	}
	else{
		start = min(reg$Position) - 100
	}
	
	if(nrow(bigger) > 0){
		diff = min(bigger$Position) - max(reg$Position)
		end = max(reg$Position) + min(100, diff / 2)
	}
	else{
		end = max(reg$Position) + 100
	}
	
	box = data.frame(start = start, end = end)
	
	return(list(df = df, box = box))
}

draw_seqlmplot = function(df, box, ylim, expand){
	plot = qplot(x = Position, y = value, geom = c("line", "point"), colour = Annotation, group = Sample, data = df) + geom_rect(aes(xmin = box$start, xmax = box$end, ymin = -Inf, ymax = Inf), colour = "grey20", fill = "grey95") + geom_point() + geom_line() +  geom_jitter(position = position_jitter(width = .1)) + scale_y_continuous(limits = ylim) + scale_x_continuous(limits = c(box$start - expand, box$end + expand)) + theme_bw() 
}
	
seqlmplot = function(segment, values, annotation, genome_information, expand, ylim = extendrange(values), filename = NA, ...){
	data = fortify_seqlmplot(segment = segment, values = values, annotation = annotation, genome_information = genome_information, expand = expand)
	
	plot = draw_seqlmplot(df = data$df, box = data$box, ylim, expand = expand)
	
	if(is.na(filename)){
		print(plot)
	}
	else{
		ggsave(filename, plot, ...)
	}
}
	
 
#' Visualise the regions
#' 
#' Generate plots about the seqlm results
#' 
#' The number of results from \code{\link{seqlm}} can be large 
#' and visualising all these regions might not be desirable. 
#' Therefore, it is advisable to filter the results befor 
#' plotting.  
#'
#' @param segments selection of significant regions by \code{\link{seqlm}} function 
#' @param values same values matrix that was used in \code{seqlm}
#' @param genome_information same genome_information object that was used in \code{seqlm}
#' @param annotation same annotation vector that was used in \code{seqlm}
#' @param expand number of basepairs to extend the region on plot
#' @param ylim two element vector giving the lower and higher limit of the y axis
#' @param dir  directory where to put the images, of NA then plots are drawn into the plotting window
#' @param filetype picture filetype 
#' @param  ... extra parameters to \code{\link{ggsave}}
#' @author  Raivo Kolde <rkolde@@gmail.com>
#' 
#' @export
seqlmplots = function(segments, values, genome_information, annotation, expand = 100, ylim = extendrange(values), dir = NA, filetype = "png", ...){
	# Match values and genome_information
	mp = match_positions(values, genome_information)
	values = mp$values
	genome_information = mp$genome_information
	
	# Draw pictures
	for(i in 1:length(segments)){
		if(is.na(dir)){
			filename = NA
		}
		else{
			filename = file.path(dir, sprintf("%d.%s", i, filetype))
		}
		seqlmplot(segment = segments[i], values = values, annotation = annotation, genome_information = genome_information,  expand = expand, ylim = ylim, filename = filename, ...)
		
	}
}
	

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
	res = paste(sprintf(annotation_template, "Coefficient", round(x$coef, 3)), 
		sprintf(annotation_template, "FDR", sprintf("%.3g", x$fdr)),
		sprintf(annotation_template, "Bonferroni", sprintf("%.3g", x$bonferroni)), 
		sprintf(annotation_template, "Length in probes", x$length),
		sprintf(annotation_template, "Length in bp", end(x) - start(x))
	)
	
	xx = as.data.frame(elementMetadata(x))
	n = which(colnames(xx) == "bonferroni")
	
	if(!(n == ncol(xx))){
		for(i in (n + 1):ncol(xx)){
			res = paste(res, sprintf(annotation_template, colnames(xx)[i], xx[1, i]), sep = "\n")
		}
	}
	
	return(res)
}

 
#' Generate the HTML report for the seqlm results
#' 
#' Generate the HTML report for the seqlm results
#'
#' @param segments selection of significant regions by \code{\link{seqlm}} function 
#' @param values same values matrix that was used in \code{seqlm}
#' @param genome_information same genome_information object that was used in \code{seqlm}
#' @param annotation same annotation vector that was used in \code{seqlm}
#' @param expand number of basepairs to extend the region on plot
#' @param ylim two element vector giving the lower and higher limit of the y axis
#' @param dir directory where to put the page, if the directory does not exist it will be created
#' @param width picture width in inches
#' @param height picture height in inches
#' @param dpi dots per inch, to calibrate the picture size in pixels
#' @param main title for the report
#' 
#' @author  Kaspar Martens <kmartens@@ut.ee> Raivo Kolde <rkolde@@gmail.com>
#' 
#' @export
seqlmreport = function(segments, values, genome_information, annotation, ylim = extendrange(values), dir = NA, expand = 100, width = 8, height = 5, dpi = 100, main = "seqlm results"){
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
	seqlmplots(segments, values, genome_information, annotation, ylim = ylim, dir = img_dir, expand = expand, width = width, height = height, dpi = dpi)
	
	# Create HTML file
	chunks = ''
	
	for(i in 1:length(segments)){
		location = sprintf(location_template, seqnames(segments[i]), start(segments[i]), end(segments[i]))
		
		chunk = sprintf(chunk_template, i, location, annotation_table(segments[i]), sprintf("img/%d.png", i))
		
		chunks = paste(chunks, chunk, sep = "\n\n")
	}
	
	page = sprintf(raport_template, main, main, chunks)
	
	cat(page, file = file.path(dir, "index.html"))
}


##
