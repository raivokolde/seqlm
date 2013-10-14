# ## Functions to visualize the regions
# 
# plotTrackList = function(extended.segment, segment, model, modelAnnotation, contr, values, sampleAnnotation, genomeAnnotation, geneAnnotation = NULL, cpgAnnotation = NULL, addDataTrack = FALSE, showRegion = FALSE, centeredData = FALSE, extend.nonzero = FALSE, ylim = extendrange(values)){
# 	chr = as.character(segment@seqnames@values)
# 	
# 	plotTrackList = list()
# 	if(!(is.null(geneAnnotation))){
# 		which.gene.overlap = findOverlaps(geneAnnotation, extended.segment)@queryHits
# 		gene.overlap = geneAnnotation[which.gene.overlap]
# 		geneTrack = GeneRegionTrack(rstarts = gene.overlap$exonStarts, rends = gene.overlap$exonEnds, symbol = gene.overlap$symbol, chromosome = chr, showId = TRUE, name="UCSC gene", cex=1.5)
# 		plotTrackList = c(plotTrackList, "geneTrack" = geneTrack)
# 	}
# 	if(!(is.null(cpgAnnotation))){
# 		which.cpg.overlap = findOverlaps(cpgAnnotation, extended.segment)@queryHits
# 		cpg.overlap = cpgAnnotation[which.cpg.overlap, ]
# 		cpgTrack = AnnotationTrack(cpg.overlap, name="CpG", stacking = "dense", Shelf = "#C6DBEF", Shore = "#4292C6", Island = "#08306B", showId = TRUE)
# 		feature(cpgTrack) = cpg.overlap$type
# 		plotTrackList = c(plotTrackList, "cpgTrack" = cpgTrack)
# 	}
# 	
# 	if(showRegion){
# 		regionTrack = AnnotationTrack(segment, name="DMR", fill = "#7f7f7f")
# 		plotTrackList = c(plotTrackList, "regionTrack" = regionTrack)
# 	}
# 	
# 	if(extend.nonzero){
# 		which.rows = findOverlaps(genomeAnnotation, extended.segment)@queryHits
# 	}
# 	else{
# 		which.rows = (segment$startIndex):(segment$endIndex)
# 	}
# 	genomeAnnotation = genomeAnnotation[which.rows, ]
# 	values = values[which.rows, ]
# 	if(centeredData) values = as.data.frame(t(apply(values, 1, function(x){x - mean(x)})))
# 	
# 	axisTrack = GenomeAxisTrack(cex=1.5, labelPos="below")
# 
# 	coefs = apply(values, 1, function(x){
# 		x = matrix(x, nrow=1)
# 		m = row_fit(x, model, modelAnnotation, 1, n0 = NA, m0 = NA, sig0 = NA, return_fit = TRUE)
# 		glht.m = glht(m, linfct = do.call("mcp", contr))
# 		glht.m$vcov[glht.m$vcov == 0] = .Machine$double.eps
# 		coef = summary(glht.m)$test$coefficients
# 		return(coef)
# 	})
# 	coefTrack = DataTrack(range = genomeAnnotation@ranges, data = coefs, chromosome = chr, genome="hg19", name="Slope", type=c("p", "polygon"), pch=16, ylim = c(-max(abs(coefs)), max(abs(coefs))))
# 	plotTrackList = c(plotTrackList, "coefTrack" = coefTrack)
# 
# 	if(addDataTrack){
# 		dataTrack = DataTrack(range = genomeAnnotation@ranges, data = values, chromosome = chr, genome="hg19", name="Beta values", groups = sampleAnnotation, legend=TRUE, pch=16, lwd=1, type=c("a", "p"), ylim = ylim)
# 		plotTrackList = c(plotTrackList, "dataTrack" = dataTrack)	
# 	}
# 	
# 	plotTrackList = c(plotTrackList, "axisTrack" = axisTrack)
# 	
# 	return(plotTrackList)
# }
# 	
# seqlmplot = function(segment, model, modelAnnotation, contr, values, sampleAnnotation, genomeAnnotation, geneAnnotation = NULL, cpgAnnotation = NULL, addDataTrack = FALSE, centeredData = FALSE, extend.left = 0, extend.right = 0, ylim = extendrange(values), filename = NA, graphics_device = png, width = 480, height = 480, units = "px", pointsize = 12, ...){
# 	# Find new start and end
# 	from = start(segment) - extend.left
# 	to = end(segment) + extend.right
# 	extended.segment = segment
# 	extended.segment@ranges = IRanges(start = from, end = to)
# 	
# 	# In order to make calculations faster if extend.left == 0 & extend.right == 0
# 	extend.nonzero = (extend.left != 0) | (extend.right != 0)
#   
# 	plotTrackList = plotTrackList(extended.segment, segment, model, modelAnnotation, contr, values, sampleAnnotation, genomeAnnotation, geneAnnotation, cpgAnnotation, addDataTrack = addDataTrack, showRegion = extend.nonzero, centeredData, extend.nonzero, ylim)
# 	
# 	if(is.na(filename)){
# 		plotTracks(plotTrackList, from = from, to = to, ...)
# 	}
# 	else{
# 		graphics_device(filename, width = width, height = height, units = units, pointsize = pointsize)
# 		plotTracks(plotTrackList, from = from, to = to, ...)
# 		dev.off()
# 	}
# }
# 	
#  
# #' Plot the significant segments that 
# #' 
# #' Function that draws the plots of significant segments. 
# #' 
# #' See examples in \code{\link{seqlm}}
# #'
# #' @param contr.res a table corresponding one contrast from the results of \code{\link{seqlm.contrast}}
# #' @param seqlmresults original seqlm results file
# #' @param min.length minimal length of region to draw
# #' @param coefficient minimal absolute coefficient value
# #' @param adjusted.pvalue threshold for adjusted pvalues
# #' @param annotation_column the figure is the column in the sample annotation file that we use to group the data
# #' @param genomic_coordinates boolean determining if the sites are drawn equidistantly or based on genomic coordinates
# #' @param expand number of extra probes to show on both sides of the region
# #' @param ylim two element vector giving the lower and higher limit of the y axis
# #' @param dir directory where to put the images, of NA then plots are drawn into the plotting window
# #' @param filetype picture filetype 
# #' @param ... extra parameters to \code{\link{ggsave}}
# #' 
# #' @author  Kaspar Martens <kmartens@@ut.ee> Raivo Kolde <rkolde@@gmail.com>
# #' 
# #' 
# #' 
# #' @export
# seqlmplots = function(contr.res, model, contr, seqlmresults, min.length = 10, coefficient = 0.2, adjusted.pvalue = 0.05, annotation_column = 1, geneAnnotation = NULL, cpgAnnotation = NULL, addDataTrack = FALSE, centeredData = FALSE, extend.left = 0, extend.right = 0, dir = NA, graphics_device = png, filetype = "png", width = 480, height = 480, units = "px",  pointsize = 12, ...){
# 	# Select only relevant regions
# 	contr.res = subset(contr.res, length >= min.length & abs(coef) >= coefficient & adjusted.p <= adjusted.pvalue)
# 		
# 	foreach(i = 1:length(contr.res)) %dopar% {
# 		if(is.na(dir)){
# 			filename = NA
# 		}
# 		else{
# 			filename = file.path(dir, sprintf("%d.%s", i, filetype))
# 		}
# 	
# 		seqlmplot(segment = contr.res[i, ], model = model, modelAnnotation = seqlmresults$description_length_par$annotation, contr = contr, values = seqlmresults$data$values, sampleAnnotation = seqlmresults$description_length_par$annotation[, annotation_column], genomeAnnotation = seqlmresults$data$genome_information, geneAnnotation = geneAnnotation, cpgAnnotation = cpgAnnotation, addDataTrack = addDataTrack, centeredData = centeredData, extend.left = extend.left, extend.right = extend.right, ylim = extendrange(seqlmresults$data$values), filename = filename, graphics_device = graphics_device, width = width, height = height, units = units, pointsize = pointsize, ...)
# 	}
# }
# 	
# # seqlmplots(pval[[1]], seqlmresults, annotation_column = "Factor2", coefficient = 0, min.length = 0, dir = "~/Desktop/")
# ##
# 
# ## seqlm report
# raport_template = '
# <!DOCTYPE html>
# <html>
# <head>
# <style type="text/css">.knitr.inline {
# 	background-color: #f7f7f7;
# 	border: solid 0px #b0b0b0
# }
# .message {
# 	font-style: italic
# }
# .source,.output,.warning,.error,.message {
# 	padding: 0em 1em;
# 	border: solid 1px #f7f7f7
# }
# .source {
# 	background-color: #f7f7f7
# }
# .rimage.left {
# 	text-align: left
# }
# .rimage.right {
# 	text-align: right
# }
# .rimage.center {
# 	text-align: center
# }
# .source {
# 	color: #333
# }
# .background {
# 	color: #f7f7f7
# }
# </style>
# <title>%s</title>
# </head>
# <body>
# 
# <code class="knitr inline">
# <h1> %s </h1>
# 
# %s
# </code>
# </body>
# </html>
# '
# chunk_template = '
# <h2> Segment %d </h2> 
# 
# <table>
# 	<tr>
# 		<td><b>Location</b></td>
# 		<td>%s</td>
# 	</tr>
# 	%s
# </table>
# 
# <div class="rimage default"><img src="%s" class="plot"/></div>
# '
# 
# 
# summary_template = '
# <h2> Chosen thresholds for this report </h2> 
# 
# <table>
# 	%s
# </table>
# 
# 
# <h2> Summary about all such regions </h2> 
# 
# <table>
# 	%s
# </table>
# 
# <h3> Summary about all such regions with positive coefficients </h3> 
# 
# <table>
# 	%s
# </table>
# 
# <h3> Summary about all such regions with negative coefficients </h3> 
# 
# <table>
# 	%s
# </table>
# '
# 
# annotation_template = '
# <tr>
# 	<td><b>%s</b></td>
# 	<td>%s</td>
# </tr>
# '
# 
# location_template = '%s:%d-%d'
# 
# annotation_table = function(x){
# 	n = which(colnames(x) == "adjusted.p")
# 	
# 	res = paste(sprintf(annotation_template, "Coefficient", round(x[1, "coef"], 3)), sprintf(annotation_template, "t-statistic", round(x[1, "tstat"], 3)), sprintf(annotation_template, "Adjusted p-value &thinsp;", round(x[1, "adjusted.p"], 3)), sprintf(annotation_template, "No. probes", x[1, "length"]), sprintf(annotation_template, "Length in bp", x[1, "width"]))
# 	
# 	if(!(n == ncol(x))){
# 		for(i in (n + 1):ncol(x)){
# 			res = paste(res, sprintf(annotation_template, colnames(x)[i], x[1, i]), sep = "\n")
# 		}
# 	}
# 	
# 	return(res)
# }
# 
# conditions_table = function(length, coefficient, adjusted.pvalue){
# 	s = paste(sprintf(annotation_template, "Condition: number of probes", ifelse(is.na(length), "all", paste("at least", length))), 
# 		sprintf(annotation_template, "Condition: absolute value of coefficient &thinsp;", ifelse(is.na(coefficient), "all", paste("at least", coefficient))), 
# 		sprintf(annotation_template, "Condition: adjusted.pvalue", ifelse(is.na(adjusted.pvalue), "all", paste("at most", adjusted.pvalue))))
# 	return(s)
# }
# 
# summary_table = function(df){
# 	if(nrow(df)==0) return(sprintf(annotation_template, "Number of such regions", 0))
# 	s = paste(sprintf(annotation_template, "Number of such regions", nrow(df)),
# 		sprintf(annotation_template, "&nbsp;", "&nbsp;"),
# 		sprintf(annotation_template, "Mean length (probes)", round(mean(df$length))), 
# 		sprintf(annotation_template, "Min length (probes)", min(df$length)), 
# 		sprintf(annotation_template, "Max length (probes)", max(df$length)),
# 		sprintf(annotation_template, "&nbsp;", "&nbsp;"),
# 		sprintf(annotation_template, "Mean length (bp)", round(mean(df$width))),
# 		sprintf(annotation_template, "Min length (bp)", min(df$width)),
# 		sprintf(annotation_template, "Max length (bp)", max(df$width)),
# 		sprintf(annotation_template, "&nbsp;", "&nbsp;"),
# 		sprintf(annotation_template, "Coefficient with max abs value &thinsp;", round(df$coef[which.max(abs(df$coef))], 3)))
# 	return(s)
# }
# 
#  
# #' Generate the HTML report for the seqlm results
# #' 
# #' Generate the HTML report for the seqlm results
# #'
# #' @param contr.res a table corresponding one contrast from the results of \code{\link{seqlm.contrast}}
# #' @param seqlmresults original seqlm results file
# #' @param length minimal length of region to draw
# #' @param coef minimal absolute coefficient value
# #' @param adjusted.pvalue threshold for adjusted pvalues
# #' @param orderby name of the contrast results table column to order results by
# #' @param annotation_column the figure is the column in the sample annotation file that we use to group the data
# #' @param genomic_coordinates boolean determining if the sites are drawn equidistantly or based on genomic coordinates
# #' @param expand number of extra probes to show on both sides of the region
# #' @param ylim two element vector giving the lower and higher limit of the y axis
# #' @param dir directory where to put the page, if the directory does not exist it will be created
# #' @param width picture width in inches
# #' @param height picture height in inches
# #' @param dpi dots per inch, to calibrate the picture size in pixels
# #' @param main title for the raport
# #' 
# #' @author  Kaspar Martens <kmartens@@ut.ee> Raivo Kolde <rkolde@@gmail.com>
# #' 
# #' @export
# seqlmreport = function(contr.res, seqlmresults, min.length = NA, coefficient = NA, adjusted.pvalue = 0.05, annotation_column = 1, geneAnnotation = NULL, cpgAnnotation = NULL, addDataTrack = FALSE, extend.left = 0, extend.right = 0, centeredData = FALSE, dir = NA, graphics_device = png, filetype = "png", width = 480, height = 480, units = "px", pointsize = 12, main = "Report: seqlm results"){
# 	# At the moment the option genomic_coordinates TRUE/FALSE has been removed
# 	
# 	length0 = ifelse(is.na(min.length), 0, min.length)
# 	coefficient0 = ifelse(is.na(coefficient), 0, coefficient)
# 	adjusted.pvalue0 = ifelse(is.na(adjusted.pvalue), 1, adjusted.pvalue)
# 
# 	# Select only relevant regions
# 	contr.res = subset(contr.res, length >= length0 & abs(coef) >= coefficient0 & adjusted.p <= adjusted.pvalue0)
# 	contr.df = as.data.frame(contr.res)
# 	
# 	if(length(contr.res) == 0) stop("Error. No such regions exist") else cat(sprintf("Number of such regions found: %s\n", length(contr.res)))
# 	
# 	# Create main directory 
# 	if(!file.exists(dir)){
# 		dir.create(dir)
# 	}
# 	
# 	# Create up and down regulated gene lists
# 	# gene_list = FALSE
# 	# if(!is.null(contr.df$UCSC_RefGene_Name)){
# 		# up = unlist(strsplit(as.character(contr.df$UCSC_RefGene_Name[contr.df$coef > 0]), ";"))
# 		# down = unlist(strsplit(as.character(contr.df$UCSC_RefGene_Name[contr.df$coef < 0]), ";"))
# 		# write.table(up, file = file.path(dir, "genes_up.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
# 		# write.table(down, file = file.path(dir, "genes_down.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
# 		# gene_list = TRUE
# 	# }
# 	
# 	# Create image directory
# 	img_dir = file.path(dir, "img")
# 	if(!file.exists(img_dir)){
# 		dir.create(img_dir)
# 	}
# 
# 	# Create images
# 	model = attr(contr.res, "model")
# 	contr = attr(contr.res, "contrast")
# 	seqlmplots(contr.res, model, contr, seqlmresults, min.length = 0, coefficient = 0, adjusted.pvalue = 1, annotation_column = annotation_column, geneAnnotation = geneAnnotation, cpgAnnotation = cpgAnnotation, addDataTrack = addDataTrack, extend.left = extend.left, extend.right = extend.right, centeredData = centeredData, dir = img_dir, graphics_device = graphics_device, filetype = filetype, width = width, height = height, units = units, pointsize = pointsize)
# 	
# 	# Create HTML file
# 	
# 	chunks = sprintf(summary_template, conditions_table(min.length, coefficient, adjusted.pvalue), summary_table(contr.df), summary_table(subset(contr.df, coef>0)), summary_table(subset(contr.df, coef<0)))
# 	
# 	# if(gene_list){
# 		# chunks = paste(chunks, '<h2> Gene lists </h2>',
# 						# '<p>', '<a href="genes_up.txt">Gene list (up regulated)</a> <br>', '</p>', '\n',
# 						# '<p>', '<a href="genes_down.txt">Gene list (down regulated)</a> <br>', '</p>', '\n')
# 	# }
# 	
# 	for(i in 1:length(contr.res)){
# 		location = sprintf(location_template, contr.df$seqnames[i], contr.df$start[i], contr.df$end[i])
# 		
# 		chunk = sprintf(chunk_template, i, location, annotation_table(contr.df[i, ]), sprintf("img/%d.png", i))
# 		
# 		chunks = paste(chunks, chunk, sep = "\n\n")
# 	}
# 	
# 	page = sprintf(raport_template, main, main, chunks)
# 	
# 	cat(page, file = file.path(dir, "index.html"))
# }
# 
# 
# ##
