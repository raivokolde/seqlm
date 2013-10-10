segmentation = function(inputlist){
	dl = inputlist[[1]]
	max_block_length = inputlist[[2]]
	m = nrow(dl)
	
	# Dynamic programming
	S = rep(0, m + 1)
	I = rep(0, m)

	for(j in 1:m){
		start = ifelse(j>max_block_length, j+1-max_block_length, 1)
		costs = S[start:j] + dl[start:j, j]
		S[j+1] = min(costs)
		i = which.min(costs) + start - 1
		I[j] = i - 1
	}
	
	# Identify Regions
	k = m
	pairs = numeric(0)
	while(k > 0){
		pairs = rbind(pairs, c(I[k] + 1, k))
		k = I[k]
	}
	pairs = as.data.frame(pairs[nrow(pairs):1, , drop=FALSE])
	names(pairs) = c("start", "end")
	return(pairs)
}


group_by_dist = function(pos, max_dist_cpg){
	dist = data.frame("d" = pos[-1] - pos[-length(pos)])
	dist$temp = (dist$d > max_dist_cpg | dist$d < 0)+0
	indicator = as.factor(c(0, cumsum(dist$temp)))
	return (indicator)
}

match_positions = function(values, genome_information){
	genome_information$pos = as.numeric(as.character(genome_information$pos))
	genome_information = genome_information[!is.na(genome_information$pos),]
	
	int = intersect(rownames(genome_information), rownames(values))
	genome_information = genome_information[int, ]
	values = values[int, ]
	
	genome_information = genome_information[order(genome_information$chr, genome_information$pos), ]
	values = values[rownames(genome_information), ]
	
	return(list(values = values, genome_information = genome_information))
}

#' Fit the model
#' 
#' Fit the model
#' 
#' The \code{genome_information} variable is expected to be a 
#' \code{data.frame} that gives genomic position and optional extra 
#' information for every row in the \code{values}. To identify the 
#' matrix rows then the rownames in the \code{values} and 
#' \code{genome_information} have to match. The \code{genome_information}
#'  has to have at least two columns: "chr" for chromosome and "pos" showing 
#' the coordinate. Additional columns in the \code{genome_information} table
#'  are used to describe the regions afterwards.    
#'
#' @param values a matrix where columns are samples and rows correspond to the sites
#' @param genome_information \code{data.frame} giving the genomic coordinates and optionally additional 
#' description for the sites.
#' @param max_dist maximal genomic distance between the sites to be considered the same region
#' @param max_block_length maximal length for a block 
#' @param max_block_length_second_stage maximal block length for the second stage of search if two-stage search is used 
#' for speeding up the analysis. Second stage is initiated only if \code{max_block_length_second_stage} > 
#' \code{max_block_length}
#' @param description_length_fun function for calculating the description length
#' @param description_length_par additional parameters for \code{description_lengths_fun}
#' 
#' @return  A list containing the input data, parameters and the segmentation.
#' 
#' @author  Kaspar Martens <kmartens@@ut.ee> Raivo Kolde <rkolde@@gmail.com>
#' 
#' @export
seqfit = function(values, genome_information, max_dist, max_block_length, description_length_fun, description_length_par){
	# Center the rows of "values"
	values = t(scale(t(values), center=TRUE, scale=FALSE))
	
	mp = match_positions(values, genome_information)
	values = mp$values
	genome_information = mp$genome_information
	gr = GRanges(seqnames = genome_information$chr, ranges = IRanges(start = genome_information$pos, width = 1, names = rownames(genome_information)))
	if(ncol(genome_information) > 2){
		which.cols = which(colnames(genome_information) %in% c("chr", "pos"))
		elementMetadata(gr) = as(genome_information[, -which.cols], "DataFrame")
	}
	
	# Divide the genome into initial segments based on genomic coordinate
	chr = genome_information$chr
	pos = genome_information$pos
	indicator = group_by_dist(pos, max_dist)
	pieces = split(seq_along(indicator), indicator)
	lengths = sapply(pieces, length)
	
	ord = order(-lengths)
	pieces = pieces[ord]
	lengths = lengths[ord]

	maxlength = max(lengths)
	which_pieces1 = which(lengths == 1)
	which_pieces2 = which(lengths >= 2)
	
	cat("Finding the best segmentation\n")
	flush.console()
	pb <- txtProgressBar(min = 0, max = length(which_pieces2), style = 3)

	# Segment based on the model
	segmentlist = foreach(i = seq_along(which_pieces2), .export=c("description_length_fun", "prior_var", "row_fit", "row_model_spec", "row_subsets", "row_bics", "vars", "segmentation", "calculate_rss", "number_of_coefs")) %dopar% {
		# Which rows from the matrix "values" belong to this piece
		u = pieces[[which_pieces2[i]]]
		
		# Input parameters
		par = description_length_par
		par$values = values[u, , drop=FALSE]
		par$maxlength = maxlength
		par$max_block_length = max_block_length
		
		# Call the partial description length function on whole region
		a = do.call("description_length_fun", par)
		result = segmentation(a)
		
		if (i%%100==0){
			setTxtProgressBar(pb, i)
		}
		# In the first column there are "startIndexes", in the second column "endIndexes"
		matrix(c(u[result$start], u[result$end]), ncol=2)
	}
	setTxtProgressBar(pb, length(which_pieces2))
	close(pb)

	index_pieces1 = unlist(pieces[which_pieces1])
	
	res = do.call("rbind", segmentlist)
	res = rbind(res, matrix(c(index_pieces1, index_pieces1), ncol=2))
	colnames(res) = c("startIndex", "endIndex")
	res = as.data.frame(res)
	res = res[order(res$startIndex), ]
	segments = cbind("chr" = chr[res$startIndex], "startPos" = pos[res$startIndex], "endPos" = pos[res$endIndex], 
					 "length" = res$endIndex - res$startIndex + 1, res)

	
	# Compile output
	output = list(
		data = list(
			values = values,
			genome_information = gr
		),
		description_length_par = description_length_par,
		segments = segments
	)
	
	# Add segment information if necessary
	# if (ncol(genome_information) > 2){
		# cat("Adding annotations to segments\n")
		# output$segment_annotation = additional_segment_information(output)
	# }
	
	return (output)
}
