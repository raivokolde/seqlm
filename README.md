seqlm
=======

An R package for identification of differentially methylated regions (DMRs) from high density chip, for example Illumina 450K, data. 

Installation
------------
The most convenient way to install the package is by using the `devtools` package.

```s
library(devtools)
install_github("raivokolde/seqlm")
```
To start using the package just load it as any other package.

```s
library(seqlm)
```

Usage
-----
For running the algorithm one has to have three objects:

* a matrix with methylation values;
* a vector specifying the classes of columns (only two-class case is supported currently) or a continuous variable;
* location information about the methylation probes in GRanges format (a file for Illumina 450K platform can be downloaded from [here](http://biit.cs.ut.ee/~kolde/seqlm/genome_information.RData)). 

All the work is done by one command `seqlm` that takes as input all the objects described above and also parameters:

* `max_block_length` - that determines the maximal number of CpG sites in a region and is used to speed up the computations (default 50)
* `max_dist` - that determines the maximal allowed genomic distance between two consecutive probes in a region (default 1000bp)

For more information use `?seqlm` in R.

Optionally one can visualise the results using `seqlmreport` command that visualises the identified regions and creates a html page to show them. The main input for this function is the result of `seqlm` and the three objects described above. 

Example
-------

An example dataset `tissue_small` is included in the package. It contains data comparing adipose tissue and brainstem, from chromosomes 17 and 18. Finding regions in this data can be accomplished using commands:

```s
data(tissue_small)
segments = seqlm(values = tissue_small$values, genome_information = tissue_small$genome_information, annotation =  tissue_small$annotation)
```

The result of the analysis is a GRanges object containing the locations of the regions and associated statistics. 

The analysis can be time consuming, if the whole genome is analysed at once. If the computer has multicore capabilities it is easy to parallelize the calculations. We use the [foreach](http://cran.r-project.org/web/packages/foreach/index.html) framework by Revolution Computing for parallelization. To enable the parallelization one has to register the parallel backend before and this will be used by seqlm. Ideally, the next commands should take roughly half the time compared to the previous.

```s
library(doParallel)
registerDoParallel(cores = 2)
segments = seqlm(values = tissue_small$values, genome_information = tissue_small$genome_information, annotation =  tissue_small$annotation)
```

To visualise the results it is possible to plot the most imortant sites and generate a HTML report

```s
temp = tempdir()
seqlmreport(segments[1:10], tissue_small$values, tissue_small$genome_information, tissue_small$annotation, dir = temp)
```

[Here](http://biit.cs.ut.ee/~kolde/seqlm/seqlmreport/index.html) is an example of the resulting file.

Method
------
The method is described in a [poster](https://raw.github.com/raivokolde/seqlm/gh-pages/pics/Poster_uus.pdf), that was presented at the Epigenetics of Common Diseases conference. Briefly, the seqlm method works in three stages. 

**Stage 1:** The genome is divided into smaller pieces based on a genomic distance cutoff. 

**Stage 2:** In each piece probes are segmented into regions that have approximately constant difference between the groups of interest. Example of the segmentation and its process is shown in [schema].

* In sliding windows with variable sizes we fit a linear models to the data.
* For each model we record the description length - the amount of bits needed to describe the data using the model
* Using dynamic programming we find the segmentation that minimizes total description length

**Stage 3:** We assess the relevance of each segment, by using a mixed model where the classes are a fixed effect and a sample is a random effect. This model takes into account the repeated nature of the consecutive methylation measurements. The segments are ordered by their significance.

![Example of seqlm segmentation][schema]

[schema]: https://raw.github.com/raivokolde/seqlm/gh-pages/pics/schema.png "Example of seqlm segmentation"






