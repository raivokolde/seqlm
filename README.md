seqlm
=======

An R package for identification of differentially methylated regions (DMRs) from high density chip, for example Illumina 450K, data. 

Method
------
The seqlm method works in three stages. 

**Stage 1:** The genome is divided into smaller pieces based on a genomic distance cutoff. 

**Stage 2:** In each piece probes are segmented into regions that have approximately constant difference between the groups of interest. Example of the segmentation and its process is shown in [schema].

* In sliding windows with variable sizes we fit a linear models to the data.
* For each model we record the description length - the amount of bits needed to describe the data using the model
* Using dynamic programming we find the segmentation that minimizes total description length

**Stage 3:** We assess the relevance of each segment, by using a mixed model where the classes are a fixed effect and a sample is a random effect. This model takes into account the repeated nature of the consecutive methylation measurements. The segments are ordered by their significance.

![Example of seqlm segmentation][schema]

[schema]: https://raw.github.com/raivokolde/seqlm/gh-pages/pics/schema.png "Example of seqlm segmentation"

Installation
------------
The most convenient way to install the package is by using the `devtools` package.

```r
library(devtools)
install_github("seqlm", "raivokolde")
```







