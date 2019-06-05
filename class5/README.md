-   [load required libararies](#load-required-libararies)
-   [2A. Line plot](#a.-line-plot)
-   [2B. Barplot](#b.-barplot)
-   [2C. Histograms](#c.-histograms)

load required libararies
========================

``` r
library(dplyr)
```

2A. Line plot
=============

The file weight\_chart.txt from the example data you downloaded above
contains data for a growth chart for a typical baby over the first 9
months of its life.

load data and plot line plot based on the info

``` r
#load data
weight <- read.table("./bimm143_05_rstats/weight_chart.txt",header=TRUE)

#doing plot 
plot(weight,pch=15,cex=1.5,lwd=2,ylim=c(2,10),xlab="Age (months)",
     ylab="Weight (kg)",main="Baby Weight with Age")
lines(weight,lwd=2)
```

![](README_files/figure-markdown_github/unnamed-chunk-2-1.png)

2B. Barplot
===========

The file feature\_counts.txt contains a summary of the number of
features of different types in the mouse GRCm38 genome.

load data and plot line barplot based on the info

``` r
# load data
feature <- read.csv("./bimm143_05_rstats/feature_counts.txt",sep = "\t") %>% 
  as.data.frame(.)

# check function documentation for further modification of plot
?barplot

#plot barplot
par(mar=c(5,10,4,2)+0.1,mgp=c(5,1,0)) #change margin
barplot(feature$Count,names.arg = feature$Feature,cex.names=0.8,
        xlab = "Count",las = 2,
        horiz=TRUE,las = 1,main="Num_Features_GRCm38")
```

![](README_files/figure-markdown_github/unnamed-chunk-3-1.png)

2C. Histograms
==============

Plot out the distribution of 10000 points sampled from a standard normal
distribution along with another 10000 points sampled from the same
distribution but with an offset of 4.

``` r
par(mar=c(5.1,4.1,4.1,2.1))
hist(c(rnorm(10000),rnorm(10000)+4),main = "Histogram",breaks=100,xlab = "x")
```

![](README_files/figure-markdown_github/unnamed-chunk-4-1.png)
