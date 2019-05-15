Untitled
================
Xuerui HUang
5/15/2019

# Get data

``` r
library(seqinr)
library(gtools)
phred <- asc( s2c("DDDDCDEDCDDDDBBDDDCC@") ) - 33
phred
```

    ##  D  D  D  D  C  D  E  D  C  D  D  D  D  B  B  D  D  D  C  C  @ 
    ## 35 35 35 35 34 35 36 35 34 35 35 35 35 33 33 35 35 35 34 34 31

``` r
test <- read.csv("https://bioboot.github.io/bggn213_W19/class-material/rs8067378_ENSG00000172057.6.txt",sep = " ")

summary(test)
```

    ##      sample     geno          exp        
    ##  HG00096:  1   A/A:108   Min.   : 6.675  
    ##  HG00097:  1   A/G:233   1st Qu.:20.004  
    ##  HG00099:  1   G/G:121   Median :25.116  
    ##  HG00100:  1             Mean   :25.640  
    ##  HG00101:  1             3rd Qu.:30.779  
    ##  HG00102:  1             Max.   :51.518  
    ##  (Other):456

``` r
#ggplot(test[c(2,3)],a)
boxplot(test$exp,group = test$geno)
```

![](Untitled_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
inds <- (test$geno=="A/G")
summary(test[inds,]$exp)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   7.075  20.626  25.065  25.397  30.552  48.034

``` r
boxplot(exp~geno,data = test,notch = TRUE)
```

![](Untitled_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->
