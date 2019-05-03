lecture7\_proj: R: Functions and packages
================
Xuerui HUang
4/24/2019

R Markdown
----------

``` r
source("http://tinyurl.com/rescale-R")
```

``` r
#test the rescale function
rescale(1:10)
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

``` r
#test function
x <- c(1:10,"string")
is.numeric(x)
```

    ## [1] FALSE

Function for recognize NA in two vectors
----------------------------------------

``` r
#sample input
x <- c( 1, 2, NA, 3, NA)
y<-c(NA,3,NA,3, 4)

sum(is.na (x)&is.na(y))
```

    ## [1] 1

``` r
#function for checking whether two elements on the same position of two vectors have NA
both_na <- function(x,y){
  sum(is.na(x)&is.na(y))
}
```

function for grading homework
-----------------------------

``` r
student1 <- c(rep(100,time= 8),90)
student2 <- c(100,NA,rep(90,time = 4),97,80)

#forloop
overall_grade <- function(stud_vec){
  stud_vec <- stud_vec[!is.na(stud_vec)]
  min_val <- min(stud_vec)
  for (i in 1:length(stud_vec)){
    if (stud_vec[i] == min_val){
      stud_vec <- stud_vec[-i]
      break
    }
  }
  mean_score <- mean(stud_vec)
  return (mean_score)
}

#better way
overall_grade_2 <- function(stud_vec){
  if (sum(is.na(stud_vec) >0)){
    stud_vec <- stud_vec[!is.na(stud_vec)]
  }
  
  (sum(stud_vec)-min(stud_vec))/(length(stud_vec)-1)
}

overall_grade(student1)
```

    ## [1] 100

``` r
overall_grade(student2)
```

    ## [1] 92.83333

``` r
overall_grade_2(student1)
```

    ## [1] 100

``` r
overall_grade_2(student2)
```

    ## [1] 92.83333

``` r
#library(bio3d)

#lbio3d()

# demo(package="bio3d")
# demo("nma")
```
