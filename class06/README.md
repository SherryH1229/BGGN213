Class\_6
================
Xuerui HUang
4/19/2019

``` r
library(dplyr)
```

# Section 1: A

Try to load files wit different
deliminators

``` r
table1 <- read.table("test1.txt",sep = ",",header = TRUE) # seperated by ,
table1
```

    ##   Col1 Col2 Col3
    ## 1    1    2    3
    ## 2    4    5    6
    ## 3    7    8    9
    ## 4    a    b    c

``` r
table2 <- read.table("test2.txt",sep = "$",header = TRUE) #seperated by $
table2
```

    ##   Col1 Col2 Col3
    ## 1    1    2    3
    ## 2    4    5    6
    ## 3    7    8    9
    ## 4    a    b    c

``` r
table3 <- read.table("test3.txt",header = TRUE)
table3
```

    ##   X1 X6 a
    ## 1  2  7 b
    ## 2  3  8 c
    ## 3  4  9 d
    ## 4  5 10 e

write the first function to add 1 to the input value

``` r
add <- function(x,y = 1){
  x+y
}

add(1)
```

    ## [1] 2

write the function for rescalling the input vector

``` r
rescale <- function(x) {
   rng <-range(x)
   (x - rng[1]) / (rng[2] - rng[1])
}

rescale(1:10)
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

# Section1:B

Q1: read.pdb() function return sa large pdb object

Q2: trim.pdb() function trims a PDB object to a subset of atoms, which
produces a new smaller PDB object, containing a subset of atoms, from a
given larger PDB object.

Q3: To turn off the marginal black and grey rectangles in the plots, we
could set the top = FALSE and bot = FALSE\> The rectangles represents
beta strands and alpha helics

Q4: To put plots together

Q5: S1 and S3 are closer

## Q6

**Input**: Two inputs
(FileName,chain\_name)

``` 
   FileName: a single element character vector containing the name of the PDB file 
               to be read, or the four letter PDB identifier for online file access
   chain_name: a single letter string to indicate which chain you want to focus on
```

**Function**: Visualizing and analyzing the specific protein drug
interactions by inputing the specific PBD data indicators

**Ouput**: A plot object for the specified protein

``` r
#Input: Two inputs (FileName,chain_name)
#       FileName: a single element character vector containing the name of the 
#                 PDB file to be read, or the four letter PDB identifier for 
#                 online file access
#       chain_name: a single letter string to indicate which chain you want 
#                   to focus on
#Function: Visualizing and analyzing the specific protein drug interactions by i
#          nputing the specific PBD data indicators
#Ouput: A plot object for the specified protein
plot_drugProteinInteract <- function (fileName,chain_name){
  require(bio3d)
  #load file
  pdb_pf <- read.pdb(fileName)
  #select specific chain and info from chain
  pf_chain <- trim.pdb(pdb_pf,chain = chain_name,elety="CA")
  pf.b <- pf_chain$atom$b
  
  #output plot
  plotb3(pf.b, sse=pf_chain, typ="l", ylab="Bfactor")
}

#test
plot_drugProteinInteract("4AKE","A")
```

    ## Loading required package: bio3d

    ##   Note: Accessing on-line PDB file

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->
