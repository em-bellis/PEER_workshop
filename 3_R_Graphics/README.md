# Session 3: R Graphics 
This session is designed in two parts. First, we will provide an brief overview of R as a programming language based on the Chapter 2 Lab from [Introduction to Statistical Learning](http://faculty.marshall.usc.edu/gareth-james/ISL/code.html) by James, Witten, Hastie and Tibshirani, available as a free pdf.  There is also a [Programming with R](http://swcarpentry.github.io/r-novice-inflammation/) lesson from [Software Carpentry](https://software-carpentry.org).

Then, we will introduce R as a software environment, learning to visualize data with `ggplot2`.  We will follow pretty closely the [Data Analysis and Visualization in R for Ecologists Lesson](https://datacarpentry.org/R-ecology-lesson/index.html) from [Data Carpentry](https://datacarpentry.org/lessons/) for this second part. There are many other parts of this lesson you can check out if you want to. Today we will just focus on the 5th episode: [Data Visualization with ggplot2](https://datacarpentry.org/R-ecology-lesson/04-visualization-ggplot2.html). 

If time permits, we will test our skills by importing the vcf file we made yesterday, and visualizing coverage and other quality information using the [vcfR](https://cran.r-project.org/web/packages/vcfR/vignettes/intro_to_vcfR.html) package.

## Prior to the session: 
1. **Download and install R and R Studio**.  Follow the instructions [here](https://datacarpentry.org/R-ecology-lesson/#Install_R_and_RStudio). 

*Even if you already have R, it is recommended to ensure you are working in R version 4.0 or later. You can check by running the command `version` in the R console.*

---
All the following commands take place within the console of R/R Studio.

## Part I: R as a programming language
### 3a: Basic commands
R uses functions to perform operations. For example, the function `c()` can be used to create a vector of numbers. The assignment operator `<-` tells R to store the vector of numbers as an object named `x`.
```
> x <- c(1,3,2,5)
> x
```
To pull up the documentation for a function, in R we use `?`.
```
> ?c
```
We can tell R to add numbers, or even vectors together
```
> 4 + 5
> y <- c(1,4,3)
> x + y
```
We can use `ls()` to list all objects that we have saved so far, and `rm()` to remove
```
> ls()
> rm(x,y)
```
We can create a matrix of numbers to learn how R creates objects like matrices
```
> x <- matrix(data=c(1,2,3,4), nrow=2, ncol=2)
> x # the hashtag precedes any comments that are not interpreted by R
> x <- matrix(c(1,2,3,4), 2, 2, byrow=2) # we don't always need to specificy the names of the arguments passed in; byrow fills the matrices in order of rows
> sqrt(x)
> x^2
```
It can be useful to generate a vector of random normal variables. For example, here we can create two correlated sets of numbers, and calculate some basic statistics
```
> x <- rnorm(50) # by default, mean = 0, sd = 1
> y <- rnorm(50, mean = 50, sd=0.1)
> cor(x,y)
> mean(x)
> var(x)
```
If we want to reproduce the exact same set of random numbers, we can use `set.seed()`
```
> set.seed(50)
> rnorm(50)
```

### 3b. Graphics
We can graph data very easily in base R. We will get more into visualization later on with `ggplot`...
```
> plot(x, y)
> plot(x, y, xlab="this is the x", ylab="this is the y", main="plot of x vs y", col="red")
```
You can save figures as a pdf
```
> pdf("Figure1.pdf")
> plot(x,y, col="green")
> dev.off()
```

### 3c. Indexing Data
```
> A <- matrix(1:16, 4, 4)
> A
> A[2,3] # selects element in 2nd row and 3rd column
> A[1:3, 2:4] # select multiple rows at a time
> A[1:2, ] # select first two rows and all columns
> A[, 1:2] # select all rows and first two columns
> A[-c(1,3),] # keep all rows except 1 and 3 and all columns
> dim(A) # the dimensions of A
```
### 3d. Loading Data 
```
> getwd()
> setwd("/path/to/folder/") # optional
> surveys_complete <- read.csv("surveys_complete.csv")
```

Once you have loaded a data file, there are a few ways to inspect it.
```
> names(surveys_complete) # what variables are included?
> str(surveys_complete # this is an important one! 99% of R issues are caused by errors in numbers being interpreted as factors...
> summary(surveys_complete) # missing data? can be removed with na.omit
> dim(surveys_complete) # what are the dimensions?
> class(surveys_complete) # what kind of object is this?
```

With data frames, you can select a particular variable with `$`
```
> surveys_complete$species
> surveys_complete$species[2] # select just the second element of this vector
> unique(surveys_complete$species) # unique elements
```
## Part II: R as a software environment
### 3e: Load required packages and data
We will be using `ggplot2` which is part of the 'tidyverse'. We only need to install the package once, and then we can load the library each time we open R Studio. The `tidyverse` packages includes several packages which are part of the tidyverse ecosystem.
```
> install.packages('tidyverse')
> library(tidyverse)
```
 
Since we are not following the Data Carpentry lesson from the beginning, we also need to download the example file separately. I have provided the example file in this repository, `surveys_complete.csv`, so carrying out the previous steps in the Data Carpentry lesson is not necessary.
```
> surveys <- read_csv("surveys_complete.csv") # part of tidyverse package, a bit faster than read.csv for large files
```

We will be spending most of our time today on plotting, but there are some excellent packages in the tidyverse ecosystem (e.g. `tidyr`, `dplyr`) for data wrangling, e.g.:
```
> surveys %>%
   filter(weight < 5) %>%
   select(species_id, sex, weight)
```
[Section 4](https://datacarpentry.org/R-ecology-lesson/03-dplyr.html) of the Data Analysis and Visualization in R for Ecologists lesson has a great intro to data wrangling with `dplyr` and `tidyr`. See also the free [R for Data Science)[https://r4ds.had.co.nz] book.

## 3f: Intro to `ggplot2`
We can now follow pretty closely the [Data Visualization with ggplot2 tutorial](https://datacarpentry.org/R-ecology-lesson/04-visualization-ggplot2.html) to practice making scatter plots, box plots, and time series plots with `ggplot2`.

For more tips on making graphs with `ggplot`, [The R Graphics Cookbook](https://r-graphics.org) is available online for free and is really great!
