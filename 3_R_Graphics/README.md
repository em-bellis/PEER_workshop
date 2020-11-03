# Session 3: R Graphics 
This session will provide an introduction to visualizing data with `ggplot2`.  We will follow pretty closely the [Data Analysis and Visualization in R for Ecologists Lesson](https://datacarpentry.org/R-ecology-lesson/index.html) from [Data Carpentry](https://datacarpentry.org/lessons/) which provides a much better and in-depth view of these topics. There are many other parts of this lesson you can check out if you want to. Today we will just focus on the 5th episode: [Data Visualization with ggplot2](https://datacarpentry.org/R-ecology-lesson/04-visualization-ggplot2.html). 

## Prior to the session: 
1. **Download and install R and R Studio**.  Follow the instructions [here](https://datacarpentry.org/R-ecology-lesson/#Install_R_and_RStudio). 

*Even if you already have R, it is recommended to ensure you are working in R version 4.0 or later. You can check by running the command `version` in the R console.*

---

## 3a: Load required packages and data
We will be using `ggplot2` which is part of the 'tidyverse'. We only need to install the package once, and then we can load the library each time we open R Studio.

All the following commands take place within the Console of R Studio.
```
> install.packages('tidyverse')
> library(tidyverse)
```
 
Since we are not following the Data Carpentry lesson from the beginning, we also need to download the example file separately. I have provided the example file in this repository, `surveys_complete.csv`, so carrying out the previous steps in the Data Carpentry lesson is not necessary.
```
> surveys_complete <- read_csv("surveys_complete.csv")
```

## 3b: Intro to `ggplot2`
We can now follow pretty closely the [Data Visualization with ggplot2 tutorial](https://datacarpentry.org/R-ecology-lesson/04-visualization-ggplot2.html) to practice making scatter plots, box plots, and time series plots with `ggplot2`.
