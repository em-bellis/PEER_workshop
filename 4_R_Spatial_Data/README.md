# Session 4: Spatial Data 
This session will provide an introduction to visualizing data with `ggplot2`.  We will follow pretty closely the [Data Analysis and Visualization in R for Ecologists Lesson](https://datacarpentry.org/R-ecology-lesson/index.html) from [Data Carpentry](https://datacarpentry.org/lessons/) which provides a much better and in-depth view of these topics. There are many other parts of this lesson you can check out if you want to. Today we will just focus on the 5th episode: [Data Visualization with ggplot2](https://datacarpentry.org/R-ecology-lesson/04-visualization-ggplot2.html). 

## Prior to the session: 
1. **Download and install R and R Studio**.  Follow the instructions [here](https://datacarpentry.org/R-ecology-lesson/#Install_R_and_RStudio). 

*Even if you already have R, it is recommended to ensure you are working in R version 4.0 or later. You can check by running the command `version` in the R console.*

---

5. Get western Kenya county boundaries for map:
# partially following tutorial here https://rpubs.com/spoonerf/countrymapggplot2
library(maptools)
library(raster)
library(plyr)
library(ggplot2)
library(rgdal)
# get polygons for western Kenyan counties
Kenya1&lt;-getData(&quot;GADM&quot;, country=&quot;KE&quot;, level=1)
Kenya1_UTM&lt;-spTransform(Kenya1, CRS(&quot;+init=EPSG:32737&quot;))
counties &lt;-Kenya1_UTM[Kenya1_UTM@data$NAME_1 == &quot;Busia&quot;| Kenya1_UTM@data$NAME_1 == &quot;Kisumu&quot; |
Kenya1_UTM@data$NAME_1==&quot;Homa Bay&quot;| Kenya1_UTM@data$NAME_1 == &quot;Kericho&quot;|Kenya1_UTM@data$NAME_1 ==
&quot;Nandi&quot;|Kenya1_UTM@data$NAME_1 == &quot;Vihiga&quot;|Kenya1_UTM@data$NAME_1 == &quot;Siaya&quot;|Kenya1_UTM@data$NAME_1
== &quot;Kakamega&quot;|Kenya1_UTM@data$NAME_1 == &quot;Kisii&quot;,]
counties.ll &lt;- spTransform(counties, CRS(&quot;+proj=longlat&quot;))
