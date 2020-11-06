# Session 4: Spatial Data 
In this session, we will create a map of Kenya, with outlines of the counties in western Kenya. On this map we will overlay a set of points, in this case locations where *Striga hermonthica* was collected in 2018.  Much of the code is adapted from the tutorial [here](https://rpubs.com/spoonerf/countrymapggplot2). 

For continued learning with more in-depth resources, Data Carpentry also has a tutorial for geospatial data in R [here](https://datacarpentry.org/r-raster-vector-geospatial/), though we won't really be drawing on it much today. My favorite tutorials for using R for geospatial data analysis come from [rspatial.org](https://rspatial.org/raster/index.html).

## Prior to the session: 
1. **Background reading on geospatial data**.  The material [here](https://datacarpentry.org/organization-geospatial/) has some great information to introduce vector and raster data, coordinate reference systems, and software for working with geospatial data. 

---

## 4a. Get data for western Kenya county boundaries:
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
