# Session 4: Spatial Data 
In this session, we will create a map of Kenya, with outlines of the counties in western Kenya. On this map we will overlay a set of points, in this case locations where leaves of *Striga hermonthica* were collected in 2018 for genetic sequencing.  Much of the code is adapted from the tutorial [here](https://rpubs.com/spoonerf/countrymapggplot2). 

For continued learning, Data Carpentry also has a tutorial for geospatial data in R [here](https://datacarpentry.org/r-raster-vector-geospatial/), though we won't really be drawing on it much today. My favorite tutorials for using R for geospatial data analysis come from [rspatial.org](https://rspatial.org/raster/index.html).

## Prior to the session: 
1. **Background reading on geospatial data**.  The material [here](https://datacarpentry.org/organization-geospatial/) has some great information to quickly introduce vector and raster data, coordinate reference systems, and software for working with geospatial data. 

2. **Install some of the required R packages**. This can be done directly within the R Studio environment through the 'Packages' tab of the bottom left panel or running the following:
```
> install.packages()
```
---

## 4a. Load required packages:
```
library(maptools)
library(raster)
library(plyr)
library(ggplot2)
library(rgdal)
```

## 4b. Get polygons for Kenya county boundaries:
It can take a couple minutes to download unfortunately...
```
Kenya1<-getData("GADM", country="KE", level=1)
```

Also get a map just of the outline of Kenya:
```
```

```
Kenya1_UTM<-spTransform(Kenya1, CRS("+init=EPSG:32737"))  
```

## 4c. Choose a subset of counties to include in the map:
```
counties <-Kenya1_UTM[Kenya1_UTM@data$NAME_1 == "Busia"| Kenya1_UTM@data$NAME_1 == "Kisumu" | Kenya1_UTM@data$NAME_1=="Homa Bay"|  Kenya1_UTM@data$NAME_1 == "Kericho"|Kenya1_UTM@data$NAME_1 == "Nandi"|Kenya1_UTM@data$NAME_1 == "Vihiga"|Kenya1_UTM@data$NAME_1 == "Siaya"|Kenya1_UTM@data$NAME_1 == "Kakamega"|Kenya1_UTM@data$NAME_1 == "Kisii",]
````

```
counties.ll <- spTransform(counties, CRS("+proj=longlat"))
```

## 4d. Plot a base map of Kenya, with outlines of the counties of interest

## 4e. Add points to show sampling locations


