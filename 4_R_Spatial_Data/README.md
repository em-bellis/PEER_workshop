# Session 4: Spatial Data 
In this session, we will create a map of Kenya, with outlines of the counties in western Kenya. On this map we will overlay a set of points, in this case locations where leaves of *Striga hermonthica* were collected in 2018 for genetic sequencing.  Much of the code is adapted from the tutorial [here](https://rpubs.com/spoonerf/countrymapggplot2). 

For continued learning, Data Carpentry also has a tutorial for geospatial data in R [here](https://datacarpentry.org/r-raster-vector-geospatial/), though we won't really be drawing on it much today. My favorite tutorials for using R for geospatial data analysis come from [rspatial.org](https://rspatial.org/raster/index.html).

## Prior to the session: 
1. **Background reading on geospatial data**.  The material [here](https://datacarpentry.org/organization-geospatial/) has some great information to quickly introduce vector and raster data, coordinate reference systems, and software for working with geospatial data. 

2. **Install some of the required R packages**. This can be done directly within the R Studio environment through the 'Packages' tab of the bottom left panel or running the following:
```
> install.packages(c("rgdal", "raster", "sp", "sf", "tmap","tmaptools"), dependencies = T)
```
---

## 4a. Load required packages
```
> library(raster)
> library(rgdal)
> library(sp)
> library(sf)
> library(ggplot2)
> library(tmap)
> library(tmaptools)
```

## 4b. Get polygons for Kenya county boundaries
It can take a couple minutes to download unfortunately...
```
> Kenya1<-getData("GADM", country="KE", level=1)
```

Also get a polygon just of the outline of Kenya
```
> Kenya<-getData("GADM", country="KE", level=0)
```

Reproject to UTM coordinate reference system. UTM is a good choice for a regional map covering a relatively small area.
```
> Kenya1_UTM<-spTransform(Kenya1, CRS("+init=EPSG:32737")) 
> Kenya0_UTM<-spTransform(Kenya, CRS("+init=EPSG:32737")) 
```

You can check that the reprojection worked by inspecting each Spatial Polygon Data Frame
```
> Kenya1_UTM
> Kenya1
```

## 4c. Choose a subset of counties to include in the map:
Look at the names of the counties and choose some we'd like to plot
```
> Kenya1_UTM@data$NAME_1
```

We can use the `[` and `]` bracket characters to subset the `Kenya1_UTM` object to just Busia county. We want to include all columns in the data frame associated with Busia where `data$NAME_1` matches "Busia" (not just individual elements of the object), so we should be sure to include the `,` before the right bracket.
```
> Kenya1_UTM[Kenya1_UTM@data$NAME_1 == "Busia",]
```

We can also use the logical OR `|` character to match multiple counties 
```
counties <-Kenya1_UTM[Kenya1_UTM@data$NAME_1 == "Busia"| Kenya1_UTM@data$NAME_1 == "Kisumu",]
```

## 4d. Plot a base map of Kenya, with outlines of the counties of interest
[tmaps](https://cran.r-project.org/web/packages/tmap/vignettes/tmap-getstarted.html)
```
counties_sf <- as(counties, Class = "sf")  #convert to sf class
kenya_sf <- as(Kenya0_UTM, Class="sf")

kenya_smooth <- simplify_shape(kenya_sf, 0.001) # reduce number of coordinates in polygon to 'smooth' shape and speed up plotting.

tm_shape(kenya_smooth) +
  tm_polygons() +
tm_shape(counties) + 
  tm_polygons(col="tomato4") + 

```

## 4e. Add points to show sampling locations


