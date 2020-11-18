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

## 4c. Choose a subset of counties to include in the map
Look at the names of the counties and choose some we'd like to plot
```
> Kenya1_UTM@data$NAME_1
```

We can use the `[` and `]` bracket characters to subset the `Kenya1_UTM` object to just Busia county. We want to include all variables in the data frame associated with Busia where `data$NAME_1` matches "Busia" (not just individual elements of the object), so we should be sure to include the `,` before the right bracket.
```
> Kenya1_UTM[Kenya1_UTM@data$NAME_1 == "Busia",]
```

We can use the logical OR character, `|`, to match multiple counties 
```
> counties <-Kenya1_UTM[Kenya1_UTM@data$NAME_1 == "Busia"| Kenya1_UTM@data$NAME_1 == "Kisumu",]
```

## 4d. Plot a base map of Kenya, with outlines of the counties of interest
There are many strategies for making maps in R! I was recently introduced to the [tmap](https://cran.r-project.org/web/packages/tmap/vignettes/tmap-getstarted.html) package and really like it so far compared to the ways I used to make maps in R :).

While `tmap` can handle SpatialPolygonsDataFrames like the `Kenya1_UTM` object we created, the `sf` class of objects is better supported. We can use the `as` function to convert to `sf` class 
```
> counties_sf <- as(counties, Class = "sf")
> kenya_sf <- as(Kenya0_UTM, Class="sf")
```

An additional benefit of conversion to the `sf` class is that we can use the `simplify_shape` function from `tmaptools` to reduce the number of coordinates in our polygon objects. This will 'smooth' the outline of the shapes and make for much faster plotting!
```
> kenya_smooth <- simplify_shape(kenya_sf, 0.01)
```

Now, make a basic `tmap`, just showing the 'smoothed' outline of Kenya. The syntax of `tmap` is somewhat similar to `ggplot2`. For each 'shape' object (in this case, `kenya_smooth` is the shape object) we can plot multiple layers (in this case, a polygon layer showing the lines of Kenya)
```
> tm_shape(kenya_smooth) +
    tm_polygons() 
```

It is possible to use multiple shapes in one plot. For example, we can include the shape object for Kenya, and shape object for `counties`, which just has the counties of interest we specified earlier. For each `tmap` element, we can customize the appearance for example, filling the Busia and Kisumu counties with color.
```
> tm_shape(kenya_smooth) +
    tm_polygons() +
  tm_shape(counties) + 
    tm_polygons(col="tomato4") 
```

## 4e. Plotting a raster file
What if, instead of solid colored shapes we want to include gridded data (e.g. a raster image showing annual rainfall)? We can also use tmap to plot layers based on raster files. 

First, let's import an example raster file, `Sther_ENM.tif` available in this repository. This file shows habitat suitability values ouput from an environmental niche model for *Striga hermonthica* described in [Bellis *et al.* (2020)](https://www.pnas.org/content/117/8/4243.short). It has been cropped to approximately the extent of Kenya to make it a bit smaller.
```
> enm <- raster('Sther_ENM.tif')
> enm
> plot(enm)
```

We can see that the shape of the raster file is a rectangle. We can clip it to the shape of Kenya. But first we must project to the same coordinate reference system as our Kenya outline object.
```
> enm_utm <- projectRaster(enm, crs="+init=EPSG:32737 +proj=utm +zone=37 +south +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
```
Now, mask the raster using the outline of Kenya, setting all grid cells outside the polygon to 'NA'.

```
> enm_kenya <- mask(enm_utm, kenya_smooth)
> plot(enm_kenya)
```

And voila!
```
> tm_shape(enm_kenya) +
     tm_raster() +
 tm_shape(kenya_smooth) + 
     tm_borders() +
 tm_layout(legend.outside=T) +
 tm_scale_bar(width=0.1)
```

## 4f. Add points to show sampling locations
As another example of how we can include multiple shape objects in one plot, lets add points to show the sampling locations of *Striga hermonthica*.

First, download and read in the `Striga_GPS.csv` file (included in this repository). Convert to a SpatialPointsDataFrame. Don't forget to transform the coordinates to the UTM projection, as we did for our base mapping layers!
```
> points <- read.csv('~/Documents/GitHub/PEER_workshop/4_R_Spatial_Data/Striga_GPS.csv', header=T)
> points_df <- SpatialPointsDataFrame(cbind.data.frame(points$Lon, points$Lat),points, proj4string = CRS("+proj=longlat"))
> points_UTM<-spTransform(points_df, CRS("+init=EPSG:32737")) 
```
Then, add a `tm_shape` for the sampling location data, and a `tm_dots` layer to plot the points. An alternative to setting the layer to a single color is to map 'color' to one of the variables in the `points_UTM` object, in this case the column indicating elevation (in feet).  We can also add a scale bar!
```
> tm_shape(kenya_smooth) +
    tm_polygons() +
  tm_shape(counties) + 
    tm_polygons() +
  tm_shape(points_UTM) +
    tm_dots(col="Elevation_ft", size=0.15) +
  tm_scale_bar(width=0.2)
```

We can change the look of the plot in many ways. To see some of the possibilities, you can always pull up the documentation for a function with, e.g.
```
> ?tm_polygons
```
