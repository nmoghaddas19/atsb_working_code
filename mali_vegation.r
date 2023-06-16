install.packages("ncdf4")
install.packages("raster")
install.packages("rgdal")
install.packages("ggplot2")
library(ncdf4)
library(raster)
library(rgdal)
library(ggplot2)
library(RColorBrewer)
library(lattice)
nc_data <- nc_open("~/Documents/GitHub/c_gls_FCOVER_201606130000_GLOBE_PROBAV_V1.5.1.nc", )
print(nc_data)

lon <- ncvar_get(nc_data, "lon")
lon <- nc_data$dim$lon$vals
lat <- ncvar_get(nc_data, "lat", verbose = F)
t <- ncvar_get(nc_data, "time")
ncatt_get(nc_data, "time", "units")
head(t)
dim(t)

HDF5_USE_FILE_LOCKING=FALSE
crs.array <- ncvar_get(nc_data, 'FCOVER') # store the data in a 3-dimensional array
dim(ndvi.array) 

nc_close()

nc_data_2 <- nc_open("~/Documents/GitHub/c_gls_FCOVER_201606230000_GLOBE_PROBAV_V1.5.1.nc")
print(nc_data_2)
lon <- ncvar_get(nc_data_2, "lon")
lat <- ncvar_get(nc_data_2, "lat", verbose = F)
t <- ncvar_get(nc_data_2, "time")
ncatt_get(nc_data, "time", "units")
head(t)
dim(t)
fillvalue <- ncatt_get(nc_data_2, "FCOVER", "_FillValue")
HDF5_USE_FILE_LOCKING=FALSE
fcover.array <- ncvar_get(nc_data_2, 'FCOVER') # store the data in a 3-dimensional array
fcover.array[fcover.array == fillvalue$value] <- NA
dim(fcover.array)

length(na.omit(as.vector(fcover.array)))
image(lon,lat,crs.array, col = brewer.pal(10,"RdBu"))
grid <- expand.grid(lon=lon, lat=lat)
cutpts <- c(-50,-40,-30,-20,-10,0,10,20,30,40,50)
levelplot(crs.array ~ lon * lat, data=grid, at=cutpts, cuts=11, pretty=T, 
          col.regions=(rev(brewer.pal(10,"RdBu"))))
nc_close(nc_data_2)

r <- raster(t(fcover.array), xmn=min(lon), xmx=max(lon), ymn=min(lat), 
            ymx=max(lat), crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))
plot(r)

extract(r, SpatialPoints(cbind(-12.565, 8.923)), method = 'simple')

install.packages('stars')
library(stars)
filename = "~/Documents/GitHub/MOD13A2.A2016353.h16v07.061.2021363143829.hdf"
sd <- gdal_subdatasets(filename)
read_stars("~/Documents/GitHub/MOD13A2.A2016097.h17v08.061.2021348033000.hdf")
r <- read_stars(filename)



