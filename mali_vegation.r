install.packages("ncdf4")
install.packages("raster")
#install.packages("rgdal")
#install.packages("ggplot2")
install.packages('stars')
library(ncdf4)
library(raster)
#library(rgdal)
#library(ggplot2)
library(RColorBrewer)
library(lattice)
library(stars)
library(dplyr)
library(ncdf4.helpers)
library(tibble)

nc_data <- nc_open("~/Documents/GitHub/c_gls_FCOVER_201606130000_GLOBE_PROBAV_V1.5.1.nc", )
print(nc_data)

lon <- ncvar_get(nc_data, "lon")
lon <- nc_data$dim$lon$vals
lat <- ncvar_get(nc_data, "lat", verbose = F)
t <- ncvar_get(nc_data, "time")
ncatt_get(nc_data, "time", "units")
head(t)
dim(t)

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

filename = "~/Documents/GitHub/MOD13A2.A2016353.h16v07.061.2021363143829.hdf"
sd <- gdal_subdatasets(filename)
read_stars("~/Documents/GitHub/MOD13A2.A2016097.h17v08.061.2021348033000.hdf")
r <- read_stars(filename)

setwd("~/Desktop/C0184943/")
filenames <- paste("~/Desktop/C0184943", list.files("~/Desktop/C0184943/"), list.files(list.files("~/Desktop/C0184943/"))[1:27], sep = "/")
coords <- data.frame(Longitude = c(8.51395, 8.373600, 8.399230, 8.468310,
                                   8.372850, 8.308260, 8.258717, 8.551460,
                                   8.446560, 8.314414, 8.424590, 8.397900,
                                   8.396860, 8.316880)*-1,
                     Latitude = c(11.985750, 12.05229, 12.04576, 11.96599, 
                                  12.09628, 12.15466, 12.073598, 11.98836,
                                  11.95466, 12.068577, 12.03865, 11.96737, 
                                  12.13444, 12.14454))
ndvi_villages <- matrix(0, nrow=14, ncol=27)
for (i in 1:27) {
  nc_data <- nc_open(filenames[i])
  # print(nc_data)
  lon <- ncvar_get(nc_data, "lon")
  lat <- ncvar_get(nc_data, "lat", verbose = F)
  t <- ncvar_get(nc_data, "time")
  fillvalue <- ncatt_get(nc_data, "NDVI", "_FillValue")
  ndvi.array <- ncvar_get(nc_data, 'NDVI') # store the data in a 3-dimensional array
  ndvi.array[ndvi.array == fillvalue$value] <- NA
  
  # grid <- expand.grid(lon=lon, lat=lat)
  # cutpts <- c(-50,-40,-30,-20,-10,0,10,20,30,40,50)
  # levelplot(crs.array ~ lon * lat, data=grid, at=cutpts, cuts=11, pretty=T, 
  #           col.regions=(rev(brewer.pal(10,"RdBu"))))
  proj4 <- nc.get.proj4.string(nc_data, "NDVI")
  nc_close(nc_data)
  rm(nc_data)
  
  r <- raster(t(ndvi.array), xmn=min(lon), xmx=max(lon), ymn=min(lat), 
              ymx=max(lat), crs=CRS(proj4))
  #plot(r)

  values <- extract(r, coords)
  ndvi_villages[,i] <- values
  rm(r, ndvi.array)
  print(i)
  # attributes(nc_data$var)
  # dim(ndvi.array)
  # x <- SpatialPoints(coords)
  # lat_lon <- as.matrix(expand.grid(lat,lon))
  # ndvi_vector <- as.vector(ndvi.array)
  # head(ndvi_vector)
  # ndvi_df <- data.frame(cbind(lat_lon, ndvi_vector))
  # ndvi_df <- data.frame(latitude = lat_lon[,1],
  #                       longitude = lat_lon[,2],
  #                       ndvi = ndvi_vector)
  # rm(nc_data, ndvi.array, lat_lon, ndvi_vector)
  # head(ndvi_df)
  # dim(ndvi_df)
  # # ndvi_df <- na.omit(ndvi_df)
  # tail(ndvi_df)
  # ndvi_df |>
  #   filter(latitude > 11.5 & latitude < 12.5 & longitude > -8.7 & longitude < -8.0) -> ndvi_filtered
  # 
  # which(abs(ndvi_filtered$latitude-12.14454) == min(abs(ndvi_filtered$latitude-12.14454))
  #       & abs(ndvi_filtered$longitude + 8.316880) == min(abs(ndvi_filtered$longitude + 8.316880)))
  # ndvi_filtered[4856,]
  }

ndvi_villages <- data.frame(ndvi_villages)
mean_nona <- function(x) {
  mean(x, na.rm=T)
}
mean_ndvi_village <- apply(ndvi_villages, MARGIN=1, FUN = mean_nona)

mali <- read.csv("~/Documents/GitHub/atsb_working_code/DB monthly natural sugar and ASB feeding/day 2-Table 1.csv")[-(64:65),1:16]
mali$dyed_fraction <- mali$females.ASB.positive/mali$TOTAL.Sample.females.Day.2
mali$total_asb_positive <- mali$females.ASB.positive 
mali$total_sampled <- mali$TOTAL.Sample.females.Day.2 
mali$total_catch <- mali$CDC.total.females.day.2
mali$days <- (mali$month-1)*30+mali$Day

par(las=1, mfrow=c(1,1), mar=c(8,4,4,1)+0.1)
mali |>
  group_by(Village) |>
  summarise(dyed_fraction=weighted.mean(dyed_fraction,total_sampled),
            total_catch=sum(total_catch)) -> mali_grouped
par(las=2)
plot(mali_grouped$dyed_fraction,
     col=brewer.pal(7,"Set2")[factor(mali_grouped$Village)],
     pch=19,
     cex=2*mali_grouped$total_catch/mean(mali_grouped$total_catch),
     frame.plot = F,
     ylim = c(0,1),
     ylab = "Dyed fraction",
     xaxt = "n",
     xlab = "")
axis(1, at=1:7, labels=mali_grouped$Village)
grid()
villages <- mali_grouped$Village
cluster_quantiles <- matrix(0,nrow = length(villages), ncol=2)
for (i in 1:length(villages)) {
  this_village <- villages[i]
  mali_filtered <- filter(mali, Village==this_village)
  replicates <- matrix(0,nrow=5000,ncol=nrow(mali_filtered))
  for (j in 1:5000) {
    replicates[j,] <- sample(x = mali_filtered$total_asb_positive/mali_filtered$total_sampled, 
                             replace = TRUE, 
                             size = nrow(mali_filtered), 
                             prob = mali_filtered$total_sampled)
  }
  cluster_means <- apply(replicates, MARGIN = 1, FUN = mean)
  cluster_quantiles[i,] <- quantile(cluster_means, probs = c(0.025,0.975))
}
colnames(cluster_quantiles) <- c("feed_rate_lower", "feed_rate_upper")
mali_grouped <- cbind(mali_grouped, cluster_quantiles)

arrows(1:7,
       mali_grouped[,4],
       1:7,
       mali_grouped[,5],
       code=0,
       lwd=2,
       col = brewer.pal(7,"Set2")[factor(mali_grouped$Village)])
vegetation <- c(mean_ndvi_village[4],
                mean_ndvi_village[1],
                mean_ndvi_village[5],
                mean_ndvi_village[3],
                mean_ndvi_village[2],
                mean_ndvi_village[6],
                mean_ndvi_village[7])
mali_grouped$ndvi <- vegetation

par(las=1, mar=c(5,4,4,1)+0.1)
plot(mali_grouped$ndvi,
     mali_grouped$dyed_fraction,
     col=brewer.pal(7,"Set2")[factor(mali_grouped$Village)],
     cex=mali_grouped$total_catch/mean(mali_grouped$total_catch),
     frame.plot = F,
     ylab = "Dyed fraction",
     xlab = "Vegetation index",
     pch=19,
     xlim=c(0,0.6),
     ylim=c(0,0.6))
grid()
arrows(mali_grouped$ndvi,
       mali_grouped[,4],
       mali_grouped$ndvi,
       mali_grouped[,5],
       code=0,
       lwd=2,
       col=brewer.pal(7,"Set2")[factor(mali_grouped$Village)])
ndvi_villages$village <- c("Balandougou", "Madina", "Korea", "Balala",
                           "Cissebougou", "Nianguanabougou", "Trekrouba",
                           "Krekrelo", "Sirakele", "Trekrou", "Farabale", 
                           "Kignele", "Tiko", "Sambadani")

villages <- mali_grouped$Village
cluster_quantiles <- matrix(0,nrow = length(villages), ncol=2)
for (i in 1:length(villages)) {
  this_village <- villages[i]
  ndvi_filtered <- filter(ndvi_villages, village==this_village)
  ndvi_filtered <- as.double(ndvi_filtered)[!is.na(as.double(ndvi_filtered))]
  replicates <- matrix(0,nrow=5000,ncol=length(ndvi_filtered))
  for (j in 1:5000) {
    replicates[j,] <- sample(ndvi_filtered, 
                             replace = TRUE, 
                             size = length(ndvi_filtered))
  }
  cluster_means <- apply(replicates, MARGIN = 1, FUN = mean)
  cluster_quantiles[i,] <- quantile(cluster_means, probs = c(0.025,0.975))
}
colnames(cluster_quantiles) <- c("ndvi_lower", "ndvi_upper")
mali_grouped <- cbind(mali_grouped, cluster_quantiles)
arrows(mali_grouped$ndvi_lower,
       mali_grouped$dyed_fraction,
       mali_grouped$ndvi_upper,
       mali_grouped$dyed_fraction,
       code=0,
       lwd=2,
       col=brewer.pal(7,"Set2")[factor(mali_grouped$Village)])

t_ndvi <- t(ndvi_villages[1:7,1:27])
colnames(t_ndvi) <- ndvi_villages$village[1:7]
t_ndvi <- data.frame(t_ndvi)
t_ndvi <- t_ndvi[ , c(4,1,5,3,2,6,7)]
t_ndvi <- add_column(t_ndvi, Day = seq(10,270,10)+90, .before = 1)
par(las = 1)
plot(
  NA,
  xlim = c(100,400),
  ylim = c(0.2, 0.8),
  xlab = "Days",
  ylab = "NDVI",
  frame.plot = FALSE
)
grid()
for (i in 1:7) {
  lines(
    t_ndvi$Day,
    t_ndvi[, 1+i],
    lwd = 2,
    col = brewer.pal(7,"Set2")[i]
  )
}
par(las = 2, mar = c(8,4,4,1)+0.1)
plot(NA, 
     frame.plot = F,
     ylim = c(0.2,0.8),
     ylab = "NDVI",
     xaxt = "n",
     xlab = "",
     xlim = c(1,7))
grid()
points(
  1:7,
  mali_grouped$ndvi,
  col=brewer.pal(7,"Set2")[factor(mali_grouped$Village)],
  cex=2,
  pch = 19
)
arrows(
  1:7,
  mali_grouped$ndvi_lower,
  1:7,
  mali_grouped$ndvi_upper,
  col = brewer.pal(7,"Set2")[factor(mali_grouped$Village)],
  code = 0,
  lwd = 2
)
axis(1, at=1:7, labels=mali_grouped$Village)

par(las=1, mar=c(5,4,4,1)+0.1)
plot(mali$days,
     mali$dyed_fraction,
     col=brewer.pal(7, "Set2")[factor(mali$Village)],
     cex=mali$total_sampled/mean(mali$total_sampled),
     pch=19,
     frame.plot=F,
     ylim=c(0,1),
     xlim=c(0,400),
     xlab = "Day",
     ylab = "Dyed fraction")
grid()
