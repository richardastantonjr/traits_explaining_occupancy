
##------------- Mapping ------------------
#libraries
library (tidyverse)
library(metafor)
library(maps)
library(mapdata)
library(rgdal)

load("./outputs/DR_metacommunity_richnesses.RData") 
coords <- read.csv("./data/LandUse.csv") 
LandUse <- read.csv("./data/LandUse.csv")
siteCovs <- siteCovs   

## one of the coords in the orginal file is extraneous and needs to be pulled
which(!(coords$GPSpt %in% samplingCovs$GPSpt )) ## #289 is not in the canonical list
coords <- coords[-289,]

## ---------------------------PROCESSING ------------------------------
siteCovs$site <- coords[, 1]  ## add a "site" column to the site covariates
siteCovs$site <- as.numeric(siteCovs$site) ## match data type with column w same name in "coords"

## subset coords to be siteCol (GPSpt), y_proj, x_proj
coords<- coords[,1:3]
colnames(coords)[1] <-"site"
#coords$y_proj<-coords$y_proj*-1     ## southern hemisphere

siteCovs <- cbind(coords, siteCovs)
siteCovs <- siteCovs[,c(1:3,5,7:10)]
## group survey points by grid and randomly select one point from each grid to map
grid_locations <- select(siteCovs, y_proj ,x_proj, LandUse, grid) %>% 
  group_by(grid) %>% 
  sample_n(1)
grid_locations <- data.frame(grid_locations)
grid_locations$y_proj <- grid_locations$y_proj*-1     ## southern hemisphere
xy <- grid_locations[,c(2,1)]    
## add the proper proj4string metadata from the gps and make a spdf
grid_locations<-SpatialPointsDataFrame(coords = xy, data = grid_locations,
                                       proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))


## create a study area map with each grid displayed, different symbols for each land use
tiff("./outputs/SwaziInSAfrica.tif", width = 4, height = 6, units = 'in',res=600, compression = "lzw")
par(mar=c(0, 0, 0, 0))
map(database = "world",regions=c("Mozambique", "South Africa", "Botswana", "Namibia", "Zambia", "Zimbabwe"), mar = rep(0, 4), 
    col = "moccasin", fill = T )
map(database = "world",regions = "Swaziland",add = T,fill = T,col = "olivedrab4")
dev.off()

tiff("./outputs/Swazi2014gridsZoomed.tif", width = 7, height =  10 , units = 'in',
     res = 600, compression = "lzw")
map(database = "world",regions = "Swaziland", xlim = c((min(xy$x_proj) - 0.3),
                (max(xy$x_proj) + 0.8)), mar = rep(0,4))
box()
points(as.matrix(grid_locations@coords), pch = 20, cex = 3,  
       col = c(3,5,4,7)[as.numeric(grid_locations$LandUse)])
map.scale(relwidth = 0.50, cex = 2) 
colors <- c(4, 3, 5, 7)
siteCovs_tiffs <- siteCovs  ## copy siteCovs and relevel LandUse for plotting purposes
siteCovs_tiffs$LandUse <- relevel(siteCovs$LandUse,"Protected")
legend("topright", # position
       legend = levels(plyr::revalue(siteCovs_tiffs$LandUse, 
                c("Protected" = "protected","ComPasture" = "pasture", 
                "Homestead" = "homestead", "SugarEstate" = "plantation"))),
       title = "Land Use",
       fill = colors,
       cex = 2.5,
       bty = "n") # border
dev.off()

## create shrub cover by land use figure
tiff("./outputs/shrub_LU.tif", width = 8.5, height = 3.5, units = 'in', res = 600,
     compression = "lzw")
par(oma = c(0, 0.9, 0, 0), cex.axis = 1.5)
boxplot(ShrubMean~LandUse, data = siteCovs_tiffs, outline = FALSE, las = 1, 5, cex.lab = 1.5, ylim = c(0, 100), 
        names= c("protected","pasture", "homestead", "plantation"), ylab = "Mean shrub cover (%)",
        col=c(4, 3, 5, 7))
abline(h = mean(siteCovs$ShrubMean),lty = "dotted")
stripchart(ShrubMean~LandUse, data = siteCovs_tiffs, vertical = TRUE,method = "jitter", 
           add = TRUE, pch = 1, cex = 0.50)
dev.off()
