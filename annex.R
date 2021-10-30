#############################################################
# Just a bunch of code snippets not used in the final product
#############################################################

# load and prep igc data ----

  lf <-  data.table(filenames=list.files("igc_aerothlon_2021"))
  lf[,athlete:=substring(gsub("\\..*","",filenames), 10)]
  lf[,I:=.I]
  igcs_raw <- list()
  for (i in 1:lf[,.N]) {
    #for (i in 1:6) {
    message(paste0(i, "/", lf[,.N]),"\r",appendLF=FALSE)
    flush.console()
    igcs_raw[[i]] <- geoviz::read_igc(paste0("igc_aerothlon_2021/",lf[i,filenames])) %>% data.table()
  }

  igc_prep <- function(igc, reduce_res = 1, dem, zscale) {
    igc[,time := as.ITime(time_char) ]
    alt_delta <- igc[,mean(altitude-altitude_pressure)] #pressure is more accurate and smooth but miscalibrated
    igc <- igc[time > as.ITime("9:00:00") & #start of the race
                 second(time) %% reduce_res == 0, #reduce resolution by 1 RoM
               .(time, lat, long, alt = altitude, alt_pressure = altitude_pressure,
                 alt_corr = altitude_pressure + alt_delta)]

    #elevate each igc so that none of them crosses the ground
    sp_gps <- sp::SpatialPoints(cbind(igc$long, igc$lat), proj4string = sp::CRS('+init=epsg:4326'))
    sp_gps <- sp::spTransform(sp_gps, sp::CRS(as.character(raster::crs(dem))))
    igc[,alt_ground := raster::extract(dem, sp_gps)]
    igc[,alt_floor := pmax(alt_pressure, alt_ground + 5)]

  }
  igcs <- lapply(igcs_raw, igc_prep, reduce_res = 1, dem = dem, zscale = zscale)
  
  igcs <- igcs[-c(2,3,25,36)] #DNF
  igcs <- igcs[-c(4,5,7,9,14,16,18:24,27,29:35)] #team
  
  saveRDS(igcs, "igcs.rds")
  igcs <- readRDS("igcs.rds")
  toc()
  

# Prepare Digital Elevation Model ----

  # alternative dem source  
  mapbox_key <- "pk.eyJ1IjoiYmFsdXgiLCJhIjoiY2tzeGxsb3JrMWd1NjJxbHNqMXQyYXVqYiJ9.k58ku_HsKNqneGzkEeNyrQ"
  dem <- mapbox_dem(lat, long, square_km, api_key = mapbox_key) 
  dem <- raster::raster("dem") #275k data points 
  
# Prepare Overlay Satellite Image -----
  tic()
  #geoviz ----
  # overlay_image <-
  #   slippy_overlay(dem,
  #                  image_source = "mapbox", image_type = "satellite",
  #                  #image_source = "stamen", image_type = "terrain",
  #                  #image_source = "stamen", image_type = "watercolor",
  #                  max_tiles = 500,
  #                  api_key = mapbox_key,
  #                  return_png = T,
  #                  png_opacity = 0.5)
  # saveRDS(overlay_image, "ol_mb_sat_500.rds")
  
  #rayshaderanimate ----
  # #library(rayshaderanimate)
  # devtools::load_all("E:/Mega Sync/R/rayshaderanimate")
  # #ez más logika szerint vágja körbe és csak a gpx-et, nekem most inkább az egész völgy kell
  # #bbox_arcgis <- get_bbox_from_gpx_table(igc, arcgis = TRUE)  
  
  # bbox_dem <- sp::bbox(dem)
  # bbox_arcgis <- list()
  # bbox_arcgis$p1 <- list(long = bbox_dem[1,1], lat = bbox_dem[2,1])
  # bbox_arcgis$p2 <- list(long = bbox_dem[1,2], lat = bbox_dem[2,2])
  # overlay_image <- get_image_overlay(bbox_arcgis = bbox_arcgis,
  #                                    output_file_loc = "ol_arcgis_test.png",
  #                                    flipflop = F,
  #                                    major_dim = 1000)
  
  #load ----
  message("Load overlay")
  #overlay_image <- readPNG("ol_arcgis_600.png") #major_dim 600 8MB 1M
  #overlay_image <- readPNG("ol_arcgis_1000.png") #major_dim 1000 21MB 2.7M
  if (!hq) overlay_image <- readRDS("ol_mb_sat_10.rds") #mapbox satelite maxtile 10
  # overlay_image <- readRDS("ol_mb_sat_100.rds") #mapbox satelite maxtile 100
  if ( hq) overlay_image <- readRDS("ol_mb_sat_500.rds") #mapbox satelite maxtile 500
  #overlay_image <- readRDS("overlay_image.rds") #stamen watercolor 10
  #overlay_image <- png::readPNG("overlay.png") #from mapporn_rsa.R (rayshaderanimate)
  
  #Optionally, turn mountainous parts of the overlay transparent ----
  # overlay_image <-
  #   elevation_transparency(overlay_image,
  #                          dem,
  #                          pct_alt_high = 0.5,
  #                          alpha_max = 0.9)
  toc()