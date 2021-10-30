# 0. Settings ----

    # Set this 'high quality' switch to FALSE at first so you can play with the script on low resolution data 
    hq <- FALSE 
    
    # Adjust this once you installed ffmpeg
    ffmpeg_install_folder = "C:\\ffmpeg\\bin;"
    
    # Set your mapbox api key here (it is free but requires registration at https://account.mapbox.com/)
    mapbox_key <- "FIXME"

# 1. Setup ----
  
  library(tictoc);tic(); 
  library(rayshader)
  library(raster)
  library(geoviz)
  library(data.table)
  
  igcs <- readRDS("data/igcs_demo.rds")
  igccolors <- rainbow(length(igcs))
  
# 2. Prepare Digital Elevation Model (DEM) ----
  tic()
  
  lat_cent  <- 47.274559 #Lat of the race start / finish (center of the vis)
  lon_cent <- 13.318525 #Lon of the race start / finish (center of the vis)
  square_km <- 5 #Length of the scene's edge
  dem_file  <- paste0("data/dem",if(hq) "_hq")
    
  # DEM file is loaded if exists or downloaded and saved at first run
    if( paste0(dem_file,".gri") %>% file.exists() ) {
      
      message("loading DEM")
      dem <- raster::raster(dem_file)
      
    } else {
      
      message("downloading and saving DEM")
      dem <- geoviz::mapzen_dem(lat_cent, 
                                lon_cent, 
                                square_km, 
                                max_tiles = if(hq) 100 else 10)
      raster::writeRaster(x = dem, dem_file, overwrite = TRUE)
      
    }
      
    toc()
  
  # Downloaded DEM tiles cover a larger area so we crop them to fit our scene
    tic()
    message("cropping DEM")
    dem <- geoviz::crop_raster_square(dem,lat_cent,lon_cent,square_km)
    toc()

# 3. Calculate Elevation Matrix ----
    tic()
    message("Calculate elevation matrix") 
    elmat = matrix(
      raster::extract(dem, raster::extent(dem), method = 'bilinear'),
      nrow = ncol(dem),
      ncol = nrow(dem)
    )
    toc()
      
# 4. Prepare Overlay Satellite Image -----

    tic()
    oi_file  <- paste0("data/overlay_image",if(hq) "_hq", ".rds")
    
    # Overlay Image is loaded if exists or downloaded and saved at first run
    if( oi_file %>% file.exists() ) {
      
      message("loading Overlay Image")
      overlay_image <- readRDS(oi_file)
      
    } else {
      
      message("downloading and saving Overlay Image")
      overlay_image <- geoviz::slippy_overlay(
             dem,
             image_source = "mapbox",
             image_type = "satellite",
             max_tiles = if(hq) 500 else 10,
             api_key = mapbox_key,
             return_png = T,
             png_opacity = 0.5)
      saveRDS(overlay_image, oi_file)
    }
    
    toc()
  
# 5. Calculate Ambient Shadows ----
  
  tic()
  as_file <- paste0("data/ambientshadows",if(hq) "_hq", ".rds")
  
  if( as_file %>% file.exists() ) {
    
    message("loading Ambient Shadows data")
    ambientshadows <- readRDS(as_file)
    
  } else {
    
    message("Calculating and saving Ambient Shadows data")
    ambientshadows <- rayshader::ambient_shade(elmat)
    saveRDS(ambientshadows, as_file)
  }
  
  toc()

# 6. Prepare the Scene ----
  
  tic()
  sc_file <- paste0("data/scene",if(hq) "_hq", ".rds")
  sunangle <- 250 #rayshader parameter
  sunaltitude <- 75 #rayshader parameter
  zscale <- if (hq) 4 else 16 # horizontal vs vertical ratio of the scene
    #(approximation alternative: geoviz::raster_zscale())
  
  if( sc_file %>% file.exists() ) {
    
    message("loading the Scene")
    scene <- readRDS(sc_file)
    
  } else {
    
    message("Calculating and saving the Scene")
    
    #all functions here are from the rayshader package
    scene <- elmat %>%
      sphere_shade(sunangle = sunangle) %>%
      #add_water(detect_water(elmat), color = "lightblue") %>%
      add_shadow(ray_shade(elmat,
                           sunangle = sunangle,
                           sunaltitude = sunaltitude,
                           zscale = zscale,
                           multicore = FALSE), 
                 max_darken = 0.2) %>%
      add_shadow(ambientshadows, max_darken = 0.1) %>%
      add_shadow(lamb_shade(elmat,zscale = zscale, sunaltitude = 3), max_darken = 0.5) %>% 
      add_overlay(overlay_image, alphalayer = .8)
    
    saveRDS(scene, sc_file)
  }
  
  toc()
  
# 7. Check the Scene by plotting it ----

  tic()
  message("Plotting the Scene")
  
  # 2D plot
    #plot_map(scene)
  
  # Interactive 3D plot
  rayshader::plot_3d(
    scene,
    elmat,
    zscale = zscale,
    #baseshape = 'circle',
    zoom = 0.8,
    fov = 10,
    mouseMode = c("none", "trackball", "zoom", "none", "zoom"),
    windowsize = c(854,480),
    shadowcolor = 'grey10',
    background = 'lightskyblue1',
    triangulate = F
    )
  toc()
  
  tic()
  message("Add tracklogs to the Scene")

  for (i in 1:length(igcs)) {
    geoviz::add_gps_to_rayshader(
      dem,
      igcs[[i]]$lat,
      igcs[[i]]$long,
      igcs[[i]]$alt,
      line_width = 3,
      lightsaber = FALSE,
      alpha = .5,
      colour = igccolors[i],
      zscale = zscale,
      ground_shadow = TRUE
    )
  }
  toc()


# 8. Create animation frames -----

  #' This function will create the frames of a rayshader animation with tracklogs
  #'
  #' @param dem digital elevation model
  #' @param elmat elevation matrix
  #' @param scene a rayshader scene
  #' @param frames_folder arbitrary folder to store the frames
  #' @param zscale horizontal vs vertical ratio of the scene
  #' @param theta camera rotation start and end in degrees (vector of 2) 
  #' @param zoom camera zoom start and end (vector of 2) 
  #' @param phi camera azimuth start and end in degrees (vector of 2) 
  #' @param duration duration of the animation in seconds
  #' @param intro_portion tracklogs will start at this point of the duration
  #' @param fps frames per seconds
  #' @param windowsize resolution of the animation (vector of 2)
  #' @param igcs list of tracklogs 
  #' @param igccolors color palette for the tracklogs
  #' @param ground_shadow set to TRUE to render the shadows of the tracklogs
  #' @param line_width width of the rendered tracklogs
  #' @param alpha transparency of the tracklogs
  #' @param hq "high quality switch" - if FALSE the resolution of the tracklogs is reduced
  #'
  #' @return vector of the location and names of the frames (to be feed into render_rayshader_animantion())
  create_rayshader_animantion_frames <- function(
    dem,
    elmat,
    scene,
    frames_folder,
    zscale = 1,
    theta = c(-160,20), 
    zoom = c(1.5, 0.25), 
    phi = c(90, 30), 
    duration = 3, 
    intro_portion = .2,
    fps = 24,
    windowsize = c(854,480), #480p
    igcs = NULL,
    igccolors = 'blue',
    ground_shadow = FALSE,
    line_width = 2,
    alpha = .8,
    hq = FALSE) {
    
    start <- Sys.time(); print(start);
    
    message("Setting animation parameters")
    
    dir.create(frames_folder, showWarnings = FALSE)
    unlink(paste0(frames_folder,"*")) 
    frame_names <-  paste0(frames_folder,"/frame_", 1:(duration * fps), ".png") 
    vapply(frame_names, file.remove, logical(1))
    num_of_frames <- length(frame_names)
    
    theta_angles <- theta[2] + (theta[1] - theta[2]) * 1/(1 + seq(0, 8, length.out = num_of_frames)^2)
    zoom_scale   <- zoom[2]  + (zoom[1]  - zoom[2] ) * 1/(1 + seq(0, 8, length.out = num_of_frames)^2)
    phi_angles   <- phi[2]   + (phi[1]   - phi[2]  ) * 1/(1 + seq(0, 8, length.out = num_of_frames)^2)
    
    message("Reducing igc resolution to fit the video length")
    
    if (!is.null(igcs)) {
      
      #first portion of video is just intro, no igc rendered
      num_of_frames_intro <- floor(num_of_frames*intro_portion)
      
      #get the longest of the igc-s to fit its length to the video duration
      igcs_max_nrow <- lapply(igcs,nrow) %>% unlist() %>% max()
      
      #if hq == FALSE then reduce resolution of igc-s by dropping rows
      if (!hq) {
        for(l in 1:length(igcs)){
          igc <- igcs[[l]]
          igc[,filter := .I %% round(igcs_max_nrow/(num_of_frames - num_of_frames_intro)) - 1]
          igc<-igc[filter == 0]
          igc[,filter := NULL]
          igcs[[l]] <- igc
        }
      }
    }
    
    message("Plotting the scene")
    
    rayshader::plot_3d(
      scene,
      elmat,
      zscale = zscale,
      #baseshape = 'circle',
      windowsize = windowsize,
      shadowcolor = 'grey10',
      background = 'lightskyblue1',
      triangulate = F
    )  
    
    message("Generating frames:")
    
    for (i in 1:length(frame_names)) {
      
      #status message for monitoring the process
      time_spent <- difftime(Sys.time(),start,units = "mins")
      message(paste0(" Frame: ", i, "/", length(frame_names), 
                     " // Time spent: ", round(time_spent[[1]], 2), " ", units(time_spent), "  "),
              "\r",appendLF = ifelse(i == length(frame_names), TRUE, FALSE))
      
      #set camera position
      render_camera(theta = theta_angles[i], zoom = zoom_scale[i], phi = phi_angles[i], fov = 50)
      
      #render igcs after the intro
      if (!is.null(igcs)) {
        if (i > num_of_frames_intro) {
          for (l in 1:length(igcs)) {
            
            igc <- igcs[[l]]
            
            j <- i - num_of_frames_intro
            
            #for high quality output full resolution igcs are rendered
            if (hq) j <- j * round(igcs_max_nrow/(num_of_frames - num_of_frames_intro))
            
            #if an igc reached its end keep its full length until video finishes
            j <- min(j,igc[,.N])
            
            #render igc
            geoviz::add_gps_to_rayshader(
              dem,
              igc$lat[1:j],
              igc$long[1:j],
              igc$alt[1:j],
              line_width = line_width,
              lightsaber = FALSE,
              alpha = alpha,
              colour = igccolors[l],
              zscale = zscale,
              ground_shadow = ground_shadow
            ) 
          }
        }
      }
      
      #save frame  
      render_snapshot(filename = frame_names[i])
      
      #remove rendered igcs from the scene after each frame
      tmp <- rgl::ids3d()
      rgl::pop3d(id = tmp[tmp[[2]]=="linestrip",1])
    }
    
    message(paste0("Total Time spent: ", round(time_spent[[1]], 2), " ", units(time_spent)))
    
    rgl::rgl.close()
    
    return(frame_names)
  }
  
  frame_names  <- create_rayshader_animantion_frames(
    elmat = elmat,
    scene = scene,
    dem=dem,
    frames_folder = paste0(getwd(),"/tmp"),
    zscale = zscale,
    duration = if (hq) 30 else 1, 
    igcs = igcs,
    igccolors=igccolors,
    hq = hq
  )
  
# 9. Retouch each frame with Photoshop ----

  # This step was done outside the R universe!
  
  # I created a Photoshop batch process to automatically retouch each frame 1-by-1:
    # Lightning effect for the background
    # Adjustment of Levels
    # Added Aerothlon logo

# 10. Render animation ----
  
  #' Create video from individual frames using ffmpeg
  #'
  #' @param frame_names vector of the location and names of the frames 
  #' (returned by create_rayshader_animantion_frames())
  #' @param outputfile name of the rendered video
  #' @param overwrite to overwrite the outputfile
  #' @param ffmpeg_install_folder location on your computer where ffmpeg is installed
  #'
  #' @return name of the rendered video
  render_rayshader_animantion <- function(frame_names,
                                          outputfile = "animation.mp4",
                                          overwrite = TRUE,
                                          ffmpeg_install_folder = "C:\\ffmpeg\\bin;") {
    message("Rendering Video")
    all_paths <-  "tmp.txt"
    writeLines(con = all_paths,
               paste0( "file '" , gsub("\\\\","/", frame_names),"'"))
    Sys.setenv("PATH" = paste0(Sys.getenv("PATH"), ffmpeg_install_folder))
    system(intern = TRUE,
           paste0("ffmpeg ",
                  ifelse(overwrite, "-y", "-n"),
                  " -f concat -r 24 -safe 0 -i \"",
                  all_paths,
                  "\" -vf \"fps=24,format=yuv420p\" ", outputfile))
    
    return(outputfile)
  }
  
  render_rayshader_animantion(frame_names, ffmpeg_install_folder = ffmpeg_install_folder)

   
   
  

  