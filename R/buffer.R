#' Return the neighbors of a reference cell 
#'
#' @description 
#' Extracts raster cells inside a buffer. Coordinates of the center of the buffer 
#' can be set by user argument (cell_id).
#' The world map (only displayed if plot is TRUE) represents the species richness
#' corresponding to the taxa.
#' Note that the taxa and associated max dispersion capabilities (buffer radius)
#' must be set as global options using the function `options()`.
#
#' @param plot a logical. If `TRUE` plot selected cells
#' 
#' @param buff_rad an integer : max disp distance (define the radius of the buffer in meters)
#'
#' @return A `data.frame` with three columns: `reference`, the reference 
#' community (cell id of the buffer center), `neighbor`, neighbor communities
#' (other cells inside the buffer), and `dist`, the distance between the 
#' reference cell and the neighbor.
#' 
#' @export
#'
#' @author Nicolas Casajus, \email{nicolas.casaju@@fondationbiodiversite.fr},
#'  


buffer <- function(cell_id,plot,buff_rad,distance=FALSE) {
  
   # cell_id <- 48903
   # plot=TRUE
   # buff_rad=600000
   
  ## Check inputs ----
  
  if (is.null(options()$"taxa")) {
    stop("Please set taxa in global options - ", 
         "`options(taxa = \"mammals\")`")
  }
  
 
  if (!file.exists(here::here("data","Loiseau2021", "grid_area_with_cell_ids.tif"))) {
    stop("Missing file: '", file.path("data", "Loiseau2021","grid_area_with_cell_ids.tif"), 
         "'")
  }
  
  if (!file.exists(here::here("data", "Loiseau2021", 
                              paste0(options()$"taxa","_richness.tif")))) {
    stop("Missing file: '", file.path("data", "Loiseau2021", 
                                      paste0(options()$"taxa", 
                                             "_richness.tif'")))
  }
  
  user_par <- graphics::par()$mar
  
  
  ## Import study area (grid) ----
  
    grd <- suppressWarnings(
      raster::raster(here::here("data", "Loiseau2021","grid_area_with_cell_ids.tif")))

  ## Convert to sf ----

  xy <- get_cell_coords(cell_id)[,2:3]
  #xy <- data.frame(x=focal_cell$x,y=focal_cell$y)
  xy_sf <- sf::st_as_sf(xy, coords = 1:2)
  xy_sf <- sf::st_set_crs(xy_sf, raster::projection(grd))
  
  
  ## Define buffer ----
  
  xy_buf <- sf::st_buffer(xy_sf, dist = buff_rad)
  
  
  ## Extract communities ----
  
  new_grd <- raster::crop(grd, xy_buf)
  new_grd <- raster::mask(new_grd, xy_buf)
  
  
  ## Plot ----
  
  if (plot) {
    
    ## Layer to map ----
    
    richness <- suppressWarnings(
      raster::raster(here::here("data", "Loiseau2021", 
                                paste0(options()$"taxa", "_richness.tif"))))
    
    richness <- raster::crop(richness, xy_buf)
    richness <- raster::mask(richness, xy_buf)
    
    color_s <- RColorBrewer::brewer.pal(name = "YlOrRd", 9)
    color_s <- grDevices::colorRampPalette(color_s)(255)
    
    graphics::par(mar = rep(1, 4))
    
    plot(sf::st_geometry(xy_buf), asp = 1, border = NA)
    
    raster::plot(richness, asp = 1, axes = FALSE, legend = FALSE,
                 legend.mar = 0, legend.width = 0, col = color_s, add = TRUE)
    
    plot(sf::st_geometry(xy_buf), add = TRUE, border = "white", lwd = 2)
    plot(sf::st_geometry(xy_buf), add = TRUE)
    
    graphics::points(xy, pch = 19)
    
    graphics::title(main = "Focal cell", family = "Arial")
    
  }
  
  graphics::par(mar = user_par)
  
  
  ## Return ----
  
  commnities <- as.integer(stats::na.omit(new_grd[]))
  commnities <- commnities[commnities != cell_id]
  
  commnities <- data.frame("reference" = rep(cell_id, length(commnities)),
                           "neighbor"  = commnities)
  
  if (distance==FALSE) {
    
    rbind(data.frame(reference=cell_id,neighbor=cell_id),commnities)
  } else {
    
    xy <- get_cell_coords(commnities$"neighbor")
    
    xy_sf <- sf::st_as_sf(xy, coords = 2:3)
    xy_sf <- sf::st_set_crs(xy_sf, raster::projection(grd))
    
    xy_ref <- get_cell_coords(commnities$"reference"[1])
    
    xy_ref_sf <- sf::st_as_sf(xy_ref, coords = 2:3)
    xy_ref_sf <- sf::st_set_crs(xy_ref_sf, sf::st_crs(xy_sf))
    
    d_ist <- as.matrix(as.numeric(sf::st_distance(xy_sf, xy_ref_sf)))
    
    xy_ref <- sf::st_drop_geometry(xy_sf)
    
    temp <- data.frame(commnities, "dist" = d_ist)
    
    rbind( data.frame(reference=cell_id,neighbor=cell_id,dist=0),temp)
    
  }
  
 
}

# 
# buffer <- function(cell_id,plot,buff_rad,distance=FALSE) {
#   
#   cell_id <- 10047
#   #plot=TRUE
#   buff_rad=600000
#   
#   ## Check inputs ----
#   
#   if (is.null(options()$"taxa")) {
#     stop("Please set taxa in global options - ", 
#          "`options(taxa = \"mammals\")`")
#   }
#   
#   
#   if (!file.exists(here::here("data", "grid_area_with_cell_ids.tif"))) {
#     stop("Missing file: '", file.path("data", "grid_area_with_cell_ids.tif"), 
#          "'")
#   }
#   
#   if (!file.exists(here::here("data", options()$"taxa", 
#                               paste0(options()$"taxa","_richness.tif")))) {
#     stop("Missing file: '", file.path("data", options()$"taxa", 
#                                       paste0(options()$"taxa", 
#                                              "_richness.tif'")))
#   }
#   
#   user_par <- graphics::par()$mar
#   
#   
#   ## Import study area (grid) ----
#   
#   grd <- suppressWarnings(
#     raster::raster(here::here("data", "grid_area_with_cell_ids.tif")))
#   
#   ## Convert to sf ----
#   
#   xy <- get_cell_coords(cell_id)[,2:3]
#   #xy <- data.frame(x=focal_cell$x,y=focal_cell$y)
#   xy_sf <- sf::st_as_sf(xy, coords = 1:2)
#   xy_sf <- sf::st_set_crs(xy_sf, raster::projection(grd))
#   
#   
#   ## Define buffer ----
#   
#   xy_buf <- sf::st_buffer(xy_sf, dist = buff_rad)
#   
#   
#   ## Extract communities ----
#   
#   new_grd <- raster::crop(grd, xy_buf)
#   new_grd <- raster::mask(new_grd, xy_buf)
#   
#   
#   ## Plot ----
#   
#   if (plot) {
#     
#     ## Layer to map ----
#     
#     richness <- suppressWarnings(
#       raster::raster(here::here("data", options()$"taxa", 
#                                 paste0(options()$"taxa", "_richness.tif"))))
#     
#     richness <- raster::crop(richness, xy_buf)
#     richness <- raster::mask(richness, xy_buf)
#     
#     color_s <- RColorBrewer::brewer.pal(name = "YlOrRd", 9)
#     color_s <- grDevices::colorRampPalette(color_s)(255)
#     
#     graphics::par(mar = rep(1, 4))
#     
#     plot(sf::st_geometry(xy_buf), asp = 1, border = NA)
#     
#     raster::plot(richness, asp = 1, axes = FALSE, legend = FALSE,
#                  legend.mar = 0, legend.width = 0, col = color_s, add = TRUE)
#     
#     plot(sf::st_geometry(xy_buf), add = TRUE, border = "white", lwd = 2)
#     plot(sf::st_geometry(xy_buf), add = TRUE)
#     
#     graphics::points(xy, pch = 19)
#     
#     graphics::title(main = "Focal cell", family = "Arial")
#     
#   }
#   
#   graphics::par(mar = user_par)
#   
#   
#   ## Return ----
#   
#   commnities <- as.integer(stats::na.omit(new_grd[]))
#   commnities <- commnities[commnities != cell_id]
#   
#   commnities <- data.frame("reference" = rep(cell_id, length(commnities)),
#                            "neighbor"  = commnities)
#   
#   if (distance==FALSE) {
#     
#     rbind(data.frame(reference=cell_id,neighbor=cell_id),commnities)
#   } else {
#     
#     xy <- get_cell_coords(commnities$"neighbor")
#     
#     xy_sf <- sf::st_as_sf(xy, coords = 2:3)
#     xy_sf <- sf::st_set_crs(xy_sf, raster::projection(grd))
#     
#     xy_ref <- get_cell_coords(commnities$"reference"[1])
#     
#     xy_ref_sf <- sf::st_as_sf(xy_ref, coords = 2:3)
#     xy_ref_sf <- sf::st_set_crs(xy_ref_sf, sf::st_crs(xy_sf))
#     
#     d_ist <- as.matrix(as.numeric(sf::st_distance(xy_sf, xy_ref_sf)))
#     
#     xy_ref <- sf::st_drop_geometry(xy_sf)
#     
#     temp <- data.frame(commnities, "dist" = d_ist)
#     
#     rbind( data.frame(reference=cell_id,neighbor=cell_id,dist=0),temp)
#     
#   }
#   
#   
# }
