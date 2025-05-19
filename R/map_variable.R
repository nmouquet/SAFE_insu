#' Map a variable (raster)
#'
#' @param data a `data.frame` with at least two columns. The first column must 
#' be raster cell ids (communities) and the first cell must be the reference
#' community. The variable to plot (argument `variable` must be in `data`).
#' 
#' @param variable a character of length 1. The column name of the variable to 
#' map.
#' 
#' @param center a logical. If `TRUE` show the center cell
#' 
#' @param col a color palette
#'
#' @return NULL
#' 
#' @export
#'
#' @examples
#' \dontrun{
#' 
#' xy <- get_cell_coords(60972)
#' datas <- select_communities(xy$x, xy$y)
#' 
#' ## Fictive data ----
#' cells <- c(datas[1, "reference"], datas[ , "neighbor"])
#' datas <- data.frame(cell_id = cells, cluster_1 = runif(length(cells)))
#' 
#' ## Map variable ----
#' map_variable(datas, "cluster_1", reference = FALSE)
#' }


map_variable <- function(data, variable, disp_max=disp_max,reference = FALSE, col = grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(255), ...) {
  
  # data=I_buffer[,-c(1:3)]
  # variable="Div_insurer"
  # reference=TRUE
  # col <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(name = "RdBu", 9)))(255)

  if (missing(data)) {
    stop("Argument 'data' is required")
  }
  
  if (!is.data.frame(data)) {
    stop("Argument 'data' must be a data frame")
  }
  
  if (missing(variable)) {
    stop("Argument 'variable' is required")
  }
  
  if (missing(reference)) {
    stop("Argument 'reference' is required")
  }
  
  if (!is.character(variable)) {
    stop("Argument 'variable' must be a character of length 1 (column name)")
  }
  
  if (length(variable) != 1) {
    stop("Argument 'variable' must be a character of length 1 (column name)")
  }
  
  if (!file.exists(here::here("data", "Loiseau2021","grid_area_with_cell_ids.tif"))) {
    stop("Missing file: '", file.path("data", "Loiseau2021","grid_area_with_cell_ids.tif"), 
         "'")
  }
  
  
  grd <- suppressWarnings(
    raster::raster(here::here("data", "Loiseau2021","grid_area_with_cell_ids.tif")))
  
  user_par <- graphics::par()$mar
  
  
  ## Center of area ----
  
  xy <- get_cell_coords(data[1, 1])
  
  
  ## Convert to sf ----
  
  xy_sf <- sf::st_as_sf(xy, coords = 2:3)
  xy_sf <- sf::st_set_crs(xy_sf, raster::projection(grd))
  
  
  ## Define buffer ----
  
  xy_buf <- sf::st_buffer(xy_sf, dist = disp_max)

  
  ## Rasterize variable to map ----
  
  ras <- grd
  ras[] <- NA
  
  xy <- get_cell_coords(data[ , 1])
  cells <- raster::cellFromXY(grd, xy[ , -1])
  
  ras[][cells] <- data[ , variable]

  
  ## Crop raster ----
  
  ras <- raster::crop(ras, xy_buf)
  ras <- raster::mask(ras, xy_buf)

  
  ## Plot ----
  
  graphics::par(mar = rep(1, 4))
  
  plot(sf::st_geometry(xy_buf), asp = 1, border = NA)
  
  raster::plot(ras, asp = 1, axes = FALSE, legend = FALSE, legend.mar = 0, 
               legend.width = 0, add = TRUE, col = col)
  
  plot(sf::st_geometry(xy_buf), add = TRUE, border = "white", lwd = 2)
  plot(sf::st_geometry(xy_buf), add = TRUE)
  
  xy <- get_cell_coords(data[1, 1])
  if (reference) graphics::points(xy[,2:3], pch = 19)
  
  graphics::title(main = variable, family = "Arial")
  
  graphics::par(mar = user_par)
  
  invisible(NULL)
}
