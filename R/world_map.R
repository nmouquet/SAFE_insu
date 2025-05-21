#'
#' @author Nicolas Casajus, \email{nicolas.casaju@@fondationbiodiversite.fr},
#'

world_map <- function(
  id = NA,
  color_s = RColorBrewer::brewer.pal(name = "YlOrRd", 9)
) {
  #id=NA
  color_s = RColorBrewer::brewer.pal(name = "YlOrRd", 9)
  ## Layer to map ----

  richness <- suppressWarnings(
    raster::raster(here::here(
      "data",
      "Loiseau2021",
      paste0(options()$"taxa", "_richness.tif")
    ))
  )

  ## Plot richness gradient

  user_par <- graphics::par()$mar
  on.exit(graphics::par(mar = user_par))

  ### Color palette ----

  color_s <- grDevices::colorRampPalette(color_s)(255)

  graphics::par(mar = rep(1, 4))

  raster::plot(
    richness_robinson,
    asp = 1,
    axes = FALSE,
    legend = FALSE,
    legend.mar = 0,
    legend.width = 0,
    col = color_s
  )

  if (!is.na(id)) {
    for (i in id) {
      ## Center of area ----

      xy <- get_cell_coords(i)

      if (!is.na(xy[1, "x"])) {
        ## Convert to sf ----

        xy_sf <- sf::st_as_sf(xy, coords = 2:3)
        xy_sf <- sf::st_set_crs(xy_sf, raster::projection(richness))

        ## Define buffer ----

        xy_buf <- sf::st_buffer(xy_sf, dist = options()$"dispersion_max")

        plot(
          sf::st_geometry(xy_buf),
          asp = 1,
          border = "black",
          col = "#00000044",
          add = TRUE
        )
      }
    }
  }

  invisible(NULL)
}

world_map_robin <- function(
  id = NA,
  color_s = RColorBrewer::brewer.pal(name = "YlOrRd", 9)
) {
  ## Layer to map ----

  richness <- suppressWarnings(
    terra::rast(here::here(
      "data",
      "Loiseau2021",
      paste0(options()$"taxa", "_richness.tif")
    ))
  )
  # Define the Robinson projection
  robinson_crs <- "ESRI:54030" # EPSG code for Robinson projection

  # Reproject the raster to Robinson projection using terra's project function
  richness_robinson <- terra::project(richness, robinson_crs)

  ## Plot richness gradient

  user_par <- graphics::par()$mar
  on.exit(graphics::par(mar = user_par))

  ### Color palette ----

  color_s <- grDevices::colorRampPalette(color_s)(255)

  graphics::par(mar = rep(1, 4))

  terra::plot(
    richness_robinson,
    asp = 1,
    axes = FALSE,
    legend = FALSE,
    legend.mar = 0,
    legend.width = 0,
    col = color_s
  )
  if (!is.na(id)) {
    for (i in id) {
      ## Center of area ----

      xy <- get_cell_coords(i)

      if (!is.na(xy[1, "x"])) {
        ## Convert to sf ----

        xy_sf <- sf::st_as_sf(xy, coords = 2:3)
        xy_sf <- sf::st_set_crs(xy_sf, raster::projection(richness))

        ## Define buffer ----

        xy_buf <- sf::st_buffer(xy_sf, dist = options()$"dispersion_max")

        plot(
          sf::st_geometry(xy_buf),
          asp = 1,
          border = "black",
          col = "#00000044",
          add = TRUE
        )
      }
    }
  }
  invisible(NULL)
}

#' Map a variable (World map)
#'
#' @param data a `data.frame` with two columns. The first column is the cell id
#'   and the second the variable to map.
#'
#' @param palname a string with RColorBrewer color pal used (ex. YlORRd, BuGn )
#'
#' @return No return value.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' vals <- get(load(here::here("results" , "buffer_mammals.RData")))
#' world_map2(vals[ , c("cell_focal", "RI_Hui")])
#' }
#'
#'
world_map2 <- function(data, color_s) {
  #data <- birds_insurance[ , c("cell_focal", "SS")]
  #color_s=color_s

  if (!is.data.frame(data)) {
    stop("Argument 'data' must be a data.frame")
  }

  if (nrow(data) == 0) {
    stop("Argument 'data' must have at least 1 row")
  }

  if (ncol(data) != 2) {
    stop(
      "Argument 'data' must have exactly 2 columns (cell id and variable ",
      "to map)"
    )
  }

  colnames(data)[1] <- "old_cell_id"

  ## Layer to map ----

  richness <- suppressWarnings(
    raster::raster(here::here(
      "data",
      "Loiseau2021",
      paste0(options()$"taxa", "_richness.tif")
    ))
  )

  ras <- richness

  ## Reset values ----

  ras[] <- NA

  ## Retrieve cell ids ----

  grd <- suppressWarnings(
    raster::raster(here::here(
      "data",
      "Loiseau2021",
      "grid_area_with_cell_ids.tif"
    ))
  )

  values <- data.frame(
    "real_cell_id" = 1:raster::ncell(grd),
    "old_cell_id" = grd[]
  )
  values <- na.omit(values)

  data <- merge(
    values,
    data,
    by.x = "old_cell_id",
    by.y = "old_cell_id",
    all = FALSE
  )

  ras[][data$"real_cell_id"] <- data[, 3]

  ## Plot variable ----

  user_par <- graphics::par()$mar
  on.exit(graphics::par(mar = user_par))

  ### Color palette ----

  graphics::par(mar = rep(1, 4))

  raster::plot(
    ras,
    asp = 1,
    axes = FALSE,
    legend = FALSE,
    legend.mar = 0,
    legend.width = 0,
    col = color_s
  )

  invisible(NULL)
}

world_map2_robin <- function(data, color_s) {
  #data <- birds_insurance[ , c("cell_focal", "SS")]
  #color_s=color_s

  if (!is.data.frame(data)) {
    stop("Argument 'data' must be a data.frame")
  }

  if (nrow(data) == 0L) {
    stop("Argument 'data' must have at least 1 row")
  }

  if (ncol(data) != 2L) {
    stop(
      "Argument 'data' must have exactly 2 columns (cell id and variable ",
      "to map)"
    )
  }

  colnames(data)[1] <- "old_cell_id"

  ## Layer to map ----

  richness <- terra::rast(
    here::here("data", "Loiseau2021", paste0(options()$"taxa", "_richness.tif"))
  )

  ras <- richness

  ## Reset values ----

  ras[] <- NA

  ## Retrieve cell ids ----

  grd <- terra::rast(
    here::here("data", "Loiseau2021", "grid_area_with_cell_ids.tif")
  )

  values <- data.frame(1:terra::ncell(grd), grd[])

  colnames(values) <- c("real_cell_id", "old_cell_id")

  values <- na.omit(values)

  data <- merge(
    values,
    data,
    by.x = "old_cell_id",
    by.y = "old_cell_id",
    all = FALSE
  )

  ras[][data$"real_cell_id"] <- data[, 3]

  ras <- terra::project(
    x = ras,
    y = paste0(
      "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 ",
      "+datum=WGS84 +units=m +no_defs"
    )
  )

  ## Plot variable ----

  user_par <- graphics::par()$mar
  on.exit(graphics::par(mar = user_par))

  ### Color palette ----

  graphics::par(mar = rep(1, 4))

  terra::plot(x = ras, asp = 1, axes = FALSE, legend = FALSE, col = color_s)

  invisible(NULL)
}
