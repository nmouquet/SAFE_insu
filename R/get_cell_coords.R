#' Get cells coordinates
#'
#' @param id a vector of integers. Cells ids
#'
#' @return A `data.frame` with three columns: `id`, `x`, and `y`. The number of
#' rows corresponds to the number of cells in input (argument `id`).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' get_cell_coords(75097)
#' get_cell_coords(c(60972, 75097, 77154))
#' }
#'
#' @author Nicolas Casajus, \email{nicolas.casaju@@fondationbiodiversite.fr},
#'

get_cell_coords <- function(id) {
  #id=48903

  if (missing(id)) {
    stop("Argument 'id' (cell id) is required")
  }

  if (!is.numeric(id)) {
    stop("Argument 'id' must be an integer of length > 0")
  }

  grd <- suppressWarnings(
    raster::raster(here::here(
      "data",
      "Loiseau2021",
      "grid_area_with_cell_ids.tif"
    ))
  )

  values <- grd[]

  cell_id <- unlist(lapply(1:length(id), function(x) {
    cell <- which(values == id[x])
    if (length(cell) == 0) {
      cell <- NA
    }
    names(cell) <- id[x]
    cell
  }))

  data.frame(cell_id = names(cell_id), raster::xyFromCell(grd, cell_id))
}
