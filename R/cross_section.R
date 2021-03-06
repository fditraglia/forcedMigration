#' Cross-section dataset
#'
#' Cross-section data on a large number of variables for each municipality in
#' Colombia. For now, you only need to worry about the ones we document here.
#'
#' @format A tibble with 831 rows and 96 variables:
#' \describe{
#'   \item{municipality}{municipality code}
#'   \item{lon_mean}{Degrees of longitude for municipality centroid (negative = west)}
#'   \item{lat_mean}{Degrees of latitude for municipality centroid (positive = north)}
#'   \item{department}{Department (equivalent of state) in which a municipality is located}
#'   \item{region}{Geographic region (within department) in which a municipality is located}
#' }
#' @source Fill in later!
"cross_section"
