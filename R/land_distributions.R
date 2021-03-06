#' Land distribution dataset
#'
#' Cross-section data on the pre-conflict land distribution of each municipality
#' in Colombia.
#'
#' @format A list with 856 elements, each of which is a dataframe for a particular
#' municipality. List items are named with municipality coedes. Rows correspond
#' to landholding "bins" (a discretized land sitribution) designated in hectares
#' as indicated in the row names. The columns are as follows:
#' \describe{
#'   \item{total_land}{Total amount of land owned by families in this landholding "bin"}
#'   \item{n_families}{Number of families in this landholding "bin"}
#'   \item{mean_land}{Average landholdings (in hectares) for families in this landholding "bin"}
#'   \item{frac_families}{Fraction of families in the municipality with landholdings in this bin}
#' }
#' @source Fill in later!
"land_distributions"
