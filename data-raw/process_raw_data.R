library(dplyr)
library(igraph)
library(purrr)
library(tidyr)
library(furrr)
library(sf)
library(spdep)
library(geosphere)
library(rmapshaper)

#------------------------------- Load raw data
# The cross section dataset contains information for only those municipalities
# with land distribution data
cross_section_raw <- haven::read_dta("data-raw/Cross_Section_final_RAW_1076.dta")

# The panel dataset contains information for municipalities that we cannot
# use in our models because we lack land distribution data for them
panel_raw <- haven::read_dta("data-raw/panel-raw.dta")

survey_raw <- haven::read_dta("data-raw/survey-raw.dta")

municipalities_and_regions_raw <- readr::read_csv("data-raw/DANE_municipality_codes_and_regions.csv")

#-------------------------------------------------------------------------------
# Shapefile for Colombian municipalities
#    source:   <https://data.humdata.org/dataset/cod-ab-col>
#    filename: "Col Administrative Divisions Shapefiles.zip"
#-------------------------------------------------------------------------------

# Read in shapefile.
CO_shape <- st_read('data-raw/CO-shapefiles/col_admbnda_adm2_mgn_20200416.shp')


#------------------- Merge department/province names
municipalities_and_regions <- municipalities_and_regions_raw |>
  rename(department = number_dept, municipality = Muni_code) |>
  mutate(municipality = as.numeric(municipality))

cross_section <- left_join(cross_section_raw, municipalities_and_regions)
panel <- left_join(panel_raw, municipalities_and_regions)

rm(municipalities_and_regions_raw, cross_section_raw, panel_raw,
   municipalities_and_regions)

#------------------------------- Remove two islands: 88001 and 88564
CO_islands <- c('88001', '88564')

survey_raw |>
  filter(municipality_origin %in% CO_islands) #None in the survey

cross_section <- cross_section |>
  filter(!(municipality %in% CO_islands))

panel <- panel |>
  filter(!(municipality %in% CO_islands))

CO_shape <- CO_shape |>
  mutate(municipality = as.numeric(stringr::str_remove(ADM2_PCODE, 'CO'))) |>
  filter(!(municipality %in% CO_islands))

rm(CO_islands)

#-------------------------------- Construct List of Neighbors
# construct list of neighbors: neighbors gives *indices* corresponding to
# rows of CO_shape while neighbor_codes gives *municipality codes*, i.e. the
# correponding values of the column municipality from CO_shape
neighbors <- poly2nb(CO_shape)
neighbor_codes <- map(neighbors, \(x) CO_shape$municipality[x])
names(neighbor_codes) <- CO_shape$municipality

#-------------------------------------- Select and rename panel variables
panel <- panel |>
  filter(year %in% 1996:2008) |>  # Only use panel data from 1996-2008
  select(municipality,
         year,
         V_cum = ac_vcivilparas_UR, # Cumulative paramilitary violence
         V_flow = vcivilparas_UR, # Flow of paramilitary violence
         placeboV_cum = ac_vcivilnotparas_UR, # Cumulative non-paramilitary violence
         placeboV_flow = vcivilnotparas_UR, # Flow of non-paramilitary violence
         D_AS = desplazados_paras_AS, # Four measures of displacement
         D_RUV = desplazados_total_RUV,
         D_CEDE = desplazados_CEDE,
         D_JYP = desplazados_JYP,
         pop1993 = Total_Poblacion1993)

#--------------------------- Compare municipalities in neighbors list vs panel
panel_munis <- panel |>
  pull(municipality) |>
  unique()

munis <- as.numeric(stringr::str_remove(CO_shape$ADM2_PCODE, 'CO'))

# All municipalities in panel are also in adjacency matrix
identical(sum(panel_munis %in% munis), length(panel_munis))

#---------------------- Construct average violence in neighboring municipalities

get_V_flow_neighbors <- function(muni_code) {
  V_flow_neighbors <- panel |>
    filter(municipality %in% neighbor_codes[[paste(muni_code)]]) |>
    select(municipality, year, V_flow) |>
    pivot_wider(names_from = municipality, values_from = V_flow) |>
    select(-year) |>
    apply(1, mean)

  out <- panel |>
    filter(municipality == muni_code) |>
    select(municipality, year)

  # A very small number of municipalities have neighbors, but not neighbors for
  # which we observe violence data, in which case V_flow_neighbors equals numeric(0)
  if(length(V_flow_neighbors) > 0){
    out$V_flow_neighbors <- V_flow_neighbors
  } else {
    out$V_flow_neighbors <- NA_real_
  }
  return(out)
}

plan(multisession, workers = 8)

V_flow_neighbors <- future_map(panel_munis, get_V_flow_neighbors) |>
  list_rbind()


panel <- panel |>
  left_join(V_flow_neighbors)

rm(V_flow_neighbors, get_V_flow_neighbors, munis, panel_munis)

#-------------------------------- Clean land distribution characteristics

#------------------------------------------------------------------------------
# The column landless_families in from Cross_section.dta was constructed by
# Camilo's original RA. There are some discrepancies in this column which we
# resolve as follows:
#
# 1. Calculate the number of landowners by summing the counts in each column
#   "landown_*". This information comes from the Cadastral census.
#
# 2. The number of families (from the Census) is given by num_families.
#
# 3. Call m the number of landowners and n the number of families. If we assume
#    that there is at most one landowner per family, then m should be less than
#    or equal to n. But in some cases m > n.
#
# 4. To address this we suppose that M ~ Poisson(mu) and N ~ Poisson(lambda).
#    We know that mu < lambda and want to estimate the ratio r = mu / lambda.
#    This is the share of landholding families. The share of landless families
#    is omega = (1 - r). We calculate the Bayesian posterior mean under the
#    Poisson model, assuming conditional independence given the parameters, and
#    flat priors. This immediately implies the posterior mean for omega, the
#    share of landless. Crucially *this will never be exactly zero*. This is
#    reasonable as an estimate and also avoids problems later on when we form
#    log-odds ratios with the share of landless in the denominator. If you don't
#    like the Bayesian interpretation, you can think of this as a data-informed
#    way of adding "epsilon" to the share of landless families.
#-------------------------------------------------------------------------------

get_omega <- function(m, n) {
# Estimate the share of landless families (omega) from the number of landowners
# (m) and the number of families (n) in the municipality, using the procedure
# described above. Assumes that m > 0 and n > 1 but allows for NA values in
# either. Vectorized.
  log_I_denom <- pbeta(0.5, n, m + 1, lower.tail = FALSE, log = TRUE)
  log_I_num <- pbeta(0.5, n - 1, m + 2, lower.tail = FALSE, log = TRUE)
  r <- exp(log(m + 1) - log(n - 1) + (log_I_num - log_I_denom))
  1 - r
}

cross_section <- cross_section |>
  rowwise(municipality) |>
  mutate(n_landowners = sum(c_across(starts_with('landown_')))) |>
  ungroup() |>
  rename(n_families_census = num_families) |>
  mutate(omega = get_omega(m = n_landowners, n = n_families_census)) |>
  select(-landless_families)


#-------------------------------------------------------------------------------
# We observe a discretized land distribution: we know the total amount of land
# and the total number of landholders within each of the following "bins"
# measured in hectares
#-------------------------------------------------------------------------------
# 0
# (0,1)
# [1,3)
# [3,5)
# [5,10)
# [10,15)
# [15, 20)
# [20, 50)
# [50,100)
# [100, 200)
# [200, 500)
# [500, 1000)
# >=2000
#-------------------------------------------------------------------------------
landowner_bins <- cross_section |>  # Number of families in each landholding "bin"
  select(starts_with('landown_')) |>
  rename(`(0,1)` = landown_less_1,
         `[1,3)` = landown_1_3,
         `[3,5)` = landown_3_5,
         `[5,10)` = landown_5_10,
         `[10,15)` = landown_10_15,
         `[15,20)` = landown_15_20,
         `[20,50)` = landown_20_50,
         `[50,100)` = landown_50_100,
         `[100,200)` = landown_100_200,
         `[200,500)` = landown_200_500,
         `[500,1000)` = landown_500_1000,
         `[1000,2000)` = landown_1000_2000 ,
         `>=2000` = landown_2000_plus)

landtotal_bins <- cross_section |>  # Total land in each landholding "bin"
  select(starts_with('tot_land_')) |>
  rename(`(0,1)` = tot_land_1,
         `[1,3)` = tot_land_2,
         `[3,5)` = tot_land_3,
         `[5,10)` = tot_land_4,
         `[10,15)` = tot_land_5,
         `[15,20)` = tot_land_6,
         `[20,50)` = tot_land_7,
         `[50,100)` = tot_land_8,
         `[100,200)` = tot_land_9,
         `[200,500)` = tot_land_10,
         `[500,1000)` = tot_land_11,
         `[1000,2000)` = tot_land_12,
         `>=2000` = tot_land_13)

#------------------------------------------------------------------------------

make_land_dist <- function(i) {
  rbind(total_land = landtotal_bins[i, ],
        n_landowners = landowner_bins[i, ]) |>
    t() |>
    as.data.frame() |>
    mutate(mean_land = if_else(n_landowners > 0, total_land / n_landowners, 0),
           frac_landowners = n_landowners / sum(n_landowners))
}

land_distributions <- map(seq_len(nrow(cross_section)), make_land_dist)
names(land_distributions) <- cross_section$municipality
rm(make_land_dist, landowner_bins, landtotal_bins)

# Function that computes summary statistics of the land distribution for use in
# our LINEAR IV models. The non-linear IV models rely on the entire land distribution
# stored in land_distributions

get_land_statistics <- function(x) {
  c('H1' = sum(asinh(x$mean_land) * x$frac_landowners),
    'H2' = sum(asinh(x$mean_land)^2 * x$frac_landowners))
}

land_statistics <- map(land_distributions, get_land_statistics)
land_statistics <- as.data.frame(do.call(rbind, land_statistics))
land_statistics$municipality <- as.numeric(rownames(land_statistics))
rownames(land_statistics) <- NULL
rm(get_land_statistics)

cross_section <- cross_section |>
  left_join(y = land_statistics) |>
  mutate(omegaC = 1 - omega)  # It's convenient to have (1 - omega) as a separate variable

rm(land_statistics)

# Cumulative violence in the final year of the panel, i.e. total violence, can
# be viewed as a cross-sectional characteristic. Merge this with the other
# cross_section data

merge_me <- panel |>
  filter(year == max(year)) |>
  select(municipality, V_cum)

cross_section <- cross_section |>
  inner_join(merge_me) |>
  rename(V_total = V_cum) # Better to call it V_total: total violence from 1996:2008

rm(merge_me)

#------------------------ Flag municipalities w/ land data and >=100 landholders

# The municipalities in cross_section are those for which we have land dist data
identical(sum(!is.na(cross_section$H1)), nrow(cross_section))

cross_section <- cross_section |>
  mutate(at_least_100_landowners = n_landowners >= 100)

munis_100_plus <- cross_section |>
  filter(at_least_100_landowners) |>
  pull(municipality)

panel <- panel |>
  mutate(has_land_data = municipality %in% cross_section$municipality,
         at_least_100_landowners = municipality %in% munis_100_plus)

rm(munis_100_plus)


#-------------------------------------Clean Displacement Data

#-------------------------------------------------------------------------------
# Some displacement measures are zero for all municipalities during a given year.
# These are not "true" zeros, but rather missing observations: some agencies did
# not choose to report in a given year:
#-------------------------------------------------------------------------------
# AS      did not record in: 1996, 2012
# CEDE    did not record in: 1996, 2010, 2011, 2012
#-------------------------------------------------------------------------------
# Replace these zeros with NA.
#-------------------------------------------------------------------------------
# Technically not necessary to include 2010:2012 since we've already subsetted
# to remove those years.
panel[panel$year %in% c(1996, 2012),]$D_AS <- NA
panel[panel$year %in% c(1996, 2010:2012),]$D_CEDE <- NA


#-------------------------------------------------------------------------------
# For all measures except JYP, there are some municipalities in which a single
# measure reports zero total displacement over the sample period, while all other
# measures report positive displacement. These "zeros" are really NAs.
#-------------------------------------------------------------------------------
disagreement <- panel |>
  select(municipality, year, starts_with('D_')) |>
  group_by(municipality) |>
  summarise(across(.cols = starts_with('D_'),
                   .fns = sum, na.rm = TRUE)) |>
  rowwise() |>
  filter(sum(c_across(starts_with('D_')) == 0) == 1) |>
  ungroup()


AS_replace <- disagreement |>
  filter(D_AS == 0) |>
  pull(municipality)
panel[panel$municipality %in% AS_replace,]$D_AS <- NA

RUV_replace <- disagreement |>
  filter(D_RUV == 0) |>
  pull(municipality)
panel[panel$municipality %in% RUV_replace,]$D_RUV <- NA

CEDE_replace <- disagreement |>
  filter(D_CEDE == 0) |>
  pull(municipality)
panel[panel$municipality %in% CEDE_replace,]$D_CEDE <- NA

rm(disagreement, AS_replace, RUV_replace, CEDE_replace)

#--------------------------------------- Create lagged violence flow
panel <- panel |>
  group_by(municipality) |>
  mutate(lag_V_flow = lag(V_flow, order_by = year)) |>
  ungroup()

#--------------------------- Re-scale displacement measures

# What share of total RUV and JYP displacement was reported in 1996?
# Average these two shares
share_1996 <- panel |>
  select(year, D_RUV, D_JYP) |>
  group_by(year) |>
  summarize(RUV_sum = sum(D_RUV, na.rm = TRUE),
            JYP_sum = sum(D_JYP, na.rm = TRUE)) |>
  mutate(RUV_share = RUV_sum / sum(RUV_sum),
         JYP_share = JYP_sum / sum(JYP_sum)) |>
  filter(year == 1996) |>
  select(RUV_share, JYP_share) |>
  rowMeans()

panel <- panel |>
  mutate(across(c(D_RUV, D_JYP), # Normalize total displacement to 6 million
         list(norm = ~ 6e6 * . / sum(., na.rm = TRUE)))) |>
  # Because AS and CEDE did not report in 1996 we don't normalize them to 6 million.
  # Instead we normalize them to a *slightly smaller value* constructed using
  # share_1996
  mutate(across(c(D_AS, D_CEDE),
                list(norm = ~ (1 - share_1996) * 6e6 * . / sum(., na.rm = TRUE)))) |>
  rowwise() |>
  mutate(D_med = median(c_across(ends_with('_norm')), na.rm = FALSE),
         D_avg = mean(c_across(ends_with('_norm')), na.rm = FALSE),
         D_med_any = median(c_across(ends_with('_norm')), na.rm = TRUE),
         D_avg_any = mean(c_across(ends_with('_norm')), na.rm = TRUE)) |>
  ungroup() |>
  mutate(D_avg_no_AS = (D_RUV + D_CEDE + D_JYP) / 3,
         D_avg_no_RUV = (D_AS + D_CEDE + D_JYP) / 3,
         D_avg_no_CEDE = (D_AS + D_RUV + D_JYP) / 3,
         D_avg_no_JYP = (D_AS + D_RUV + D_CEDE) / 3) |>
  # After taking medians of the normalized displacement variables (sum to 6 million) the result
  # will no longer sum to 6 million. Renormalize so it does
  mutate(across(c(starts_with('D_med'), starts_with('D_avg')), ~ 6e6 * . / sum(., na.rm = TRUE)))

rm(share_1996)


# Covariates ------------------------------------------------------------------

cross_section <- cross_section |>

  rename(mines = mine_titles_90, oil = oil_prod_98) |>

  mutate(mines = if_else(mines > 0, 1, 0),
         drug_routes = if_else(drug_routes > 1, 1, 0),
         oil = if_else(oil > 0, 1, 0),
         bureaucracy = local_bureaucracy_95 + state_bureaucracy_95,
         elec_comp = 0.5 * (ratio_votes_dif_alc + ratio_votes_dif_asa),
         offices = rowSums(across(c(notaries_95,
                                    banks_95,
                                    savings_95,
                                    health_centers_95,
                                    health_posts_95,
                                    schools_95,
                                    libraries_95,
                                    fire_95,
                                    jails_95,
                                    culture_95,
                                    instruments_95,
                                    tax_95,
                                    police_95,
                                    police_posts_95,
                                    courts_95,
                                    telecom_95,
                                    post_95,
                                    agrarian_95,
                                    hospitals_95)))) |>

  select(-starts_with('tot_land'), -starts_with('landown_'), -gini,
         -convexity_L, -starts_with('owners_'), -ends_with('_95'),
         -starts_with('ratio_votes'), -n_families_census,
         -contains('_signal'), -contains('educ_'), -dum_2001_land)


# Spatial stuff ----------------------------------------------------------------

# Read in ARCGIS-generated list of municipalities intersecting road network
roads <- readr::read_csv('data-raw/real_roads.csv') |>
  mutate(municipality = stringr::str_remove(ADM2_PCODE, 'CO'),
         municipality = as.numeric(municipality))


# Read and clean forest & elevation data --------------------------------------

forests <- readr::read_csv("data-raw/Forest_Master_50.csv") |>
  mutate(share_forested = extent_2000_ha / area_ha,
         is_forested = share_forested > 0.5) |>
  select(ADM2_PCODE, is_forested) |>
  mutate(municipality = stringr::str_remove(ADM2_PCODE, 'CO'),
         municipality = as.numeric(municipality)) |>
  select(-ADM2_PCODE)

# Municipality 70523 (Palmito, Sucre) appears twice in forests, once with
# is_forested TRUE and again with is_forested FALSE. There are no forests here!
forests <- forests |>
  filter(!((municipality == 70523) & (is_forested)))

# Elevation *difference*
elevation <- readr::read_csv("data-raw/elevationsFinal.csv") |>
  select(-`...1`) |>  # first column is a row number starting from zero
  mutate(municipality = stringr::str_remove(ADM2_PCODE, 'CO'),
         municipality = as.numeric(municipality)) |>
  select(-ADM2_PCODE) |>
  arrange(municipality)


# Create dataframe of geographic information ----------------------------------
geography <- cross_section |>
  select(municipality, ruggedness, slope, elevation) |>
  left_join(elevation) |>
  left_join(forests) |>
  mutate(has_road = if_else(municipality %in% roads$municipality, 1, 0),
         is_forested = 1 * is_forested)

rm(elevation, forests, roads)

centroid_coords <- CO_shape |>
  st_centroid() |>  # compute lat and lon of centroids
  st_coordinates() # extract lat and lon of centroids as *matrix*

get_neighbor_distances <- function(i) {
  if(all(neighbors[[i]] == 0)) {
    return(NA) # No neighbors? No distances!
  } else {
    own_coords <- centroid_coords[i, ]
    neighbor_coords <- centroid_coords[neighbors[[i]], ]
    distGeo(own_coords, neighbor_coords) / 1000 # convert meters to km
  }
}
neighbor_distances <- map(1:length(neighbors), get_neighbor_distances)

get_pairwise_distances <- function(i) {
# only return neighbors whose index is greater than i: eliminate "repeat" pairs
  i_neighbors <- neighbors[[i]]
  i_neighbor_distances <- neighbor_distances[[i]]
  keep <- i < i_neighbors
  if(any(keep)) {
    data.frame(i, j = i_neighbors[keep], dist = i_neighbor_distances[keep])
  } else {
    # return null if no neighbors (entry in neighbors is zero) or all repeats
    NULL
  }
}


pairwise_distances <- map(2:length(neighbors), get_pairwise_distances) |>
  list_rbind()
pairwise_distances <- tibble(pairwise_distances)

rm(get_neighbor_distances, get_pairwise_distances, centroid_coords)


# Impute missing values --------------------------------------------------------

# There are 41 missing values for the variables elevation_difference and
# is_forested from geography. We impute these with an average of the same
# variables for *neighboring* municipalities. If a variable is missing for all
# neighboring municipalities, we fill in the overall mean of the variable for
# geography.


# Create list of neighboring municipalities with municipality *codes* rather
# than row numbers from CO_shape
muni_codes <- as.numeric(stringr::str_remove(CO_shape$ADM2_PCODE, 'CO'))
neighbors_codes <- map(neighbors, function(x) muni_codes[x])
names(neighbors_codes) <- muni_codes

# geography only has data for 1049 municipalities while CO_shape has data for
# 1122. Add the missing municipalities to geography
geography <- geography |>
  full_join(tibble(municipality = muni_codes))


impute_missing_var <- function(var_name) {
  muni_missing_var <- geography[is.na(geography[, var_name]), 'municipality', drop = TRUE]
  municipality_indices <- which(muni_codes %in% muni_missing_var)
  neighbor_list <- neighbors_codes[municipality_indices]
  f <- function(x) {
    neighbor_geography <- pull(geography[geography$municipality %in% x, var_name])
    mean(neighbor_geography, na.rm = TRUE)
  }
  out <- map(neighbor_list, f)
  out <- do.call(rbind, out)
  out <- data.frame(row.names(out), out[,1])
  names(out) <- c('municipality', paste0('impute_', var_name))
  row.names(out) <- NULL
  out$municipality <- as.numeric(out$municipality)
  # If no neighbors, then there's an NaN in the second column. Replace with the
  # overall average
  out[is.na(out[, 2]), 2] <- mean(geography[, var_name, drop = TRUE], na.rm = TRUE)
  return(tibble(out))
}

# Impute missing values for all variables *except* municipality (first column)
# and store the result in a list
imputed_geography <- map(names(geography)[-1], impute_missing_var)

imputed_geography <- imputed_geography |>
  reduce(full_join, by = 'municipality')


# Last step of the imputation, following the logic of the comment, but
# for all the imputed variables, not just impute_forest
geography <- geography |>
  full_join(imputed_geography) |>
  mutate(ruggedness = if_else(!is.na(ruggedness), ruggedness, impute_ruggedness)) |>
  mutate(slope = if_else(!is.na(slope), slope, impute_slope)) |>
  mutate(elevation = if_else(!is.na(elevation), elevation, impute_elevation)) |>
  mutate(elevation_difference = if_else(!is.na(elevation_difference),
                                       elevation_difference,
                                       impute_elevation_difference)) |>
  mutate(is_forested =  if_else(!is.na(is_forested), is_forested, impute_is_forested)) |>
  mutate(has_road = if_else(!is.na(has_road), has_road, impute_has_road)) |>
  select(-starts_with('impute_'))


rm(impute_missing_var, imputed_geography, muni_codes, neighbors_codes)

cross_section <- cross_section |>
  full_join(geography)

# Construct dataset of all unique *pairs* of neighboring municipalities --------

# Create geography_i by renaming the columns of geography to have an '_i'
# at the end. Do the same for geography_j. Then merge with pairwise distances.
geography_i <- geography |>
  rename_with(~ stringr::str_c(., '_i'))

geography_j <- geography |>
  rename_with(~ stringr::str_c(., '_j'))

pairwise_distances <- pairwise_distances |>
  mutate(municipality_i = stringr::str_remove(CO_shape$ADM2_PCODE[i], 'CO'),
         municipality_j = stringr::str_remove(CO_shape$ADM2_PCODE[j], 'CO'),
         municipality_i = as.numeric(municipality_i),
         municipality_j = as.numeric(municipality_j)) |>
  select(-i, -j) |>
  left_join(geography_i) |>
  left_join(geography_j)

rm(geography_i, geography_j, geography)

link_stats <- pairwise_distances |>
  mutate(forest = is_forested_i + is_forested_j,
         elevation = (elevation_difference_i + elevation_difference_j) / 2,
         ruggedness = (ruggedness_i + ruggedness_j) / 2,
         slope = (slope_i + slope_j) / 2,
         elevation_diff = abs(elevation_i - elevation_j)) |>
  select(municipality_i, municipality_j, dist, forest, elevation, slope,
         ruggedness, elevation_diff)

pca <- prcomp(~ forest + elevation + ruggedness + slope + elevation_diff,
              data = link_stats, scale = TRUE, na.action = na.omit)


# Extract first PC and "rotate" so that higher values indicate less
# favorable terrain (higher "friction"). Note that this is mean zero
print(pca)

link_stats$PC1_geography <- -1 * pca$x[,1]

rm(pca, pairwise_distances)

# "Effective distance" is defined as:
#    dist * friction
# The friction is constructed from PC1 as follows:
#     z = (PC1 - min(PC1)) / (max(PC1) - min(PC1))
#     friction = 1 + (max_friction - 1) * z
#
# This gives a friction of 1 at the minimum of PC1 and a fraction of max_friction
# at the maximum of PC1. We use four values of max_friction: 1, 2, 2.5, and 3.
# Notice that max_friction = 1 means that PC1 plays no role: this is distance
# "as the crow flies" from one municipality centroid to the next.

link_stats <- link_stats |>
  mutate(z = (PC1_geography - min(PC1_geography)) / diff(range(PC1_geography)))


# Dijkstra's algorithm ---------------------------------------------------------

# First construct a graph from our dataframe of link_statistics
munigraph <- graph_from_data_frame(link_stats, directed = FALSE)

# The "Epicenter" of paramilitary violence: Valencia, Cordoba
epicenter <- '23855'
# Unweighted graph, no weights -> return "number of hops" between nodes
hops <- distances(munigraph, v = epicenter) # to argument defaults to V(graph)


get_dist_friction <- function(max_friction) {
  dist <- edge_attr(munigraph, 'dist')
  z <- edge_attr(munigraph, 'z')
  friction <- 1 + (max_friction - 1) * z
  my_weights <- dist * friction
  distances(munigraph, v = epicenter, weights = my_weights)
}

max_friction <- c(1, 2, 2.5, 3)
d_geography <- lapply(max_friction, get_dist_friction)
d_geography <- do.call(rbind, d_geography)
d_geography <- as.data.frame(t(d_geography))
names(d_geography) <- paste0('d_geography', max_friction)
d_geography$municipality <- as.numeric(rownames(d_geography))

dist_to_epicenter <- d_geography
dist_to_epicenter$d_hops <- drop(hops)
dist_to_epicenter <- as_tibble(dist_to_epicenter)
dist_to_epicenter <- dist_to_epicenter |>
  relocate(municipality) |>
  rename(d_crow = d_geography1)

# Calculate percentiles of distance measures (except hops)
dist_to_epicenter <- dist_to_epicenter |>
  mutate(across(c('d_crow', starts_with('d_geography')), ~ (rank(.x) / length(.x)),
                .names = 'p_{.col}'))

rm(get_dist_friction, epicenter, max_friction, d_geography, hops,
   link_stats, neighbor_distances, munigraph)

cross_section <- cross_section |>
  full_join(dist_to_epicenter)

rm(dist_to_epicenter)

# Simplify CO_shape polygons and remove islands --------------------------------

CO_shape <- ms_simplify(CO_shape) # defaults to retaining 5% of points

#-------------------------------------------------------------------------------
# Create abandoned land dataset.

# Note: There was an intermediate hand-cleaning step that was
# omitted between merging raw digitized with AttributeTableFinal.csv and input
# here to remove duplicates; the abandoned_land_handcleaned.csv file is the
# result of that cleaned merge.
#-------------------------------------------------------------------------------


usethis::use_data(CO_shape, overwrite = TRUE)
usethis::use_data(cross_section, overwrite = TRUE)
usethis::use_data(land_distributions, overwrite = TRUE)
usethis::use_data(neighbors, overwrite = TRUE)
usethis::use_data(neighbor_codes, overwrite = TRUE)
usethis::use_data(panel, overwrite = TRUE)

# clean up
rm(list = ls())


