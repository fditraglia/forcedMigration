library(dplyr)
library(purrr)
library(magrittr)
library(tidyr)
library(parallel)
library(sf)

#------------------------------- Load raw data
# The cross section dataset contains information for only those municipalities
# with land distribution data
cross_section_raw <- haven::read_dta("data-raw/Cross_Section_final_RAW_1076.dta")

# The panel dataset contains information for municipalities that we cannot
# use in our models because we lack land distribution data for them
panel_raw <- haven::read_dta("data-raw/panel-raw.dta")

survey_raw <- haven::read_dta("data-raw/survey-raw.dta")

municipalities_and_regions_raw <- readr::read_csv("data-raw/DANE_municipality_codes_and_regions.csv")

#------------------------------------------------------------------------------
# NOTE: we should redo everything involving the adjacency matrix to ensure that
# it agrees with the one we construct from the shapefiles below!
#------------------------------------------------------------------------------
adjacency_matrix_raw <- readxl::read_excel("data-raw/Adjacency_Matrix.xlsx")

#------------------- Merge department/province names
municipalities_and_regions <- municipalities_and_regions_raw %>%
  rename(department = number_dept, municipality = Muni_code) %>%
  mutate(municipality = as.numeric(municipality))

cross_section <- left_join(cross_section_raw, municipalities_and_regions)
panel <- left_join(panel_raw, municipalities_and_regions)

rm(municipalities_and_regions_raw, cross_section_raw, panel_raw)

#------------------------------- Remove two islands: 88001 and 88564
CO_islands <- c('88001', '88564')

survey_raw %>%
  filter(municipality_origin %in% CO_islands) #None in the survey

cross_section %<>% # Assignment pipe!
  filter(!(municipality %in% CO_islands))

panel <- panel %>%
  filter(!(municipality %in% CO_islands))

adjacency_matrix <- adjacency_matrix_raw %>%
  filter(!(codes %in% CO_islands)) # Note that we don't have to remove
                                   # columns 'border_with_88001' etc. since
                                   # these two municipalities are not in the
                                   # adjacency matrix: 1117 versus 1119 obs.

rm(adjacency_matrix_raw)

#-------------------------------------- Select and rename panel variables
panel %<>% # Assignment pipe!
  filter(year %in% 1996:2008) %>% # Only use panel data from 1996-2008
  select(municipality,
         year,
         V_cum = ac_vcivilparas_UR, # Cumulative paramilitary violence
         V_flow = vcivilparas_UR, # Flow of paramilitary violence
         placeboV_cum = ac_vcivilnotparas_UR, # Cumulative non-paramilitary violence
         placeboV_flow = vcivilnotparas_UR, # Flow of non-paramilitary violence
         D_AS = desplazados_paras_AS, # Five measures of displacement
         D_CODHES = desplazados_CODHES,
         D_RUV = desplazados_total_RUV,
         D_CEDE = desplazados_CEDE,
         D_JYP = desplazados_JYP,
         popn1993 = Total_Poblacion1993)


#--------------------------------- Construct list of neighboring municipalities
munis <- adjacency_matrix %>%
  pull(codes)

adjacency_matrix %<>% # Assignment pipe!
  select(-codes) %>%
  as.matrix()

rownames(adjacency_matrix) <- colnames(adjacency_matrix) <- munis
neighbors <- apply(adjacency_matrix, 1, function(row) munis[row == 1])

rm(adjacency_matrix)

#--------------------------- Compare municipalities in adjacency matrix vs panel

panel_munis <- panel %>%
  pull(municipality) %>%
  unique()

# All municipalities in panel are also in adjacency matrix
identical(sum(panel_munis %in% munis), length(panel_munis))

#------------------------------------------------------------------------------
# NOTE: we may need to redo this now that we've created an adjacency matrix
# directly from the shapefiles below!
#---------------------- Construct average violence in neighboring municipalities

get_V_flow_neighbors <- function(muni_code) {
  V_flow_neighbors <- panel %>%
    filter(municipality %in% neighbors[[paste(muni_code)]]) %>%
    select(municipality, year, V_flow) %>%
    pivot_wider(names_from = municipality, values_from = V_flow) %>%
    select(-year) %>%
    apply(1, mean)

  out <- panel %>%
    filter(municipality == muni_code) %>%
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

V_flow_neighbors <- mclapply(panel_munis, get_V_flow_neighbors, mc.cores = 8)
V_flow_neighbors <- do.call(rbind, V_flow_neighbors)


panel %<>% # Assignment pipe!
  left_join(V_flow_neighbors)

rm(V_flow_neighbors, get_V_flow_neighbors, munis, panel_munis)

#-------------------------------- Clean land distribution characteristics

#------------------------------------------------------------------------------
# The column landless_families in from Cross_section.dta was constructed by
# Camilo's original RA. There are some discrepancies in this column. Our
# procedure for constructing a measure of the number of landless families in a
# municipality from the underlying census data is as follows. We assume, absent
#evidence to the contrary, that there is at most one landowner per family.
#------------------------------------------------------------------------------
# 1. Calculate n_landowners by summing counts in each column "landown_*"
# 2. Set n_families = pmax(num_families, n_landowners)
# 3. Set n_landless = n_families - n_landowners
#------------------------------------------------------------------------------

cross_section %<>% # Assignment pipe!
  rowwise(municipality) %>%
  mutate(n_landowners = sum(c_across(starts_with('landown_')))) %>%
  ungroup() %>%
  rename(n_families_census = num_families) %>%
  mutate(n_families = pmax(n_families_census, n_landowners),
         n_landless = n_families - n_landowners,
         omega_n = n_landless / n_families) %>%
  select(-landless_families)


#---------------------------------------------------------------------------------------
# We observe a discretized land distribution: we know the total amount of land and the
# total number of landholders within each of the following "bins" measured in hectares
#---------------------------------------------------------------------------------------
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
#---------------------------------------------------------------------------------------
landowner_bins <- cross_section %>% # Number of families in each landholding "bin"
  select(n_landless, starts_with('landown_')) %>%
  relocate(n_landless) %>% # Make this the first column
  rename(`0` = n_landless,
         `(0,1)` = landown_less_1,
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

landtotal_bins <- cross_section %>% # Total land in each landholding "bin"
  mutate(tot_land_0 = 0) %>% # landless families are those with exactly 0 land
  select(starts_with('tot_land_')) %>%
  relocate(tot_land_0) %>% # Make this the first column
  rename(`0` = tot_land_0,
         `(0,1)` = tot_land_1,
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

# Construct a list of land distributions: one for each municipality
make_land_dist <- function(i) {
  total_land <- landtotal_bins[i,]
  n_families <- landowner_bins[i,]
  out <- data.frame(t(rbind(total_land = landtotal_bins[i,], n_families = landowner_bins[i,])))
  out$mean_land <- ifelse(out$n_families > 0, out$total_land / out$n_families, 0)
  out$frac_families <- out$n_families / sum(out$n_families)
  return(out)
}
land_distributions <- lapply(seq_len(nrow(cross_section)), make_land_dist)
names(land_distributions) <- cross_section$municipality
rm(make_land_dist, landowner_bins, landtotal_bins)

# Function that computes summary statistics of the land distribution for use in
# our LINEAR IV models. The non-linear IV models rely on the entire land distribution
# stored in land_distributions
get_land_statistics <- function(x) {
  z <- x[-1,] # summary stats for *landholders* so remove landless
  z$frac_families <- z$frac_families / sum(z$frac_families)
  c('omega' = x$frac_families[1],
    'H1' = sum(asinh(z$mean_land) * z$frac_families),
    'H2' = sum(asinh(z$mean_land)^2 * z$frac_families),
    n_families = sum(x$n_families))
}

land_statistics <- lapply(land_distributions, get_land_statistics)
land_statistics <- as.data.frame(do.call(rbind, land_statistics))
land_statistics$municipality <- as.numeric(rownames(land_statistics))
rownames(land_statistics) <- NULL
rm(get_land_statistics)

cross_section %<>% # Assignment pipe!
  left_join(y = land_statistics) %>%
  mutate(omegaC = 1 - omega) # It's convenient to have (1 - omega) as a separate variable

rm(land_statistics)

# Cumulative violence in the final year of the panel, i.e. total violence, can
# be viewed as a cross-sectional characteristic. Merge this with the other
# cross_section data

merge_me <- panel %>%
  filter(year == max(year)) %>%
  select(municipality, V_cum)

cross_section %<>%
  inner_join(merge_me) %>%
  rename(V_total = V_cum) # Better to call it V_total: total violence from 1996:2008

rm(merge_me)

#------------------------ Flag municipalities w/ land data and >=100 landholders

# The municipalities in cross_section are those for which we have land dist data
identical(sum(!is.na(cross_section$H1)), nrow(cross_section))

cross_section %<>%
  mutate(at_least_100_landowners = n_landowners >= 100)

munis_100_plus <- cross_section %>%
  filter(at_least_100_landowners) %>%
  pull(municipality)

panel %<>%
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
disagreement <- panel %>%
  select(municipality, year, starts_with('D_')) %>%
  group_by(municipality) %>%
  summarise(across(.cols = starts_with('D_'),
                   .fns = sum, na.rm = TRUE)) %>%
  rowwise() %>%
  filter(sum(c_across(starts_with('D_')) == 0) == 1) %>%
  ungroup()


AS_replace <- disagreement %>%
  filter(D_AS == 0) %>%
  pull(municipality)
panel[panel$municipality %in% AS_replace,]$D_AS <- NA

CODHES_replace <- disagreement %>%
  filter(D_CODHES == 0) %>%
  pull(municipality)
panel[panel$municipality %in% CODHES_replace,]$D_CODHES <- NA

RUV_replace <- disagreement %>%
  filter(D_RUV == 0) %>%
  pull(municipality)
panel[panel$municipality %in% RUV_replace,]$D_RUV <- NA

CEDE_replace <- disagreement %>%
  filter(D_CEDE == 0) %>%
  pull(municipality)
panel[panel$municipality %in% CEDE_replace,]$D_CEDE <- NA

rm(disagreement, AS_replace, CODHES_replace, RUV_replace, CEDE_replace)

#--------------------------------------- Create lagged violence flow
panel <- panel %>%
  group_by(municipality) %>%
  mutate(lag_V_flow = lag(V_flow, order_by = year)) %>%
  ungroup()

#--------------------------- Re-scale displacement measures

# What share of total RUV and JYP displacement was reported in 1996?
# Average these two shares
share_1996 <- panel %>%
  select(year, D_RUV, D_JYP) %>%
  group_by(year) %>%
  summarize(RUV_sum = sum(D_RUV, na.rm = TRUE),
            JYP_sum = sum(D_JYP, na.rm = TRUE)) %>%
  mutate(RUV_share = RUV_sum / sum(RUV_sum),
         JYP_share = JYP_sum / sum(JYP_sum)) %>%
  filter(year == 1996) %>%
  select(RUV_share, JYP_share) %>%
  rowMeans()

panel %<>%
  mutate(across(c(D_RUV, D_JYP), # Normalize total displacement to 6 million
         list(norm = ~ 6e6 * . / sum(., na.rm = TRUE)))) %>%
  mutate(across(c(D_AS, D_CEDE),
                list(norm = ~ (1 - share_1996) * 6e6 * . / sum(., na.rm = TRUE)))) %>%
  rowwise() %>%
  mutate(D_med = median(c_across(ends_with('_norm')), na.rm = FALSE),
         D_avg = mean(c_across(ends_with('_norm')), na.rm = FALSE),
         D_med_any = median(c_across(ends_with('_norm')), na.rm = TRUE),
         D_avg_any = mean(c_across(ends_with('_norm')), na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(D_avg_no_AS = (D_RUV + D_CEDE + D_JYP) / 3,
         D_avg_no_RUV = (D_AS + D_CEDE + D_JYP) / 3,
         D_avg_no_CEDE = (D_AS + D_RUV + D_JYP) / 3,
         D_avg_no_JYP = (D_AS + D_RUV + D_CEDE) / 3) %>%
  mutate(across(c(starts_with('D_med'), starts_with('D_avg')), ~ 6e6 * . / sum(., na.rm = TRUE)))

rm(share_1996)


# Covariates ------------------------------------------------------------------

cross_section %<>%

  rename(mines = mine_titles_90, oil = oil_prod_98) %>%

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
                                    hospitals_95)))) %>%

  select(-starts_with('tot_land'), -starts_with('landown_'), -gini,
         -convexity_L, -starts_with('owners_'), -ends_with('_95'),
         -starts_with('ratio_votes'), -n_families_census,
         -contains('_signal'), -contains('educ_'), -dum_2001_land)


# Spatial stuff ----------------------------------------------------------------

# Read in ARCGIS-generated list of municipalities intersecting road network
roads <- readr::read_csv('data-raw/real_roads.csv') %>%
  mutate(municipality = stringr::str_remove(ADM2_PCODE, 'CO'),
         municipality = as.numeric(municipality))


# Read and clean forest & elevation data --------------------------------------

forests <- readr::read_csv("data-raw/Forest_Master_50.csv") %>%
  mutate(share_forested = extent_2000_ha / area_ha,
         is_forested = share_forested > 0.5) %>%
  select(ADM2_PCODE, is_forested) %>%
  mutate(municipality = stringr::str_remove(ADM2_PCODE, 'CO'),
         municipality = as.numeric(municipality)) %>%
  select(-ADM2_PCODE)

# Municipality 70523 (Palmito, Sucre) appears twice in forests, once with
# is_forested TRUE and again with is_forested FALSE. There are no forests here!
forests %<>%
  filter(!((municipality == 70523) & (is_forested)))

elevation <- readr::read_csv("data-raw/elevationsFinal.csv") %>%
  select(-`...1`) %>% # first column is a row number starting from zero
  mutate(municipality = stringr::str_remove(ADM2_PCODE, 'CO'),
         municipality = as.numeric(municipality)) %>%
  select(-ADM2_PCODE) %>%
  arrange(municipality)


# Create dataframe of geographic information ----------------------------------
geography <- cross_section %>%
  select(municipality, ruggedness, slope, elevation) %>%
  left_join(elevation) %>%
  left_join(forests) %>%
  mutate(has_road = if_else(municipality %in% roads$municipality, 1, 0),
         is_forested = 1 * is_forested)

rm(elevation, forests, roads)


#-------------------------------------------------------------------------------
# Shapefile for Colombian municipalities
#    source:   <https://data.humdata.org/dataset/cod-ab-col>
#    filename: "Col Administrative Divisions Shapefiles.zip"
#-------------------------------------------------------------------------------

# Read in shapefile.
library(sf)
library(spdep)
library(geosphere)

CO_shape <- st_read('data-raw/CO-shapefiles/col_admbnda_adm2_mgn_20200416.shp')
centroid_coords <- CO_shape %>%
  st_centroid() %>% # compute lat and lon of centroids
  st_coordinates() # extract lat and lon of centroids as *matrix*

neighbors <- poly2nb(CO_shape) # construct list of neighbors

# Note: islands get a 0 for their neighbors
# CO_shape$ADM2_PCODE[705]
# neighbors[[705]]

get_neighbor_distances <- function(i) {
  if(all(neighbors[[i]] == 0)) {
    return(NA) # No neighbors? No distances!
  } else {
    own_coords <- centroid_coords[i, ]
    neighbor_coords <- centroid_coords[neighbors[[i]], ]
    distGeo(own_coords, neighbor_coords) / 1000 # convert meters to km
  }
}
neighbor_distances <- lapply(1:length(neighbors), get_neighbor_distances)

get_pairwise_distances <- function(i) {
# Only return neighbors whose index is greater than i: eliminate "repeat" pairs
  i_neighbors <- neighbors[[i]]
  i_neighbor_distances <- neighbor_distances[[i]]
  keep <- i < i_neighbors
  if(any(keep)) {
    data.frame(i, j = i_neighbors[keep], dist = i_neighbor_distances[keep])
  } else {
    # Return NULL if no neighbors (entry in neighbors is zero) or all repeats
    NULL
  }
}


pairwise_distances <- lapply(1:length(neighbors), get_pairwise_distances)
pairwise_distances <- do.call(rbind, pairwise_distances)
pairwise_distances <- tibble(pairwise_distances)

rm(get_neighbor_distances, get_pairwise_distances)


# Impute missing values --------------------------------------------------------

# There are 41 missing values for the variables elevation_difference and
# is_forested from geography. We impute these with an average of the same
# variables for *neighboring* municipalities. If a variable is missing for all
# neighboring municipalities, we fill in the overall mean of the variable for
# geography.


# Create list of neighboring municipalities with municipality *codes* rather
# than row numbers from CO_shape
muni_codes <- as.numeric(stringr::str_remove(CO_shape$ADM2_PCODE, 'CO'))
neighbors_codes <- lapply(neighbors, function(x) muni_codes[x])
names(neighbors_codes) <- muni_codes

# geography only has data for 1049 municipalities while CO_shape has data for
# 1122. Add the missing municipalities to geography
geography %<>%
  full_join(tibble(municipality = muni_codes))


impute_missing_var <- function(var_name) {
  muni_missing_var <- geography[is.na(geography[, var_name]), 'municipality', drop = TRUE]
  municipality_indices <- which(muni_codes %in% muni_missing_var)
  neighbor_list <- neighbors_codes[municipality_indices]
  f <- function(x) {
    neighbor_geography <- pull(geography[geography$municipality %in% x, var_name])
    mean(neighbor_geography, na.rm = TRUE)
  }
  out <- lapply(neighbor_list, f)
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
imputed_geography <- lapply(names(geography)[-1], impute_missing_var)

imputed_geography %<>%
  reduce(full_join, by = 'municipality')


# TODO : last step of the imputation, following the logic of the comment, but
# for all the imputed variables, not just impute_forest
geography %>%
  full_join(imputed_geography)

# Implement the following logic *repeatedly* for all variables that start with
# 'impute_', is that mutate_at()?
#geography <- geography %>%
#  left_join(impute_forest) %>%
#  mutate(is_forested = if_else(!is.na(is_forested),
#         is_forested, impute_is_forested))

rm(muni_missing_forest, muni_missing_elev, get_neighbor_geography_avg,
   imputed_geography, muni_codes, neighbors_codes)


# Construct dataset of all unique *pairs* of neighboring municipalities --------

# Create geography_i by renaming the columns of geography to have an '_i'
# at the end. Do the same for geography_j. Then merge with pairwise distances.
geography_i <- geography %>%
  rename_with(~ stringr::str_c(., '_i'))

geography_j <- geography %>%
  rename_with(~ stringr::str_c(., '_j'))

pairwise_distances %<>% # Assignment pipe!
  mutate(municipality_i = stringr::str_remove(CO_shape$ADM2_PCODE[i], 'CO'),
         municipality_j = stringr::str_remove(CO_shape$ADM2_PCODE[j], 'CO'),
         municipality_i = as.numeric(municipality_i),
         municipality_j = as.numeric(municipality_j)) %>%
  select(-i, -j) %>%
  left_join(geography_i) %>%
  left_join(geography_j)

rm(geography_i, geography_j, geography)

link_stats <- pairwise_distances %>%
  mutate(forest = is_forested_i + is_forested_j,
         elevation = (elevation_difference_i + elevation_difference_j) / 2,
         ruggedness = (ruggedness_i + ruggedness_j) / 2,
         slope = (slope_i + slope_j) / 2,
         elevation_diff = abs(elevation_i - elevation_j)) %>%
  select(municipality_i, municipality_j, dist, forest, elevation, slope,
         ruggedness, elevation_diff)

pca <- prcomp(~ forest + elevation + ruggedness + slope + elevation_diff,
              data = link_stats, scale = TRUE, na.action = na.omit)

summary(pca)
print(pca)



# TO DO NEXT! -----------------------------------------------------------------
# - Go back and impute the missing geographic data for the 73 or so municipalities
#     from the shapefile that aren't in our cross-section dataset.
# - Extract and store the first PC
# -----------------------------------------------------------------------------

#---------------------------------------------------
#------------------ UPDATE BELOW!!!!
#---------------------------------------------------

# Take principal components from prcomp()
pca <- prcomp(df,scale = TRUE)
pca_noroad <- as.data.frame(pca$x[,1:2])
pca_noroad <- mutate(pca_noroad,index = as.integer(rownames(pca_noroad)))
pca_noroad <- merge(pca_noroad,for_merge_df,by="index")
pca_noroad <- arrange(pca_noroad,index)


# Take principal components from prcomp()
pca <- prcomp(df,scale = TRUE)
pca_road <- as.data.frame(pca$x[,1:2])
pca_road <- mutate(pca_road, index = as.integer(rownames(pca_road)))
pca_road <- merge(pca_road, for_merge_df, by="index")
pca_road <- arrange(pca_road,index)

#-------------------------------------------------------------------------------
# Create abandoned land dataset.

# Note: There was an intermediate hand-cleaning step that was
# omitted between merging raw digitized with AttributeTableFinal.csv and input
# here to remove duplicates; the abandoned_land_handcleaned.csv file is the
# result of that cleaned merge.
#-------------------------------------------------------------------------------

#mrg <- read.csv("data-raw/abandoned_land_handcleaned.csv")

#abandoned_land <- subset(mrg,select = c(displaced,hect_abandoned_paramilitary,
#hect_abandoned_other_armed,total_hect_abandoned,ADM2_ES,ADM2_PCODE,adm2Nm,lat,lon))

#usethis::use_data(covariates, overwrite = TRUE)
#usethis::use_data(cross_section, overwrite = TRUE)
usethis::use_data(panel, overwrite = TRUE)
#usethis::use_data(displacement, overwrite = TRUE)
#usethis::use_data(Vcum, overwrite = TRUE)
#usethis::use_data(Vcum_pop, overwrite = TRUE)
usethis::use_data(land_distributions, overwrite = TRUE)
usethis::use_data(adjacent_municipalities, overwrite = TRUE)
#usethis::use_data(municipalities_and_regions, overwrite = TRUE)
usethis::use_data(geographic_covariates, overwrite = TRUE)
usethis::use_data(munigraph, overwrite = TRUE)
usethis::use_data(violence_data, overwrite = TRUE)
usethis::use_data(pca_noroad, overwrite = TRUE)
usethis::use_data(pca_road, overwrite = TRUE)
#usethis::use_data(abandoned_land,overwrite = TRUE)

# clean up
rm(list = ls())


