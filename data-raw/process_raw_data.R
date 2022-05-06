library(dplyr)
library(magrittr)
library(tidyr)
library(parallel)

#------------------------------- Load raw data
# The cross section dataset contains information for only those municipalities
# with land distribution data
cross_section_raw <- haven::read_dta("data-raw/cross-section-raw.dta")

# The panel dataset contains information for municipalities that we cannot
# use in our models because we lack land distribution data for them
panel_raw <- haven::read_dta("data-raw/panel-raw.dta")

survey_raw <- haven::read_dta("data-raw/survey-raw.dta")

municipalities_and_regions_raw <- readr::read_csv("data-raw/DANE_municipality_codes_and_regions.csv")

adjacency_matrix_raw <- readxl::read_excel("data-raw/Adjacency_Matrix.xlsx")
# Note: the adjacency matrix predated Parker. Adjacency_Matrix.xlsx doesn't
# perfectly match the adjacency matrix generated using the shapefile.
# (HOW DO THEY DIFFER AND WHY? Is it just the islands?)


#------------------------------- Merge department/province names into cross-section
cross_section <- municipalities_and_regions_raw %>%
  rename(department = number_dept, municipality = Muni_code) %>%
  mutate(municipality = as.numeric(municipality)) %>%
  left_join(cross_section_raw, .)

rm(municipalities_and_regions_raw, cross_section_raw)

#------------------------------- Remove two islands: 88001 and 88564
CO_islands <- c('88001', '88564')

survey_raw %>%
  filter(municipality_origin %in% CO_islands) #None in the survey

cross_section %<>% # Assignment pipe!
  filter(!(municipality %in% CO_islands))

panel <- panel_raw %>%
  filter(!(municipality %in% CO_islands))

adjacency_matrix <- adjacency_matrix_raw %>%
  filter(!(codes %in% CO_islands)) # Note that we don't have to remove
                                   # columns 'border_with_88001' etc. since
                                   # these two municipalities are not in the
                                   # adjacency matrix: 1117 versus 1119 obs.

rm(panel_raw, adjacency_matrix_raw, CO_islands)


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

# There are some variables in the panel dataset that are really cross-section
# variables: e.g. 1993 population, cumulative violence in the final year of the
# panel, i.e. total violence. Merge these into the cross-section dataset

merge_me <- panel %>%
  filter(year == max(year)) %>%
  select(municipality, popn1993, V_cum)

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






#---------------------------------------------------
#------------------ PARKER STUFF BELOW!!!!
#---------------------------------------------------

## Read in violence data.
violence_data <- read.csv("data-raw/COL_muni_flows.csv")



## Create full merged cross_section data.

# Read in and clean forest data.
forests = read.csv("data-raw/Forest_Master_50.csv")
forests = dplyr::mutate(forests, share_forested = extent_2000_ha / area_ha)
forests = dplyr::mutate(forests, is_forested = (share_forested>0.5))
forest_merge = dplyr::select(forests,c("ADM2_PCODE","is_forested"))

# Read in elevation data.
elevation = read.csv("data-raw/elevationsFinal.csv")

# Define helper function to create merge column.
addCO <- function(muni){
  ifelse(nchar(paste(muni)) == 4,paste("CO0",muni,sep = ""),paste("CO",muni,sep = ""))
}

# Helper function to convert ADM2_PCODE to municipality code.
removeCO <- function(muni){ return(paste(as.numeric(substring(muni,3,7))))
}

# Merge forest and elevation into existing geographic covariate dataframe.
full_geography = dplyr::select(cross_section,municipality,lat_mean,lon_mean,ruggedness,slope,elevation,V_cum)
full_geography <- dplyr::mutate(full_geography,ADM2_PCODE = purrr::map_chr(municipality,addCO))
full_geography <- merge(full_geography,elevation,by = "ADM2_PCODE")
full_geography <- merge(full_geography,forest_merge,by="ADM2_PCODE")

# Read in web-scraped municipality dataset with latitudes and longitudes.
AttributeTableFinal <- read.csv("data-raw/AttributeTableFinal.csv")
colnames(AttributeTableFinal)[1] <- "ADM2_ES"
AttributeTableFinal <- dplyr::mutate(AttributeTableFinal,latnum = as.numeric(lat), lonnum = as.numeric(lon))

# Read in ARCGIS-generated list of municipalities intersecting road network and create
# dummy indicating whether each municipality contains a road.
roads = read.csv("data-raw/real_roads.csv")
road_code_list <- roads$ADM2_PCODE
AttributeTableFinal <- dplyr::mutate(AttributeTableFinal, has_road = (ADM2_PCODE %in% road_code_list))

colnames(AttributeTableFinal)[1] <- "ADM2_ES"
full_municipalities <- subset(AttributeTableFinal,select = c(ADM2_ES,ADM2_PCODE,
admn2Nm,lat,lon))

cross_section_merged <- merge(AttributeTableFinal,full_geography,by = "ADM2_PCODE")
cross_section_merged <- dplyr::select(cross_section_merged,c("ADM2_PCODE","admn2Nm","latnum","lonnum","has_road","municipality","ruggedness","slope","elevation","elevation_difference","is_forested"))
geographic_covariates <- cross_section_merged


## Create spatial map object and graph object.(The graph and adjacency matrix produced is identical to those produced by
## ArcGis using this shapefile, but this method is easier to document).

# Read in shapefile.
muni_pol <- st_read("data-raw/col muni polygons/col_admbnda_adm2_mgn_20200416.shp")

# Remove island municipalities.
muni_pol <- subset(muni_pol,ADM2_PCODE != "CO88001")
muni_pol <- subset(muni_pol,ADM2_PCODE!="CO88564")
muni_pol <- mutate(muni_pol,municipality = purrr::map_chr(ADM2_PCODE,removeCO))

# Using spdep, transform spatial polgygons into adjacency matrix.
polygon <- left_join(muni_pol,AttributeTableFinal,by=c("ADM2_PCODE"))
nb <- spdep::poly2nb(polygon)
listw <- spdep::nb2listw(nb,zero.policy = TRUE)
mat2 <- spdep::listw2mat(listw)

# Create graph object from generated adjacency matrix.
munigraph <- igraph::graph_from_adjacency_matrix(mat2,mode ="undirected",weighted = TRUE)


## Generate principal components. This was originally done in a separate file.

# To allow us to take the maximum over crow-flies distances for all pairs of municipalities, we save the crow-flies
# distance for each adjacent pair.
adjacent_pairs <- list()
d1_list <- list()

for(i in 1:nrow(cross_section_merged)){
  for(j in 1:nrow(cross_section_merged)){
    if(igraph::are.connected(munigraph,i,j) & i != j){

      # Add in adjacent pairs.
      adjacent_pairs <- append(adjacent_pairs,c(i,j))

      # Calculate d1.
      d1 <- geosphere::distGeo(c(AttributeTableFinal$latnum[i],AttributeTableFinal$lonnum[i]),c(AttributeTableFinal$latnum[j],AttributeTableFinal$lonnum[j]))/1000
      d1_list <- append(d1_list,d1)
    }
  }
}

# Initialize row counter, and matrices.
count <- 1
pca_matrix <- matrix(0,nrow = 10000,ncol = 6)
for_merge <- matrix(0,nrow = 10000,ncol = 5)
for(i in 1:nrow(cross_section_merged)){
  for(j in 1:nrow(cross_section_merged)){
    if(igraph::are.connected(munigraph,i,j) & i != j){

      # Calculate crow-flies distance.
      d1 <- geosphere::distGeo(c(AttributeTableFinal$latnum[i],AttributeTableFinal$lonnum[i]),c(AttributeTableFinal$latnum[j],AttributeTableFinal$lonnum[j]))

      # Find the maximum for each covariate.

      max_elevation <- max(cross_section_merged$elevation)
      max_ruggedness <- max(cross_section_merged$ruggedness)
      max_slope <- max(cross_section_merged$slope)
      max_forest <- 1
      max_withinmuni_elevation_difference <- max(cross_section_merged$elevation_difference)

      # Get the max of crow-flies distance over adjacent municipalities from our previously-generated list.
      max_d1 <- max(unlist(d1_list))


      # Calculate the average of each covariate across the two municipalities.

      elevation_difference <- max((cross_section_merged$elevation[i]-cross_section_merged$elevation[j]),0)
      average_ruggedness <- 0.5*cross_section_merged$ruggedness[i]+0.5*cross_section_merged$ruggedness[j]
      average_slope <- (0.5*cross_section_merged$slope[i]+0.5*cross_section_merged$slope[j])
      average_forest <- 0.5*as.integer(cross_section_merged$is_forested[i])+0.5*as.integer(cross_section_merged$is_forested[j])
      average_withinmuni_elevation_difference <- 0.5*cross_section_merged$elevation_difference[i]+0.5*cross_section_merged$elevation_difference[j]

      # Get the max of crow-flies distance over adjacent municipalities from our previously-generated list.
      max_d1 <- max(unlist(d1_list))

      # Standardize averages by dividing by maximum.

      pca_matrix[count,1] <- elevation_difference / max_elevation
      pca_matrix[count,2] <- average_ruggedness / max_ruggedness
      pca_matrix[count,3] <- average_slope / max_slope
      pca_matrix[count,4] <- average_forest / max_forest
      pca_matrix[count,5] <- average_withinmuni_elevation_difference / max_withinmuni_elevation_difference
      pca_matrix[count,6] <- ifelse(is.na(d1),0,d1) / max_d1

      # To identify each municipality, save municipality name and administrative code in a separate dataframe to merge in.
      for_merge[count,1] <- AttributeTableFinal$ADM2_PCODE[i]
      for_merge[count,2] <- AttributeTableFinal$ADM2_PCODE[j]
      for_merge[count,3] <- AttributeTableFinal$ADM2_ES[i]
      for_merge[count,4] <- AttributeTableFinal$ADM2_ES[j]
      for_merge[count,5] <- count

      # Update row counter.
      count <- (count+1)
  }
 }
}

# Create our two dataframes and set column names.
for_merge_df <- as.data.frame(for_merge)
colnames(for_merge_df) <- c("muni_i","muni_j","muni_i_name","muni_j_name","index")

df <- as.data.frame(pca_matrix)
colnames(df) <- c("elev_difference","average_ruggedness","average_slope","average_forest","max_withinmuni_elevation_difference","d1")

# Remove unused rows (all 0;s) from dataframe and rescale.
df <- df[rowSums(df)>0,]
scaled <- scale(df)


# Take principal components from prcomp()
pca <- prcomp(df,scale = TRUE)
pca_noroad <- as.data.frame(pca$x[,1:2])
pca_noroad <- mutate(pca_noroad,index = as.integer(rownames(pca_noroad)))
pca_noroad <- merge(pca_noroad,for_merge_df,by="index")
pca_noroad <- arrange(pca_noroad,index)



## Same PCA as above, but run with roads.

# Initialize row counter, matrices, and list of adjacent pairs.
count <- 1
pca_matrix <- matrix(0,nrow = 10000,ncol = 7)
for_merge <- matrix(0,nrow = 10000,ncol = 5)

for(i in 1:nrow(cross_section_merged)){
  for(j in 1:nrow(cross_section_merged)){
    if(igraph::are.connected(munigraph,i,j) & i != j){

    # Calculate crow-flies distance.
    d1 <- geosphere::distGeo(c(AttributeTableFinal$latnum[i],AttributeTableFinal$lonnum[i]),c(AttributeTableFinal$latnum[j],AttributeTableFinal$lonnum[j]))

    # Find the maximum for each covariate.

    max_elevation <- max(cross_section_merged$elevation)
    max_ruggedness <- max(cross_section_merged$ruggedness)
    max_slope <- max(cross_section_merged$slope)
    max_forest <- 1
    max_withinmuni_elevation_difference <- max(cross_section_merged$elevation_difference)
    max_road <- 1

    # Calculate the average of each covariate across the two municipalities.

    elevation_difference <- max((cross_section_merged$elevation[i]-cross_section_merged$elevation[j]),0)
    average_ruggedness <- 0.5*cross_section_merged$ruggedness[i]+0.5*cross_section_merged$ruggedness[j]
    average_slope <- (0.5*cross_section_merged$slope[i]+0.5*cross_section_merged$slope[j])
    average_forest <- 0.5*as.integer(cross_section_merged$is_forested[i])+0.5*as.integer(cross_section_merged$is_forested[j])
    average_withinmuni_elevation_difference <- 0.5*cross_section_merged$elevation_difference[i]+0.5*cross_section_merged$elevation_difference[j]
    average_road <- 0.5*cross_section_merged$has_road[i]+0.5*cross_section_merged$has_road[j]


    # Standardize averages by dividing by maximum.

    pca_matrix[count,1] <- elevation_difference / max_elevation
    pca_matrix[count,2] <- average_ruggedness / max_ruggedness
    pca_matrix[count,3] <- average_slope / max_slope
    pca_matrix[count,4] <- average_forest / max_forest
    pca_matrix[count,5] <- average_withinmuni_elevation_difference / max_withinmuni_elevation_difference
    pca_matrix[count,6] <- ifelse(is.na(d1),0,d1) / max_d1

    # Take the negative of road.
    pca_matrix[count,7] <- -1*average_road

    # To identify each municipality, save municipality name and administrative code in a separate dataframe to merge in.
    for_merge[count,1] <- AttributeTableFinal$ADM2_PCODE[i]
    for_merge[count,2] <- AttributeTableFinal$ADM2_PCODE[j]
    for_merge[count,3] <- AttributeTableFinal$ADM2_ES[i]
    for_merge[count,4] <- AttributeTableFinal$ADM2_ES[j]
    for_merge[count,5] <- count

    # Update row counter.
    count <- (count+1)
  }
 }
}


# Create our two dataframes and set column names.
for_merge_df <- as.data.frame(for_merge)
colnames(for_merge_df) <- c("muni_i","muni_j","muni_i_name","muni_j_name","index")

df <- as.data.frame(pca_matrix)
colnames(df) <- c("elev_difference","average_ruggedness","average_slope","average_forest","max_withinmuni_elevation_difference","d1")

# Remove unused rows (all 0;s) from dataframe and rescale.
df <- df[rowSums(df)>0,]
scaled <- scale(df)


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















