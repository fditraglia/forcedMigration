library(dplyr) # not used elsewhere in the package: only for data prep
library(readr) # likewise

# Load raw data in STATA format
cross_section_raw <- haven::read_dta("data-raw/Cross_Section.dta")
panel_raw <- haven::read_dta("data-raw/Panel_D_V.dta")
survey_raw <- haven::read_dta("data-raw/Survey.dta")

# The adjacency matrix for Colombian municipalities is stored as an excel file
adjacency_matrix <- readxl::read_excel("data-raw/Adjacency_Matrix.xlsx")

# Check that the codes in column 1 match the names in columns 2, 3, ...
municipality_codes <- adjacency_matrix$codes
adjacency_matrix <- adjacency_matrix[,-1]
row_indices <- as.integer(substr(colnames(adjacency_matrix), start = 13, stop = 100L))
sum(abs(municipality_codes - row_indices))

# Convert adjacency matrix into a list of neighbors for each municipality
get_neighbors <- function(row_index) {
  adj_row <- adjacency_matrix[row_index,]
  adj_row <- if_else(adj_row == 1, TRUE, FALSE)
  municipality_codes[adj_row]
}
adjacent_municipalities <- lapply(1:nrow(adjacency_matrix), get_neighbors)
names(adjacent_municipalities) <- municipality_codes
rm(municipality_codes, adjacency_matrix, get_neighbors, row_indices)

# Load csv file with municipality and region codes, merge with cross_section_raw
#-------------------------------------------------------------------------------
# NOTE: the adjacency matrix above has only 1117 municipalities, whereas this
# csv file has 1119. Why? We think this is because of two islands. Check this.
#-------------------------------------------------------------------------------
municipalities_and_regions <- read_csv("data-raw/DANE_municipality_codes_and_regions.csv")
names(municipalities_and_regions) <- c('department', 'region', 'municipality')
municipalities_and_regions$municipality <- as.numeric(municipalities_and_regions$municipality)
cross_section_raw <- left_join(cross_section_raw, municipalities_and_regions)
rm(municipalities_and_regions)

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


#------------------------------------------------------------------------------
# The column landless_families in from Cross_section.dta was constructed by
# Camilo's original RA. We can't figure out what he did, but the results don't
# make any sense. The following is our procedure for constructing a measure of
# the number of landless families in a municipality from the underlying census
# data. We assume, absent evidence to the contrary, that there is at most one
# landowner per family.
#------------------------------------------------------------------------------
# 1. Calculate n_landowners by summing counts in each column "landown_*"
# 2. Set n_families = pmax(num_families, n_landowners)
# 3. Set n_landless = n_families - n_landowners
#------------------------------------------------------------------------------

cross_section <- cross_section_raw %>%
  rowwise(municipality) %>%
  mutate(n_landowners = sum(c_across(starts_with('landown_')))) %>%
  ungroup() %>%
  rename(n_families_census = num_families) %>%
  mutate(n_families = pmax(n_families_census, n_landowners),
         n_landless = n_families - n_landowners,
         omega_n = n_landless / n_families) %>%
  select(-landless_families)

rm(cross_section_raw)

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


# Select the panel variables we'll use and give them simpler names
panel <- panel_raw %>%
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

rm(panel_raw)

# There are some variables in the panel dataset that are really cross-section
# variables: e.g. 1993 population, cumulative violence in the final year of the
# panel. Merge these into the cross-section dataset

merge_me <- panel %>%
  filter(year == max(year)) %>%
  select(municipality, popn1993, V_cum)

cross_section <- inner_join(cross_section, merge_me, by = 'municipality')
rm(merge_me)


#--------------------------------------------------------
# Keep only municipalities with at least 100 landowners
#--------------------------------------------------------
keep_municipalities <- cross_section %>%
  filter(n_landowners >= 100) %>%
  pull(municipality)
cross_section <- cross_section %>%
  filter(municipality %in% keep_municipalities)
panel <- panel %>%
  filter(municipality %in% keep_municipalities)

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

#-------------------------------------------------------------------------------
# Create lagged violence flow
#-------------------------------------------------------------------------------
panel <- panel %>%
  group_by(municipality) %>%
  mutate(lag_V_flow = lag(V_flow, order_by = year)) %>%
  ungroup()


# Arrange displacement measures into 3d array: municipality-year-measure
raw_displacement <- panel %>%
  select(starts_with('D_'))
g <- function(i) {
  temp <- data.frame(municipality = panel$municipality, year = panel$year,
                     D = raw_displacement[,i])
  temp <- reshape(temp, direction = 'wide', timevar = 'year',
                  idvar = 'municipality')
  colnames(temp)[-1] <- as.character(1996:2012)
  rownames(temp) <- temp$municipality
  temp <- as.matrix(temp)
  temp <- temp[,-1]
  return(temp)
}
displacement <- sapply(1:ncol(raw_displacement), g, simplify = 'array')
dimnames(displacement)[[3]] <- unlist(lapply(strsplit(names(raw_displacement), '_'), function(x) x[2]))


Vcum <- reshape(as.data.frame(panel[,c('municipality', 'year', 'V_cum')]),
                direction = 'wide', idvar = 'municipality',
                timevar = 'year')
colnames(Vcum)[-1] <- as.character(1996:2012)
rownames(Vcum) <- Vcum$municipality
Vcum <- Vcum[,-1]
Vcum_pop <- Vcum / cross_section$popn1993

# Covariates (cross-section data)
bureaucracy <- with(cross_section, local_bureaucracy_95 + state_bureaucracy_95)
names(cross_section[c(11:22, 25:31)])
offices <- rowSums(cross_section[,c(11:22, 25:31)])
elec_comp <- with(cross_section, 0.5 * ratio_votes_dif_alc +
                    0.5 * ratio_votes_dif_asa )
covariates <- cross_section[, c('municipality',
                                    'land_return',
                                    'rainfall',
                                    'ruggedness',
                                    'coffee', # dummy variable
                                    'predicted_signal_strength_1',
                                    'elevation',
                                    'dist_cap',
                                    'drug_routes',
                                    'mine_titles_90',
                                    'oil_prod_98',
                                    'guerrilla',
                                    'grazing',
                                    'grasses')]

names(covariates)[which(names(covariates) %in% c('predicted_signal_strength_1',
                               'mine_titles_90',
                               'oil_prod_98'))] <- c('radio', 'mines', 'oil')

covariates$offices <- offices
covariates$bureaucracy <- bureaucracy
covariates$elec_comp <- elec_comp
rm(offices, bureaucracy, elec_comp)

# Drop municipalities with fewer than 100 people
covariates <- subset(covariates, municipality %in% keep_municipalities)


# Trasform some variables to dummies:
drug_routes <- 1 * (covariates$drug_routes > 1) # at least 1km of drug routes?
mines <- 1 * (covariates$mines > 0) # any mines?
oil <- 1 * (covariates$oil > 0) # any oil production?

# Overwrite original, untransformed variables
covariates$drug_routes <- drug_routes
covariates$mines <- mines
covariates$oil <- oil
rm(drug_routes, mines, oil)

# Transform the non-dummy variables to Uniform(-0.5, 0.5)
transform_me <- setdiff(names(covariates), c('municipality',
                                             'coffee', # already a dummy
                                             'drug_routes',
                                             'mines',
                                             'oil'))
mytransform <- function(x) {
  stopifnot(all(!is.na(x)))
  (rank(x, ties.method = 'random') / length(x)) - 0.5
}
set.seed(1234) # for random tie-breaking in ranks within mytransform
covariates[,transform_me] <- apply(covariates[,transform_me], 2, mytransform)

# clean up
rm(mytransform, transform_me)

#usethis::use_data(covariates, overwrite = TRUE)
usethis::use_data(cross_section, overwrite = TRUE)
usethis::use_data(panel, overwrite = TRUE)
#usethis::use_data(displacement, overwrite = TRUE)
#usethis::use_data(Vcum, overwrite = TRUE)
#usethis::use_data(Vcum_pop, overwrite = TRUE)
usethis::use_data(land_distributions, overwrite = TRUE)
usethis::use_data(adjacent_municipalities, overwrite = TRUE)
#usethis::use_data(municipalities_and_regions, overwrite = TRUE)

# clean up
rm(list = ls())















