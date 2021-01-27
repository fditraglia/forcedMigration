library(dplyr) # not used elsewhere in the package: only for data prep

# Load raw data in STATA format
cross_section_raw <- haven::read_dta("data-raw/Cross_Section.dta")
panel_raw <- haven::read_dta("data-raw/Panel_D_V.dta")
survey_raw <- haven::read_dta("data-raw/Survey.dta")

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

landowner_bins <- cross_section_raw %>% # Number of landowners in each bin
  mutate(landown_0 = landless_families) %>% # landless families are those with exactly 0 land
  select(starts_with('landown_')) %>%
  relocate(landown_0) %>% # Make this the first column
  rename(`0` = landown_0,
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


#------------------------------------------------------------------------------
# We don't know what the RA did to construct the column landless_families in
# the cross-section dta file. Here's how we'll create our version
#------------------------------------------------------------------------------
# 1. Sum number of landholders -> n_landholding_families
# 2. If num_families >= n_landholding_families, set the difference as n_landless_families
# 3. If num_families < n_landholding_families:
#     a. Set num_families = n_landholding_families
#     b. Set n_landless_families = 0
#------------------------------------------------------------------------------

# How many landowners in each municipality?
n_landowners <- rowSums(landowner_bins)

landtotal_bins <- cross_section_raw %>% # Total land in each bin
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
land_distributions <- lapply(seq_len(nrow(cross_section_raw)), make_land_dist)
names(land_distributions) <- cross_section_raw$municipality
rm(make_land_dist)

cross_section <- cross_section_raw %>%
  mutate(omega_n = landless_families / num_families) %>%
  select(municipality, n_families, omega_n)

# Extract 1993 population from panel_raw and merge with cross_section
popn1993 <- panel_raw %>%
  filter(year == 1996) %>%
  select(municipality, Total_Poblacion1993) %>%
  rename(popn1993 = Total_Poblacion1993)


cross_section <- merge(cross_section, popn1993)



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

#-------------------------------------------------
# Update from here down! 2020-12-15
#-------------------------------------------------

cum_violence <- subset(panel, year == max(panel$year))[,c("municipality", "V_cum")]
cross_section <- tibble::as_tibble(merge(cross_section, cum_violence))
names(cross_section[, 'V_cum']) <- 'VTotal'

# Keep only municipalities with at least 100 landowners
keep_me <- n_landholders >= 100
keep_municipalities <- cross_section[keep_me,]$municipality
cross_section <- subset(cross_section, municipality %in% keep_municipalities)
panel <- subset(panel, municipality %in% keep_municipalities)
rm(keep_me)

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
# Create lagged violence flow
#-------------------------------------------------------------------------------
panel <- panel %>%
  group_by(municipality) %>%
  mutate(lag_V_flow = lag(V_flow, order_by = year)) %>%
  ungroup()


# Arrange displacement measures into 3d array: municipality-year-measure
rawZ <- as.data.frame(panel[,c('D_AS', 'D_CODHES', 'D_RUV', 'D_CEDE', 'D_JYP')])
g <- function(i) {
  temp <- data.frame(municipality = panel$municipality, year = panel$year,
                     D = rawZ[,i])
  temp <- reshape(temp, direction = 'wide', timevar = 'year',
                  idvar = 'municipality')
  colnames(temp)[-1] <- as.character(1996:2012)
  rownames(temp) <- temp$municipality
  temp <- as.matrix(temp)
  temp <- temp[,-1]
  return(temp)
}
obsZ <- sapply(1:ncol(rawZ), g, simplify = 'array')
dimnames(obsZ)[[3]] <- unlist(lapply(strsplit(names(rawZ), '_'), function(x) x[2]))


# Create matrices of land parameters and cumulative violence
land_parameters <- cross_section[, c('p', 'q', 'H_bar', 'omega_n')]

Vcum <- reshape(as.data.frame(panel[,c('municipality', 'year', 'V_cum')]),
                direction = 'wide', idvar = 'municipality',
                timevar = 'year')
colnames(Vcum)[-1] <- as.character(1996:2012)
rownames(Vcum) <- Vcum$municipality
Vcum <- Vcum[,-1]
Vcum_pop <- Vcum / cross_section$popn1993

# Covariates (cross-section data)
bureaucracy <- with(cross_section_raw, local_bureaucracy_95 + state_bureaucracy_95)
names(cross_section_raw[c(11:22, 25:31)])
offices <- rowSums(cross_section_raw[,c(11:22, 25:31)])
elec_comp <- with(cross_section_raw, 0.5 * ratio_votes_dif_alc +
                    0.5 * ratio_votes_dif_asa )
covariates <- cross_section_raw[, c('municipality',
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

usethis::use_data(covariates, overwrite = TRUE)
usethis::use_data(cross_section, overwrite = TRUE)
usethis::use_data(panel, overwrite = TRUE)
usethis::use_data(obsZ, overwrite = TRUE)
usethis::use_data(land_parameters, overwrite = TRUE)
usethis::use_data(Vcum, overwrite = TRUE)
usethis::use_data(Vcum_pop, overwrite = TRUE)

# clean up
rm(list = ls())















