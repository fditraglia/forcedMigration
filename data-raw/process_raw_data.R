cross_section_raw <- haven::read_dta("data-raw/Cross_Section.dta")
panel_raw <- haven::read_dta("data-raw/Panel_D_V.dta")
survey_raw <- haven::read_dta("data-raw/Survey.dta")

# Calculate total land in each municipality
landtotal_cols <- sapply(names(cross_section_raw), substr, 1, 9) == 'tot_land_'
landtotal_bins <- cross_section_raw[, landtotal_cols]
total_land <- rowSums(landtotal_bins)
rm(landtotal_bins, landtotal_cols)

# Construct land distribution from raw data on landowner counts binned by amount
# of land
landowner_cols <- sapply(names(cross_section_raw), substr, 1, 8) == 'landown_'
landowner_bins <- cross_section_raw[, landowner_cols]
names(landowner_bins) <- c('[0,1)', '[1,3)', '[3,5)', '[5,10)', '[10,15)',
                         '[15,20)', '[20,50)', '[50,100)', '[100,200)',
                         '[200,500)', '[500,1000)', '[1000,2000)', '>=2000')
rm(landowner_cols)

# Calculate total number of landholders and mean landholdings
n_landholders <- rowSums(landowner_bins)
mean_land_owned <- total_land / n_landholders
rm(total_land)

h <- function(x){
  x <- x / sum(x)
  n <- length(x)
  cumsum(x)[1:(n-1)]
}

land_cdfs <- t(apply(landowner_bins, 1, h))
colnames(land_cdfs) <- paste0('x=', c(1, 3, 5, 10, 15, 20, 50, 100, 200, 500,
                                     1000, 2000))

land_cdf_list <- lapply(seq_len(nrow(land_cdfs)), function(i) land_cdfs[i,])
highest_occupied <- apply(landowner_bins, 1, function(x) max(which(x > 0)))
land_cdf_list <- lapply(seq_along(land_cdf_list),
                        function(i) land_cdf_list[[i]][1:(highest_occupied[i] - 1)])
rm(h, land_cdfs, landowner_bins, highest_occupied)


fit_3_param_beta <- function(Fhat, mu){
  x_full <- c(1, 3, 5, 10, 15, 20, 50, 100, 200, 500, 1000, 2000)
  x <- x_full[1:length(Fhat)]
  f <- function(params){
    p <- params[1]
    q <- params[2]
    b <- mu * (p + q) / p
    Fmodel <- actuar::pgenbeta(x, shape1 = p, shape2 = q, shape3 = 1,
                               scale = b)
    sum((Fmodel - Fhat)^2)
  }
  try(solution <- optim(c(1, 1), f, lower = c(0.01, 0.01),
                        method = "L-BFGS-B")$par)
  return(c(shape1 = solution[1], shape2 = solution[2],
           scale = mu * (solution[1] + solution[2]) / solution[1]))
}


beta_params <- as.data.frame(t(mapply(fit_3_param_beta,
                                      Fhat = land_cdf_list,
                                      mu = mean_land_owned)))
names(beta_params) <- c('p', 'q', 'H_bar')
omega_n <- with(cross_section_raw, landless_families / num_families)

cross_section <- data.frame(municipality = cross_section_raw$municipality,
                            num_families = cross_section_raw$num_families,
                                beta_params, omega_n)

# Extract 1993 population from panel_raw and merge with cross_section
popn1993 <- subset(panel_raw, year == 1996)[,c("municipality",
                                               "Total_Poblacion1993")]
names(popn1993) <- c('municipality', 'popn1993')
cross_section <- merge(cross_section, popn1993)




panel <- panel_raw[,c("municipality",
                      "year",
                      "ac_vcivilparas_UR", # Cumulative paramilitary violence
                      "vcivilparas_UR", # Flow of paramilitary violence
                      "ac_vcivilnotparas_UR", # Cumulative non-paramilitary violence
                      "vcivilnotparas_UR", # Flow of non-paramilitary violence
                      "desplazados_paras_AS", # Five measures of displacement
                      "desplazados_CODHES",
                      "desplazados_total_RUV",
                      "desplazados_CEDE",
                      "desplazados_JYP",
                      "Total_Poblacion1993")]

names(panel) <- c("municipality",
                  "year",
                  "V_cum",
                  "V_flow",
                  "placeboV_cum",
                  "placeboV_flow",
                  "D_AS",
                  "D_CODHES",
                  "D_RUV",
                  "D_CEDE",
                  "D_JYP",
                  "popn1993")



cum_violence <- subset(panel, year == max(panel$year))[,c("municipality", "V_cum")]
cross_section <- tibble::as_tibble(merge(cross_section, cum_violence))
names(cross_section[, 'V_cum']) <- 'VTotal'

# Keep only municipalities with at least 100 landowners
keep_me <- n_landholders >= 100
keep_municipalities <- cross_section[keep_me,]$municipality
cross_section <- subset(cross_section, municipality %in% keep_municipalities)
panel <- subset(panel, municipality %in% keep_municipalities)
rm(keep_me)

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

ihs <- function(x) {
  y <- log(x + sqrt(x ^ 2 + 1))
  return(y)
}

# Transform some variables to a more sensible scale
drug_routes <- 1 * (covariates$drug_routes > 1) # at least 1km of drug routes?
mines <- 1 * (covariates$mines > 0) # any mines?
oil <- 1 * (covariates$oil > 0) # any oil production?
guerrilla <- ihs(covariates$guerrilla)
grazing <- ihs(covariates$grazing)
grasses <- ihs(covariates$grasses)
offices <- ihs(covariates$offices)
bureaucracy <- ihs(covariates$bureaucracy)

# Overwrite original, untransformed variables
covariates$drug_routes <- drug_routes
covariates$mines <- mines
covariates$oil <- oil
covariates$guerrilla <- guerrilla
covariates$grazing <- grazing
covariates$grasses <- grasses
covariates$offices <- offices
covariates$bureaucracy <- bureaucracy

rm(drug_routes, mines, oil, guerrilla, grazing, grasses, offices, bureaucracy)

# Center and standardize everything except the dummy variables
standardize_me <- setdiff(names(covariates), c('municipality',
                                               'coffee',
                                               'drug_routes',
                                               'mines',
                                               'oil'))

covariates[,standardize_me] <- scale(covariates[,standardize_me])
covariates <- subset(covariates, municipality %in% keep_municipalities)

devtools::use_data(covariates, overwrite = TRUE)
devtools::use_data(cross_section, overwrite = TRUE)
devtools::use_data(panel, overwrite = TRUE)
devtools::use_data(obsZ, overwrite = TRUE)
devtools::use_data(land_parameters, overwrite = TRUE)
devtools::use_data(Vcum, overwrite = TRUE)
devtools::use_data(Vcum_pop, overwrite = TRUE)

# clean up
rm(list = ls())















