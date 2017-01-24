cross_section <- haven::read_dta("data-raw/Cross_Section.dta")
panel_D_V <- haven::read_dta("data-raw/Panel_D_V.dta")
survey <- haven::read_dta("data-raw/Survey.dta")

# Construct land distribution from raw data on land ownership (in bins)
landowner_cols <- sapply(names(cross_section), substr, 1, 8) == 'landown_'
landowner_bins <- cross_section[, landowner_cols]
names(landowner_bins) <- c('[0,1)', '[1,3)', '[3,5)', '[5,10)', '[10,15)',
                         '[15,20)', '[20,50)', '[50,100)', '[100,200)',
                         '[200,500)', '[500,1000)', '[1000,2000)', '>=2000')
rm(landowner_cols)

h <- function(x){
  x <- x / sum(x)
  n <- length(x)
  cumsum(x)[1:(n-1)]
}
land_cdfs <- t(apply(landowner_bins, 1, h))
colnames(land_cdfs) <- paste0('x=', c(1, 3, 5, 10, 15, 20, 50, 100, 200, 500,
                                     1000, 2000))
rm(h, landowner_bins)

land_cdf_list <- lapply(seq_len(nrow(land_cdfs)), function(i) land_cdfs[i,])
rm(land_cdfs)

fit3_B1 <- function(Fhat, mu){
  x_full <- c(1, 3, 5, 10, 15, 20, 50, 100, 200, 500, 1000, 2000)
  x <- x_full[1:length(Fhat)]
  f <- function(params){
    p <- params[1]
    q <- params[2]
    b <- mu * (p + q) / p
    Fmodel <- actuar::pgenbeta(x, shape1 = p, shape2 = q, shape3 = 1, scale = b)
    sum((Fmodel - Fhat)^2)
  }
  try(solution <- optim(c(1, 1), f, lower = c(0.01, 0.01),
                        method = "L-BFGS-B")$par)
  return(c(shape1 = solution[1], shape2 = solution[2],
           scale = mu * (solution[1] + solution[2]) / solution[1]))
}

B1_params3 <- as.data.frame(t(mapply(fit3_B1, Fhat = land_cdf_list,
                                     mu = data_means)))

