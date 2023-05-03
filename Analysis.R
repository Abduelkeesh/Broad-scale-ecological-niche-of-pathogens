# bringing functions from repository
source("https://raw.githubusercontent.com/marlonecobos/host-pathogen/main/Scripts/niche_signal.R")
source("https://raw.githubusercontent.com/marlonecobos/host-pathogen/main/Scripts/plot_niche_signal.R")

# load data
load("Data/all_neon_env_data.RData")

# keep ticks and pathogens that are relevant
sum_tests <- table(neon_env$testPathogenName, neon_env$tick_name)

sum_tests <- data.frame(Pathogen = rownames(sum_tests), 
                        A_americanum = sum_tests[, 1],
                        I_scapularis = sum_tests[, 2])

## positives for pathognes
neon_env_p <- neon_env[neon_env$testResult == 1, ]

sum_pos <- table(neon_env_p$testPathogenName, neon_env_p$tick_name)


sum_pos <- data.frame(Pathogen = rownames(sum_pos), A_americanump = sum_pos[, 1],
                      I_scapularisp = sum_pos[, 2])


## getting intersect of tests and positives
sum_comb <- merge(sum_tests, sum_pos, by = "Pathogen")

sum_comb <- sum_comb[-7, ]

## get relevant ones
a_ame_rel <- sum_comb[sum_comb$A_americanum > 2000 & 
                        sum_comb$A_americanump >=1, 1]
i_sca_rel <- sum_comb[sum_comb$I_scapularis  > 2000 & 
                         sum_comb$I_scapularisp >=1, 1]


# analyses
## relevant variables
varnames <- colnames(neon_env)[12:14]

## sets of variables
sets <- list(set_1 = varnames[1:2], set_2 = varnames[c(1, 3)],
             set_3 = varnames[2:3], set_4 = varnames[1:3])

## saving objects to use them later
save(varnames, sets, a_ame_rel, i_sca_rel, sum_comb,
     file = "Data/in_objects_analysis_plots.RData")

## analyses for A. americanum
### multivariate
a_ame_mult <- lapply(a_ame_rel, function(x) {
  part_data <- neon_env[neon_env$testPathogenName == x & 
                          neon_env$tick_name == "A_americanum", ]
  
  list_aame_permanova <- lapply(1:length(sets), function(y) {
    niche_signal(data = part_data, condition = "testResult", 
                 variables = sets[[y]], method = "permanova",
                 iterations = 1000, verbose = FALSE)
  })
  
  names(list_aame_permanova) <- names(sets)
  print(x)
  
  list_aame_permanova
})

names(a_ame_mult) <- gsub(" ", "_", a_ame_rel)

save(a_ame_mult, file = "Results/aame_permanova.RData")


### univariate
a_ame_univ <- lapply(a_ame_rel, function(x) {
  part_data <- neon_env[neon_env$testPathogenName == x & 
                          neon_env$tick_name == "A_americanum", ]
  
  aame_all_univariate <- lapply(varnames, function(y) {
    niche_signal(data = part_data, condition = "testResult", 
                 variables = y, method = "univariate",
                 iterations = 1000, verbose = FALSE)
  })
  
  names(aame_all_univariate) <- varnames
  
  aame_all_univariate
})

names(a_ame_univ) <- gsub(" ", "_", a_ame_rel)

save(a_ame_univ, file = "Results/aame_univariate.RData")

## analyses for I. scapularis
### multivariate
i_sca_mult <- lapply(i_sca_rel, function(x) {
  part_data <- neon_env[neon_env$testPathogenName == x & 
                          neon_env$tick_name == "I_scapularis", ]
  
  list_aame_permanova <- lapply(1:length(sets), function(y) {
    niche_signal(data = part_data, condition = "testResult", 
                 variables = sets[[y]], method = "permanova",
                 iterations = 1000, verbose = FALSE)
  })
  
  names(list_aame_permanova) <- names(sets)
  print(x)
  
  list_aame_permanova
})

names(i_sca_mult) <- gsub(" ", "_", i_sca_rel)

save(i_sca_mult, file = "Results/iisca_multivariate.RData")

### univariate
i_sca_univ <- lapply(i_sca_rel, function(x) {
  part_data <- neon_env[neon_env$testPathogenName == x & 
                          neon_env$tick_name == "I_scapularis", ]
  
  aame_all_univariate <- lapply(varnames, function(y) {
    niche_signal(data = part_data, condition = "testResult", 
                 variables = y, method = "univariate",
                 iterations = 1000, verbose = FALSE)
  })
  
  names(aame_all_univariate) <- varnames
  
  aame_all_univariate
})

names(i_sca_univ) <- gsub(" ", "_", i_sca_rel)

save(i_sca_univ, file = "Results/iisca_univariate.RData")