# Function to calculate standardized bias
std_bias <- function(df, variavel_interesse,
                     variavel_tratamento = "mahog_area") {
  if (isTRUE(all.equal(df_probit,df)) & variavel_interesse == "distance") {
    variavel_interesse <- ".fitted"
  }
  treat_vec <- as.logical(as.vector(as.data.frame(df)[, variavel_tratamento]))
  int_vec <- as.vector(as.data.frame(df)[, variavel_interesse])
  treated <- int_vec[treat_vec]
  untreated <- int_vec[!treat_vec]
  100 * (mean(treated) - mean(untreated)) / 
    sqrt(0.5 * (var(treated) + var(untreated)))
}

# Function to create dataframe for analyzing changes in standardized bias
# (before matching vs. after matching)
analisar_std_bias <- function (df_unmatched, df_matched) {
  
  covariates <- c("ln_gdp_pc", "gdpfr_ag", "area_plant2", "und5_m", "heart_m",
                  "infecc_m", "neop_m", "suic_m", "traff_m", "pol_m", "hom_tx")
  
  tibble(
    covariada = c("propensity_score", covariates),
    
    std_bias_unmatched = 
      purrr::map_dbl(.x = c("distance", covariates),
                     ~ std_bias(df = df_unmatched, variavel_interesse = .x)),
    
    std_bias_matched = 
      purrr::map_dbl(.x = c("distance", covariates),
                     ~ std_bias(df = df_matched, variavel_interesse = .x)),
    
    reduction_perc = 
      (abs(std_bias_unmatched) - abs(std_bias_matched)) / abs(std_bias_unmatched) * 100
  )
}
