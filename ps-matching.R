library(tidyverse)
library(stargazer)
library(gridExtra)
library(DescTools)
library(MatchIt)
library(haven)
library(lfe)

# Loading auxiliary functions
source('ps-helpers.R')

# Loading clean and complete DD dataframe
df_mahogany <- readRDS('df_mahogany.rds')

# PREPARING DATA FOR PROPENSITY SCORE ESTIMATION ------------------------------
# Selecting variables and filtering years (1995 to 1998)
df_probit <- df_mahogany %>% 
  select(uf, code, year, mahog_area, hom_tx, area_plant, ln_gdp_pc, gdpfr_ag,
         heart, heart_m, infecc, infecc_m, neop, neop_m, und5_mort, und5_m, 
         traff, traff_m, suic, suic_m, pol_deaths, pol_m,
         avg_pop, pop_und5, pop) %>%
  filter(year %in% 1995:1998)

# Replacing dependent variables with averages
# Note: for ln_gdp_pc and gdpfr_ag there is only one available value
df_probit <- df_probit %>% 
  # ln_gdp_pc and gdpfr_ag: 1996 values for all municipalities
  mutate(ln_gdp_pc = rep( (df_probit$ln_gdp_pc[df_probit$year == 1996]), each = 4),
         gdpfr_ag = rep( (df_probit$gdpfr_ag[df_probit$year == 1996]), each = 4) ) %>% 
  # Averages
  group_by(code) %>% 
  mutate_at(vars(-uf, -code, -year), mean) %>% 
  ungroup() %>% 
  # Keeping only one observation per municipality
  filter(year == 1996) %>% 
  # Planted area: we will use fractions instead of percentages
  mutate(area_plant2 = area_plant / 100) %>% 
  # Removing observation with missing data
  filter(!is.na(und5_m), !is.infinite(und5_m), !is.na(area_plant) )

# PS ESTIMATION - PROBIT ------------------------------------------------------
probit <- glm(mahog_area ~  ln_gdp_pc + gdpfr_ag + area_plant2 + und5_m + 
                heart_m + infecc_m + neop_m + suic_m + traff_m + pol_m,
              family = binomial("probit"), data = df_probit) ; summary(probit)

# Including propensity scores and inverse-probability weights in df_probit
df_probit <- df_probit %>% 
  mutate(.fitted = predict.glm(probit, type = "response"),
         ip_weight = case_when( mahog_area == 1 ~ 1/.fitted,
                                mahog_area == 0 ~ 1/ (1-.fitted)))

# Including propensity scores and inverse-probability weights in df_mahogany
df_mahogany <- df_mahogany %>% 
  left_join(df_probit %>% 
              select(code, ip_weight, .fitted),
            by = c("code"))

# Confusion matrix to assess prediction quality
conf_matrix <- table(df_probit$mahog_area,
                     df_probit$.fitted > mean(df_probit$.fitted))

conf_matrix_prop <- prop.table(conf_matrix) %>% round(digits = 2)

conf_matrix
conf_matrix_prop

# MATCHING ON PS SCORES -------------------------------------------------------
# NN: 1 neighboor, without replacement, no caliper
match_near <- 
  matchit(formula = mahog_area ~  ln_gdp_pc + gdpfr_ag + area_plant2 + und5_m +
            heart_m + infecc_m + neop_m + suic_m + traff_m + pol_m,
          distance = "probit", method = "nearest", data = df_probit)

summary(match_near)

# Storing results in a dataframe
df_match_near <- match.data(match_near) %>% as_tibble()

# NN: broad caliper (0.1)
match_caliper_grande <- 
  matchit(formula = mahog_area ~  ln_gdp_pc + gdpfr_ag + area_plant2 + und5_m + 
            heart_m + infecc_m + neop_m + suic_m + traff_m + pol_m,
          distance = "probit", method = "nearest",
          caliper = 0.1, data = df_probit)

summary(match_caliper_grande)

# Storing results in a dataframe
df_match_caliper_grande <- match.data(match_caliper_grande) %>% as_tibble()

# NN: narrow caliper (0.01)
match_caliper_pequeno <- 
  matchit(formula = mahog_area ~  ln_gdp_pc + gdpfr_ag + area_plant2 + und5_m +
            heart_m + infecc_m + neop_m + suic_m + traff_m + pol_m,
          distance = "probit", method = "nearest",
          caliper = 0.01, data = df_probit)

summary(match_caliper_pequeno)

# Storing results in a dataframe
df_match_caliper_pequeno <- match.data(match_caliper_pequeno) %>% as_tibble()

# Kernel matching
# Implemented in Stata using psmatch2, the same databse (df_probit),
# and same probit specification

# Saving df_probit as a dta file to perform matching in Stata
haven::write_dta(df_probit, 'df_probit_unmatched.dta')

# Loading dta file with calculated kernel weights
df_kernel <- haven::read_dta("df_probit_matched.dta") %>% 
  mutate_at(vars(uf, code, year, `_treated`, `_support`), as_factor)

# DF containing only municipalities and respective weights
df_weights <- df_kernel %>%
  rename(peso_kernel = `_weight`) %>% 
  select(code, peso_kernel) %>%
  arrange(code)

# Joining weight data (same weights for all years)
df_probit <- left_join(df_probit, df_weights, by = c("code") )

# Removing unmatched observations
df_match_kernel <- df_probit %>% 
  filter(!is.na(peso_kernel)) 

# ASSESSING COMMON SUPPORT AND BALANCE ----------------------------------------

# Histograms ------------------------------------------------------------------
## Before matching
ph1 <- ggplot(df_probit,
              aes(x = .fitted,
                  col = factor(mahog_area),
                  fill = factor(mahog_area))) +
  geom_histogram(position = "identity", bins = 15) +
  facet_wrap(~factor(mahog_area, labels = c("Controle", "Tratamento"))) +
  theme(legend.position = "none") +
  xlab("Propensity score") +
  ylab("Numero de municipios") +
  ggtitle("A. Amostra nao-pareada") ; ph1

# NN1, simples, sem caliper
ph2 <- ggplot(df_match_near,
              aes(x = distance,
                  col = factor(mahog_area),
                  fill = factor(mahog_area))) +
  geom_histogram(position = "identity", bins = 15) +
  facet_wrap(~factor(mahog_area, labels = c("Controle", "Tratamento"))) +
  theme(legend.position = "none") +
  xlab("Propensity score") +
  ylab("Numero de municipios") +
  ggtitle("B. NN simples") ; ph2

# Caliper permissivo, 0.1
ph3 <- ggplot(df_match_caliper_grande,
              aes(x = distance,
                  col = factor(mahog_area),
                  fill = factor(mahog_area))) +
  geom_histogram(position = "identity", bins = 15) +
  facet_wrap(~factor(mahog_area, labels = c("Controle", "Tratamento"))) +
  theme(legend.position = "none") +
  xlab("Propensity score") +
  ylab("Numero de municipios") +
  ggtitle("C. NN com caliper de 0,1") ; ph3

# Caliper restritivo, 0.01
ph4 <- ggplot(df_match_caliper_pequeno,
              aes(x = distance,
                  col = factor(mahog_area),
                  fill = factor(mahog_area))) +
  geom_histogram(position = "identity", bins = 15) +
  facet_wrap(~factor(mahog_area, labels = c("Controle", "Tratamento"))) +
  theme(legend.position = "none") +
  xlab("Propensity score") +
  ylab("Numero de municipios") +
  ggtitle("D. NN com caliper de 0,01") ; ph4

# Kernel (unweighted)
ph_kernel1 <- 
  df_probit %>% 
  filter(!is.na(peso_kernel)) %>%
  ggplot(aes(x = .fitted,
             col = factor(mahog_area),
             fill = factor(mahog_area))) +
  geom_histogram(position = "identity", bins = 15) +
  facet_wrap(~factor(mahog_area, labels = c("Controle", "Tratamento"))) +
  theme(legend.position = "none") +
  xlab("Propensity score") +
  ylab("Numero de municipios") +
  ggtitle("Amostra pareada por kernel, sem considerar os pesos") +
  ylim(c(0,130)) ; ph_kernel1

# Kernel (weighted)
ph_kernel2 <- 
  df_probit %>% 
  filter(!is.na(peso_kernel)) %>% 
  ggplot(aes(x = .fitted,
             col = factor(mahog_area),
             fill = factor(mahog_area),
             weights = peso_kernel)) +
  geom_histogram(position = "identity", bins = 15) +
  facet_wrap(~factor(mahog_area, labels = c("Controle", "Tratamento"))) +
  theme(legend.position = "none") +
  xlab("Propensity score") +
  ylab("Numero de municipios") +
  ggtitle("Amostra pareada por kernel, considerando os pesos") +
  ylim(c(0,130)) ; ph_kernel2

# Putting all histograms together (except kernel)
histogramas <- gridExtra::grid.arrange(ph1, ph2, ph3, ph4, ncol = 2)

# Density plots ---------------------------------------------------------------
# General settings
density_plot <- ggplot(mapping = aes(x = .fitted, fill = factor(mahog_area))) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 4.5)) +
  xlab("Propensity Score") +
  ylab("Densidade")

# Before matching
pd1 <- density_plot +
  geom_density(data = df_probit, alpha = 0.4) +
  ggtitle("A. Sem matching") +
  scale_fill_discrete(
    name = "Grupo",
    labels = c("Controle (459 obs.)", "Tratamento (148 obs.)")
    ) ; pd1

# NN1, no caliper
pd2 <- density_plot +
  geom_density(data = df_match_near, alpha = 0.4) +
  ggtitle("B. Vizinho mais proximo, sem caliper") +
  scale_fill_discrete(
    name = "Grupo",
    labels = c("Controle (148 obs.)", "Tratamento (148 obs.)")
    ) ; pd2

# Broad caliper, 0.1
pd3 <- density_plot +
  geom_density(data = df_match_caliper_grande, alpha = 0.4) +
  ggtitle("C. Vizinho mais proximo, com caliper de 0,1") +
  scale_fill_discrete(
    name = "Grupo",
    labels = c("Controle (137 obs.)", "Tratamento (137 obs.)")
    ) ; pd3

# Narrow caliper, 0.01
pd4 <- density_plot +
  geom_density(data = df_match_caliper_pequeno, alpha = 0.4) +
  ggtitle("D. Vizinho mais proximo, com caliper de 0,01") +
  scale_fill_discrete(
    name = "Grupo",
    labels = c("Controle (93 obs.)", "Tratamento (93 obs.)")
    ) ; pd4

# Kernel (weighted)
pd5 <- density_plot +
  geom_density(data = df_probit %>% filter(!is.na(peso_kernel)),
               mapping = aes(weight = peso_kernel/147),
               alpha = 0.4) +
  ggtitle("E. Kernel (ponderado)") +
  scale_fill_discrete(
    name = "Grupo",
    labels = c("Controle (457 obs.)", "Tratamento (147 obs.)")
    ) ; pd5

# Putting all plots together
gridExtra::grid.arrange(pd1, pd2, pd3, pd4, pd5, ncol = 1)


# ASSESSING MATCHING QUALITY --------------------------------------------------
# Covariates + hom_tx
covariates <- c("ln_gdp_pc", "gdpfr_ag", "area_plant2", "und5_m", "heart_m",
                "infecc_m", "neop_m", "suic_m", "traff_m", "pol_m", "hom_tx")

# 1) Differences in means: before and after matching
tableone::CreateTableOne(vars = c(".fitted", covariates),
                         strata = "mahog_area", data = df_probit)
tableone::CreateTableOne(vars = c("distance", covariates),
                         strata = "mahog_area", data = df_match_near)
tableone::CreateTableOne(vars = c("distance", covariates),
                         strata = "mahog_area", data = df_match_caliper_grande)
tableone::CreateTableOne(vars = c("distance", covariates),
                         strata = "mahog_area", data = df_match_caliper_pequeno)

# 2) Changes in standardized bias
tabela_std_bias_near <- analisar_std_bias(df_unmatched = df_probit,
                                          df_matched = df_match_near)
tabela_std_bias_caliper_grande <- analisar_std_bias(df_unmatched = df_probit,
                                                    df_matched = df_match_caliper_grande)
tabela_std_bias_caliper_pequeno <- analisar_std_bias(df_unmatched = df_probit,
                                                     df_matched = df_match_caliper_pequeno)

comparacao_std_bias <- 
  tibble(variavel = tabela_std_bias_near$covariada,
         near = tabela_std_bias_near$reduction_perc,
         caliper_grande = tabela_std_bias_caliper_grande$reduction_perc,
         caliper_pequeno = tabela_std_bias_caliper_pequeno$reduction_perc) %>%
  mutate_at(vars(c("near", "caliper_grande", "caliper_pequeno")),
            round, digits = 2)

comparacao_std_bias

# 3) Pseudo R2

# Reestimating probit using NN1-matched sample (no caliper)
probit_matched_near <- 
  glm(mahog_area ~ ln_gdp_pc + gdpfr_ag + area_plant2 + und5_m + 
        heart_m + infecc_m + neop_m + suic_m + traff_m + pol_m,
      family = binomial("probit"), data = df_match_near)

# Reestimating probit using broad caliper-matched sample (0.1)
probit_matched_caliper_grande <- 
  glm(mahog_area ~  ln_gdp_pc + gdpfr_ag + area_plant2 + und5_m + 
        heart_m + infecc_m + neop_m + suic_m + traff_m + pol_m,
      family = binomial("probit"), data = df_match_caliper_grande)

# Reestimating probit using narrow caliper-matched sample (0.01)
probit_matched_caliper_pequeno <- 
  glm(mahog_area ~  ln_gdp_pc + gdpfr_ag + area_plant2 + und5_m + 
        heart_m + infecc_m + neop_m + suic_m + traff_m + pol_m,
      family = binomial("probit"), data = df_match_caliper_pequeno)

# Reestimating probit using kernel-matched sample (weighted)
probit_matched_kernel <- 
  glm(mahog_area ~ ln_gdp_pc + gdpfr_ag + area_plant2 + und5_m + 
        heart_m + infecc_m + neop_m + suic_m + traff_m + pol_m,
      family = binomial("probit"), data = df_match_kernel,
      weights = peso_kernel)

# Table showing each Pseudo-R2
comparacao_pseudoR2 <- 
  tibble(Amostra_utilizada = 
           c("Original", "Pareada NN", "Pareada NN caliper 0,1",
             "Pareada NN caliper 0,01", "Pareada Kernel"),
         PseudoR2 =
           map_dbl(.x = list(probit, probit_matched_near,
                             probit_matched_caliper_grande,
                             probit_matched_caliper_pequeno,
                             probit_matched_kernel),
                   .f = DescTools::PseudoR2))

comparacao_pseudoR2


# POST-MATCHING ANALYSIS - PART 1 ---------------------------------------------
# Building 2-period DD database
# (pre and post-treatment averages, considering 1999 and 2002)
df_2p <- df_mahogany %>%
  # Only municipalities used in PS estimation
  filter(code %in% df_probit$code) %>% 
  # Dummies: before/after each treatment
  mutate(apos1999 = if_else(year %in% 1999:2013, 1, 0),
         apos2002 = if_else(year %in% 2002:2013, 1, 0) ) %>% 
  # Tomando media da taxa de homicidios pre e pos 1999, para cada municipio
  group_by(code, apos1999) %>% 
  mutate(hom_tx_med1999 = mean(hom_tx)) %>% 
  ungroup() %>% 
  # Tomando media da taxa de homicidios pre e pos 2002, para cada municipio
  group_by(code, apos2002) %>% 
  mutate(hom_tx_med2002 = mean(hom_tx)) %>% 
  ungroup() %>% 
  # Filtrando para considerar apenas uma observacoes pre e pos tratamento
  filter(year %in% c(1995, 2013))

# Defining subsets
## Only municipalities matched by NN1
df_2p_near <- df_2p %>% 
  filter(code %in% df_match_near$code)
## Only municipalities matched by broad caliper
df_2p_caliper_grande <- df_2p %>% 
  filter(code %in% df_match_caliper_grande$code)
## Only municipalities matched by narrow caliper
df_2p_caliper_pequeno <- df_2p %>% 
  filter(code %in% df_match_caliper_pequeno$code)
## Only municipalities matched by kernel
df_2p_kernel <- df_2p %>% 
  filter(code %in% df_match_kernel$code) %>% 
  mutate(peso_kernel = rep(df_match_kernel$peso_kernel, each = 2) )

# Simple DD
## 1999
summary(dd_1999 <- lm(hom_tx_med1999 ~ mahog_area*apos1999,
                      data = df_2p))

summary(dd_1999_near <- lm(hom_tx_med1999 ~ mahog_area*apos1999,
                           data = df_2p_near))

summary(dd_1999_caliper_grande <- lm(hom_tx_med1999 ~ mahog_area*apos1999,
                                     data = df_2p_caliper_grande))

summary(dd_1999_caliper_pequeno <- lm(hom_tx_med1999 ~ mahog_area*apos1999,
                                      data = df_2p_caliper_pequeno))

summary(dd_1999_kernel <- lm(hom_tx_med1999 ~ mahog_area*apos1999,
                             data = df_2p_kernel, weights = peso_kernel))
## 2002
summary(dd_2002 <- lm(hom_tx_med2002 ~ mahog_area*apos2002,
                      data = df_2p))

summary(dd_2002_near <- lm(hom_tx_med2002 ~ mahog_area*apos2002,
                           data = df_2p_near))

summary(dd_2002_caliper_grande <- lm(hom_tx_med2002 ~ mahog_area*apos2002,
                                     data = df_2p_caliper_grande))

summary(dd_2002_caliper_pequeno <- lm(hom_tx_med2002 ~ mahog_area*apos2002,
                                      data = df_2p_caliper_pequeno))

summary(dd_2002_kernel <- lm(hom_tx_med2002 ~ mahog_area*apos2002,
                             data = df_2p_kernel, weights = peso_kernel))

# Simple mean comparison: avg hom_tx after treatment vs maho_area
## 1999
summary(comp_1999 <- lm(formula = hom_tx_med1999 ~ mahog_area,
                        data = df_2p,
                        subset = (apos1999 == 1)))

summary(comp_1999_near <- lm(formula = hom_tx_med1999 ~ mahog_area,
                             data = df_2p_near,
                             subset = (apos1999 == 1)))

summary(comp_1999_caliper_grande <- lm(formula = hom_tx_med1999 ~ mahog_area,
                                       data = df_2p_caliper_grande,
                                       subset = (apos1999 == 1)))

summary(comp_1999_caliper_pequeno <- lm(formula = hom_tx_med1999 ~ mahog_area,
                                        data = df_2p_caliper_pequeno,
                                        subset = (apos1999 == 1)))

summary(comp_1999_kernel <- lm(formula = hom_tx_med1999 ~ mahog_area,
                               data = df_2p_kernel,
                               subset = (apos1999 == 1),
                               weights = peso_kernel))
## 2002
summary(comp_2002 <- lm(formula = hom_tx_med2002 ~ mahog_area,
                        data = df_2p,
                        subset = (apos2002 == 1)))

summary(comp_2002_near <- lm(formula = hom_tx_med2002 ~ mahog_area,
                             data = df_2p_near,
                             subset = (apos2002 == 1)))

summary(comp_2002_caliper_grande <- lm(formula = hom_tx_med2002 ~ mahog_area,
                                       data = df_2p_caliper_grande,
                                       subset = (apos2002 == 1)))

summary(comp_2002_caliper_pequeno <- lm(formula = hom_tx_med2002 ~ mahog_area,
                                        data = df_2p_caliper_pequeno,
                                        subset = (apos2002 == 1)))

summary(comp_2002_kernel <- lm(formula = hom_tx_med2002 ~ mahog_area,
                               data = df_2p_kernel,
                               subset = (apos2002 == 1),
                               weights = peso_kernel))


# POST-MATCHING ANALYSIS - PART 2 ---------------------------------------------
## Including kernel weights in the main dataframe (weights are year-invariant)
df_mahogany <- left_join(df_mahogany, df_weights, by = c("code") )

## Defining subsets
# Only muns used in the probit (some were removed due to missing data)
df_mahogany_probit <- df_mahogany %>% 
  filter(code %in% df_probit$code)
# Only NN1-matched municipalities
df_mahogany_near <- df_mahogany %>% 
  filter(code %in% df_match_near$code)
# Only broad caliper-matched municipalities
df_mahogany_caliper_grande <- df_mahogany %>% 
  filter(code %in% df_match_caliper_grande$code)
# Only narro caliper-matched municipalities
df_mahogany_caliper_pequeno <- df_mahogany %>% 
  filter(code %in% df_match_caliper_pequeno$code)
# Only kernel-matched municipalities
df_mahogany_kernel <- df_mahogany %>% 
  filter(!is.na(peso_kernel)) %>%
  # avg_pop X kernel weights: updating regression weights in order to
  # incorporate kernel matching results
  mutate(avg_pop = avg_pop * peso_kernel)

# Benchmark specifications using matched samples
# (same specifications as those of columns 1-4 from table 2)

df_nested <- 
  tibble(
    referencia = c("df_mahogany_probit",
                   "df_mahogany_near",
                   "df_mahogany_caliper_grande",
                   "df_mahogany_caliper_pequeno",
                   "df_mahogany_kernel"),
    data = list(df_mahogany_probit,
                df_mahogany_near,
                df_mahogany_caliper_grande,
                df_mahogany_caliper_pequeno,
                df_mahogany_kernel)) %>% 
    
  # 1) Simple, uncontrolled DD
  mutate(m1 = map(.x = data, ~ felm(formula = hom_tx ~ treat1 + treat2 + treat3 
                                    + year | code | 0 | code,
                                    data = . , weights = .$avg_pop) ),
         s1 = map(.x = m1, ~ summary(.)),
         c1 = map(.x = s1, ~ coef(.)[1:3,] )) %>% 
  
  # 2) State/year fixed-effects
  mutate(m2 = map(.x = data, ~ felm(formula = hom_tx ~ treat1 + treat2 + treat3 
                                    + uf*year + year | code | 0 | code,
                                    data = . , weights = .$avg_pop)),
         s2 = map(.x = m2, ~ summary(.)),
         c2 = map(.x = s2, ~ coef(.)[1:3,])) %>% 
  
  # 3) State/year fixed-effects and baseline charct. X year
  mutate(m3 = map(.x = data,
                  ~ felm(formula = hom_tx ~ treat1 + treat2 + treat3 + uf*year +
                           year*hom_base + year*area_base + year*lngdp_base + 
                           year*gdpag_base + year*pol_base + year*und5_base + 
                           year*infecc_base + year*heart_base + 
                           year*neop_base + year*suic_base + 
                           year*traff_base + year | code | 0 | code,
                         data = . , weights = .$avg_pop) ),
         s3 = map(.x = m3, ~ summary(.)),
         c3 = map(.x = s3, ~ coef(.)[1:3,] )) %>% 
  
  # 4) State/year fixed-effects + int. between treatment and linear trend
  mutate(m4 = map(.x = data,
                  ~ felm(formula = hom_tx ~ treat1 + int1 + treat2 + int2 + 
                           treat3 + int3 + uf*year + year | code | 0 | code,
                         data = . , weights = .$avg_pop) ),
         s4 = map(.x = m4, ~ summary(.)),
         c4 = map(.x = s4, ~ coef(.)[1:3,] ))

# Results

# Loading  Chimeli & Soares' benchmark results for comparison
load('benchmark-models.RData')

# 1) Simple, uncontrolled DD
stargazer(t2c1, df_nested$m1[[2]], df_nested$m1[[3]],
          df_nested$m1[[4]], df_nested$m1[[5]],
          omit = c("year", "uf", "base"),
          omit.stat = c("ser", "adj.rsq"),
          decimal.mark = ",",
          digit.separator = ".", type = "text")

# 2) State/year fixed-effects
stargazer(t2c2, df_nested$m2[[2]], df_nested$m2[[3]],
          df_nested$m2[[4]], df_nested$m1[[5]],
          omit = c("year", "uf", "base"),
          omit.stat = c("ser", "adj.rsq"),
          decimal.mark = ",",
          digit.separator = ".", type = "text")

# 3) State/year fixed-effects and baseline charct. X year
stargazer(t2c3, df_nested$m3[[2]], df_nested$m3[[3]],
          df_nested$m3[[4]], df_nested$m1[[5]],
          omit = c("year", "uf", "base"),
          omit.stat = c("ser", "adj.rsq"),
          decimal.mark = ",",
          digit.separator = ".", type = "text")

# 4) State/year fixed-effects + int. between treatment and linear trend
stargazer(t2c4, df_nested$m3[[2]], df_nested$m3[[3]],
          df_nested$m3[[4]], df_nested$m1[[5]],
          omit = c("year", "uf", "base"),
          omit.stat = c("ser", "adj.rsq"),
          decimal.mark = ",",
          digit.separator = ".", type = "text")

# PARALLEL TRENDS -------------------------------------------------------------

# Dataframes with annual avg homicide rate for each group
# All observations, unmathced
df_par_mahogany <- df_mahogany %>% 
  group_by(mahog_area, year) %>% 
  summarise(hom = mean(hom_tx))
# Unmatched observations used to estimate the probit
df_par_probit <- df_mahogany_probit %>% 
  group_by(mahog_area, year) %>% 
  summarise(hom = mean(hom_tx))
# NN1-matched observations
df_par_near <- df_mahogany_near %>% 
  group_by(mahog_area, year) %>% 
  summarise(hom = mean(hom_tx))
# Broad caliper-matched observations
df_par_caliper_grande <- df_mahogany_caliper_grande %>% 
  group_by(mahog_area, year) %>% 
  summarise(hom = mean(hom_tx))
# Narrow caliper-matched observations
df_par_caliper_pequeno <- df_mahogany_caliper_pequeno %>% 
  group_by(mahog_area, year) %>% 
  summarise(hom = mean(hom_tx))
# Kernel-matched observations (weighted average)
df_par_kernel <- df_mahogany_kernel %>% 
  group_by(mahog_area, year) %>% 
  summarise(hom = weighted.mean(x = hom_tx, w = peso_kernel))

# Plot general settings
line_plot <- ggplot(mapping = aes(x = factor(year),
                                  y = hom,
                                  group = mahog_area,
                                  col = factor(mahog_area)) ) +
  xlab("Ano (1995-2013, ultimos digitos)") +
  ylab("Homicidios por 100 mil hab.") +
  scale_color_discrete(name = "Grupo",
                       labels = c("Controle", "Tratamento")) +
  scale_x_discrete(labels = c("95", "96", "97", "98", "99","00", "01",
                              "02", "03", "04", "05", "06", "07", "08",
                              "09", "10", "11", "12", "13") ) +
  theme(legend.position = "bottom")

# Plotting

# All observations, unmatched
par_mahogany <- line_plot +
  geom_line(data = df_par_mahogany, size = 1.5) +
  geom_vline(aes(xintercept = 5), size = 1, linetype = "dashed")

# Probit observations, unmatched
par_probit <- line_plot +
  geom_line(data = df_par_probit, size = 1.5) +
  geom_vline(aes(xintercept = 5 ), size = 1, linetype = "dashed")

# NN1-matched observations
par_near <- line_plot +
  geom_line(data = df_par_near, size = 1.5) +
  geom_vline(aes(xintercept = 5 ), size = 1, linetype = "dashed")

# Broad caliper-matched observations
par_caliper_grande <-  line_plot +
  geom_line(data = df_par_caliper_grande, size = 1.5) +
  geom_vline(aes(xintercept = 5 ), size = 1, linetype = "dashed")

# Narrow caliper-matched observations
par_caliper_pequeno <-  line_plot +
  geom_line(data = df_par_caliper_pequeno, size = 1.5) +
  geom_vline(aes(xintercept = 5 ), size = 1, linetype = "dashed")

# Kernel-matched observations
par_kernel <- line_plot +
  geom_line(data = df_par_kernel, size = 1.5) +
  geom_vline(aes(xintercept = 5 ), size = 1, linetype = "dashed")

## Putting plots in a single image
gridExtra::grid.arrange(par_mahogany, par_probit,
                        par_near, par_caliper_grande,
                        par_caliper_pequeno, par_kernel,
                        ncol = 3)
