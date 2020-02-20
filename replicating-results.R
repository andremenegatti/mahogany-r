library(tidyverse)
library(stargazer)
library(haven)
library(lfe)

# Data wrangling: setting up the DD dataframe ---------------------------------
# Importing data from dta file
df_mahogany <- haven::read_dta("AEJApp-2016-0055-Stata_Dataset.dta")

# Replacing NAs in pop with 0s
df_mahogany <- df_mahogany %>% 
  mutate(pop = replace_na(pop, 0))

# Only municipalities that already existed in 1995
df_mahogany <- df_mahogany %>% 
  filter(!code %in% 
           df_mahogany$code[df_mahogany$pop == 0 & df_mahogany$year == 1995])

# Defining control/auxiliary variables
df_mahogany <- df_mahogany %>% 
  mutate(mahogx_state = if_else(is.na(mahogx_state), 0, mahogx_state),
         otherx_state = if_else(is.na(otherx_state), 0, otherx_state),
         total_mahogx = (mahogx_state + otherx_state) / 1e+6,
         gdp_pc = gdp/pop,
         ln_gdp_pc = log(gdp_pc),
         gdpfr_ag = gdp_ag / gdp,
         area_plant = (plant / area) * 100,
         ln_otherx_state = log(otherx_state),
         trend = year - 1995)

# Defining variables related to mortality and homicides
df_mahogany <- df_mahogany %>% 
  mutate(
    hom = replace_na(hom, 0),
    
    hom_male = if_else((is.na(hom_male) & year > 1995),
                       0, hom_male),
    hom_m_prime = if_else((is.na(hom_m_prime) & year > 1995),
                           0, hom_m_prime ),
    hom_m_sing = if_else((is.na(hom_m_sing) & year > 1995),
                          0, hom_m_sing ),
    hom_m_nothome = if_else((is.na(hom_m_nothome) & year > 1995),
                             0, hom_m_nothome ),
    hom_m_firearm = if_else((is.na(hom_m_firearm) & year > 1995),
                             0, hom_m_firearm ),
    
    heart = if_else((is.na(heart) & year < 2008), 0, heart),
    infecc = if_else((is.na(infecc) & year < 2008), 0, infecc),
    neop = if_else((is.na(neop) & year < 2008), 0, neop ),
    und5_mort = if_else((is.na(und5_mort) & year < 2008), 0, und5_mort),
    traff = if_else((is.na(traff) & year < 2008), 0, traff),
    suic = if_else((is.na(suic) & year < 2008), 0, suic),
    pol_deaths = if_else((is.na(pol_deaths) & year < 2008), 0, pol_deaths),
         
    hom_tx = (hom / pop) * 1e+5,
    male_h = (hom_male / pop_male) * 1e+5,
    male_h_prime = (hom_m_prime / pop_m_prime) * 1e+5,
    male_h_nothome = (hom_m_nothome / pop_m_prime) * 1e+5,
    male_h_sing = (hom_m_sing / pop_m_prime) * 1e+5,
    male_h_firearm = (hom_m_firearm / pop_m_prime) * 1e+5,
    
    und5_m = (und5_mort / pop_und5) * 1e+3,
    heart_m = (heart / pop) * 1e+3,
    infecc_m = (infecc / pop) * 1e+3,
    neop_m = (neop / pop) * 1e+3,
    suic_m = (suic / pop) * 1e+3,
    traff_m = (traff / pop) * 1e+3,
    pol_m = (pol_deaths / pop) * 1e+3,
    pol_deaths_dummy = if_else((pol_deaths > 0 & year < 2008), 1, 0),
    pol_deaths_dummy = if_else((is.na(pol_deaths_dummy) & year < 2008),
                                0, pol_deaths_dummy)
    )

# Simple function to define baseline values
definindo_baseline <- function (variavel, ano_base = 1995) {
  rep( ( variavel[df_mahogany$year == ano_base]), each = 19 )
}

# Defining baseline-constant variables (to interact with time dummies)
df_mahogany <- df_mahogany %>% 
  mutate(
    hom_base = definindo_baseline(variavel = .$hom_tx),
    area_base = definindo_baseline(variavel = .$area_plant),
    lngdp_base = definindo_baseline(variavel = .$ln_gdp_pc, ano = 1996),
    gdpag_base = definindo_baseline(variavel = .$gdpfr_ag, ano = 1996),
    pol_base = definindo_baseline(variavel = .$pol_m),
    und5_base = definindo_baseline(variavel = .$und5_m),
    infecc_base = definindo_baseline(variavel = .$infecc_m),
    heart_base = definindo_baseline(variavel = .$heart_m),
    neop_base = definindo_baseline(variavel = .$neop_m),
    suic_base = definindo_baseline(variavel = .$suic_m),
    traff_base = definindo_baseline(variavel = .$traff_m),
    pop1995 = definindo_baseline(variavel = .$pop)
  )

# Average population for each municipality (for weighting the regressions)
df_mahogany <- df_mahogany %>% 
  group_by(code) %>% 
  mutate(avg_pop = mean(pop)) %>% 
  ungroup()

# Treatment variables
df_mahogany <- df_mahogany %>% 
  mutate(
    treat1 = if_else((mahog_area == 1 & year >= 1999 & year <= 2001), 1, 0),
    treat2 = if_else((mahog_area == 1 & year >= 2002 & year <= 2008), 1, 0),
    treat3 = if_else((mahog_area == 1 & year >= 2009), 1, 0),
    
    treatx_area1 = if_else( (mahog_area == 1 & year >= 1999 & year <= 2001),
                            mahog_exp_pre, 0),
    treatx_area2 = if_else( (mahog_area == 1 & year >= 2002 & year <= 2008),
                            mahog_exp_pre, 0),
    treatx_area3 = if_else( (mahog_area == 1 & year >= 2009),
                            mahog_exp_pre, 0),
    
    treatt_area1 = if_else((mahog_area == 1 & year >= 1999 & year <= 2001),
                            total_mahogx, 0),
    treatt_area2 = if_else((mahog_area == 1 & year >= 2002 & year <= 2008),
                            total_mahogx, 0),
    treatt_area3 = if_else((mahog_area == 1 & year >= 2009),
                            total_mahogx, 0),
    
    trend1 = if_else((year >= 1999 & year <= 2001), year - 1999, 0),
    trend2 = if_else((year >= 2002 & year <= 2008), year - 2002, 0),
    trend3 = if_else((year >= 2009), year - 2009, 0),
    
    int1 = treat1*trend1,
    int2 = treat2*trend2,
    int3 = treat3*trend3
    )

# Placebo treatment
df_mahogany <- df_mahogany %>% 
  mutate(pre = if_else((mahog_area == 1 & year > 1996 & year < 1999), 1, 0))

# Coercing year, uf and municipality codes to factor
df_mahogany <- df_mahogany %>% 
  mutate(year = factor(year),
         code = factor(code),
         uf = factor(uf) )

# Saving in rds format for future use
saveRDS(df_mahogany, 'df_mahogany.rds')

# Creating subsets required to run some specifications ------------------------
# Only municipalities from Pará
df_mahogany_para <- df_mahogany %>% 
  filter(uf == 15) %>% 
  droplevels()

# Only municipalities outside Pará
df_mahogany_not_para <- df_mahogany %>% 
  filter(uf != 15) %>% 
  droplevels()

# Only years for which there is gdp data
df_mahogany_para_gdp <- df_mahogany %>% 
  filter(uf == 15) %>% 
  filter( year %in% c(1996, 1999:2010) ) %>% 
  droplevels()

# Only municipalities in which we could calculate baseline values
df_mahogany_baseline <- df_mahogany %>% 
  filter(code != df_mahogany$code[(df_mahogany$year == 1995 & 
                                     is.na(df_mahogany$plant))],
         !code %in% df_mahogany$code[(df_mahogany$year == 1995 & 
                                        df_mahogany$pop_und5 == 0)]) %>% 
  droplevels()

# Baseline, Pará
df_mahogany_baseline_para <- df_mahogany_baseline %>% 
  filter(uf == 15) %>% 
  droplevels()

# Baseline, other states
df_mahogany_baseline_not_para <- df_mahogany_baseline %>% 
  filter(uf != 15) %>% 
  droplevels()

# Pará, years in which there' detailed homicide data
df_mahogany_para_demog <- df_mahogany_para %>% 
  filter(year %in% 1996:2013)

## Pará, years in which there's data on land land conflicts
df_mahogany_para_pol <- df_mahogany_para %>% 
  filter(year %in% 1995:2007)


# Replicating paper's main diff-in-diff specifications ------------------------
# Note: the models' names indicate the corresponding table and column
# E.g.: t2c3 corresponds to the model whose results are shown in the 3rd
# column of the 2nd table from the paper

# Table 2: benchmark specification --------------------------------------------
t2c1 <- 
  felm(formula = hom_tx ~ treat1 + treat2 + treat3 + year | code | 0 | code,
       data = df_mahogany, weights = df_mahogany$avg_pop)

t2c2 <- 
  felm(formula = hom_tx ~ treat1 + treat2 + treat3 + uf*year + year | code
       | 0 | code,
       data = df_mahogany, weights = df_mahogany$avg_pop)

t2c3 <- 
  felm(formula = hom_tx ~ treat1 + treat2 + treat3 + uf*year +
         year*hom_base + year*area_base + year*lngdp_base + year*gdpag_base +
         year*pol_base + year*und5_base + year*infecc_base + year*heart_base +
         year*neop_base + year*suic_base + year*traff_base + year | code
       | 0 | code,
       data = df_mahogany_baseline, weights = df_mahogany_baseline$avg_pop)

t2c4 <- 
  felm(formula = hom_tx ~ treat1 + int1 + treat2 + int2 + 
         treat3 + int3 + uf*year + year | code | 0 | code,
         data = df_mahogany, weights= df_mahogany$avg_pop)

t2c5 <- 
  felm(formula = hom_tx ~ treatx_area1 + treatx_area2 + treatx_area3 +
         year + uf*year | code | 0 | code,
       data = df_mahogany, weights = df_mahogany$avg_pop)

t2c6 <- 
  felm(formula = hom_tx ~ treatt_area1 + treatt_area2 + treatt_area3 +
         year + uf*year | code | 0 | code,
       data = df_mahogany, weights = df_mahogany$avg_pop)

# Showing results
t2_results <- stargazer(t2c1, t2c2, t2c3, t2c4, t2c5, t2c6,
                        omit = c("uf", "year", "_base"),
                        type = "text")

# Table 3: Pará vs. other states ----------------------------------------------
t3c1 <- 
  felm(formula = hom_tx ~ treat1 + treat2 + treat3 + year | code | 0 | code,
       data = df_mahogany_para, weights = df_mahogany_para$avg_pop)

t3c2 <- 
  felm(formula = hom_tx ~ treat1 + treat2 + treat3 +
         year*hom_base + year*area_base + year*lngdp_base + year*gdpag_base +
         year*pol_base + year*und5_base + year*infecc_base + year*heart_base +
         year*neop_base + year*suic_base + year*traff_base + year 
       | code | 0 | code,
       data = df_mahogany_baseline_para,
       weights = df_mahogany_baseline_para$avg_pop)

t3c3 <- 
  felm(formula = hom_tx ~ treat1 + int1 + treat2 + int2 + treat3 + int3 + year 
       | code | 0 | code, 
       data = df_mahogany_para, weights = df_mahogany_para$avg_pop)

t3c4 <- 
  felm(formula = hom_tx ~ treat1 + treat2 + treat3 + year | code | 0 | code,
       data = df_mahogany_not_para, weights = df_mahogany_not_para$avg_pop)

t3c5 <- 
  felm(formula = hom_tx ~ treat1 + treat2 + treat3 + year + uf*year 
       | code | 0 | code,
       data = df_mahogany_not_para, weights = df_mahogany_not_para$avg_pop)

t3c6 <- 
  felm(formula = hom_tx ~ treat1 + treat2 + treat3 + uf*year +
         year*hom_base + year*area_base + year*lngdp_base + year*gdpag_base +
         year*pol_base + year*und5_base + year*infecc_base + year*heart_base +
         year*neop_base + year*suic_base + year*traff_base + year 
       | code | 0 | code,
       data = df_mahogany_not_para, weights = df_mahogany_not_para$avg_pop)

t3c7 <- 
  felm(formula = hom_tx ~ treat1 + int1 + treat2 + int2 + treat3 + int3 +
         year + uf*year | code | 0 | code,
       data = df_mahogany_not_para, weights = df_mahogany_not_para$avg_pop)

# Showing results
t3_results <- stargazer(t3c1, t3c2, t3c3, t3c4, t3c5, t3c6, t3c7,
                        type = "text", omit = c("uf", "year", "_base"))

# Table 4: Placebo, municipality-specific trends, simultaneous changes --------
t4c1 <- 
  felm(formula = hom_tx ~ treat1 + treat2 + treat3 +pre + year 
       | code | 0 | code,
       data = df_mahogany_para, weights = df_mahogany_para$avg_pop)

t4c2 <- 
  felm(formula = hom_tx ~ treat1 + treat2 + treat3 + code*trend + year
       | code | 0 | code,
       data = df_mahogany_para, weights = df_mahogany_para$avg_pop)

t4c3 <- 
  felm(formula = hom_tx ~ treat1 + treat2 + treat3 + pre + year 
       | code | 0 | code,
       data = df_mahogany_para_gdp, weights = df_mahogany_para_gdp$avg_pop)

t4c4 <- 
  felm(formula = hom_tx ~ treat1 + treat2 + treat3 + code*trend + year
       | code | 0 | code,
       data = df_mahogany_para_gdp, weights = df_mahogany_para_gdp$avg_pop)

t4c5 <- 
  felm(formula = ln_gdp_pc ~ treat1 + treat2 + treat3 + pre + year 
       | code | 0 | code,
       data = df_mahogany_para_gdp, weights = df_mahogany_para_gdp$avg_pop)

t4c6 <- 
  felm(formula = ln_gdp_pc ~ treat1 + treat2 + treat3 + code*trend + year 
       | code | 0 | code,
       data = df_mahogany_para_gdp, weights = df_mahogany_para_gdp$avg_pop)

t4c7 <- 
  felm(formula = gdpfr_ag ~ treat1 + treat2 + treat3 + pre + year 
       | code | 0 | code,
       data = df_mahogany_para_gdp, weights = df_mahogany_para_gdp$avg_pop)

t4c8 <- 
  felm(formula = gdpfr_ag ~ treat1 + treat2 + treat3 + code*trend + year 
       | code | 0 | code,
       data = df_mahogany_para_gdp,
       weights = df_mahogany_para_gdp$avg_pop)

# Showing results
t4_results <- stargazer(t4c1, t4c2, t4c3, t4c4, t4c5, t4c6, t4c7, t4c8,
                        omit = c("uf", "year", "_base", "code", "trend"),
                        type = "text")

# Table 5: demographic characterization of homicides --------------------------
t5c1 <- 
  felm(formula = hom_tx ~ treat1 + treat2 + treat3 + year 
       | code | 0 | code,
       data = df_mahogany_para_demog, weights = df_mahogany_para_demog$avg_pop)

t5c2 <- 
  felm(formula = male_h ~ treat1 + treat2 + treat3 + year 
       | code | 0 | code,
       data = df_mahogany_para_demog, weights = df_mahogany_para_demog$avg_pop)

t5c3 <- 
  felm(formula = male_h_prime ~ treat1 + treat2 + treat3 + year 
       | code | 0 | code,
       data = df_mahogany_para_demog, weights = df_mahogany_para_demog$avg_pop)

t5c4 <- 
  felm(formula = male_h_sing ~ treat1 + treat2 + treat3 + year 
       | code | 0 | code,
       data = df_mahogany_para_demog, weights = df_mahogany_para_demog$avg_pop)

t5c5 <- 
  felm(formula = male_h_nothome ~ treat1 + treat2 + treat3 + year 
       | code | 0 | code,
       data = df_mahogany_para_demog, weights = df_mahogany_para_demog$avg_pop)

t5c6 <- 
  felm(formula = male_h_firearm ~ treat1 + treat2 + treat3 + year
       | code | 0 | code,
                      data = df_mahogany_para_demog,
                      weights = df_mahogany_para_demog$avg_pop)

t5c7 <- 
  felm(formula = pol_m ~ treat1 + treat2 + treat3 + year 
       | code | 0 | code,
       data = df_mahogany_para_pol, weights = df_mahogany_para_pol$avg_pop)

t5c8 <- 
  felm(formula = pol_deaths_dummy ~ treat1 + treat2 + treat3 + year 
       | code | 0 | code,
       data = df_mahogany_para_pol, weights = df_mahogany_para_pol$avg_pop)

# Showing results
t5_results <- stargazer(t5c1, t5c2, t5c3, t5c4, t5c5, t5c6, t5c7, t5c8,
                        type = "text", omit = c("uf", "year", "_base"))
