library(tidyverse); library(TOSTER)

# Load data
#dem_all <- read_csv(here::here("data", "clean", "ospan_set_z_scores.csv"))
dem_all <- read_csv(here::here("data", "pupurri_analysis.csv"))


# Remove participants to make the groups homogeneous in L2 proficiency (DELE)

unique(dem_all$participant)

dem_all <- dem_all %>%
  # separate(., col = participant,
  #          into = c('group', 'id'),
  #          sep = 3,
  #          remove = FALSE) %>%
  separate(., col = group,
           into = c("proficiency", "l1"), # es = EN speaker, ms = MA speaker, on = ES speaker
           sep = 1,
           remove = FALSE)

dem_all$l1 <- str_replace(dem_all$l1, "es", "en")
dem_all$l1 <- str_replace(dem_all$l1, "ms", "ma")
dem_all$l1 <- str_replace(dem_all$l1, "on", "es")


dem_all %>%
  group_by(l1) %>%
  summarize(n = n()) # en = Eglish speakers, ma = Mandarin Chinese speakers, es = Spanish speakers
#   l1        n
# 1 en       65
# 2 es       30
# 3 ma       64


dem_all %>%
  group_by(., l1) %>%
  summarise(., ospan_mean = round(mean(WM_set),2),   # ospan instead of WM_set if using z scores
            ospan_sd = round(sd(WM_set),2),          # ospan instead of WM_set if using z scores
            n = length(unique(participant))) %>%
  knitr::kable()
#   |l1 | ospan_mean| ospan_sd|  n|
#   |:--|----------:|--------:|--:|
#   |en |          0|     2.11| 65|
#   |es |          0|     2.72| 30|
#   |ma |          0|     2.14| 64|

#   |l1 | ospan_mean| ospan_sd|  n|
#   |:--|----------:|--------:|--:|
#   |en |       8.89|     2.11| 65|
#   |es |       6.20|     2.72| 30|
#   |ma |       7.78|     2.14| 64|



zsc <- dem_all %>% 
  group_by(., l1) %>%
  mutate(zscoreWM = (WM_set - mean(WM_set))/sd(WM_set)) 




bartlett.test(zscoreWM ~ l1, data = zsc)
# Bartlett's K-squared = 6.376e-15, df = 2, p-value = 1




### with sd

TOSTtwo(m1 = 0, sd1 = 2.11, n1 = 65, # EN
        m2 = 0, sd2 = 2.14, n2 = 64, # MA
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)
# The null hypothesis test was non-significant, t(126.89) = 0.000, p = 1.000, given an alpha of 0.05.

TOSTtwo(m1 = 0, sd1 = 2.11, n1 = 65, # EN
        m2 = 0, sd2 = 2.72, n2 = 30, # ES
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)
# The null hypothesis test was non-significant, t(45.75) = 0.000, p = 1.000, given an alpha of 0.05.

TOSTtwo(m1 = 0, sd1 = 2.14, n1 = 64, # MA
        m2 = 0, sd2 = 2.72, n2 = 30, # ES
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)
# The null hypothesis test was non-significant, t(46.47) = 0.000, p = 1.000, given an alpha of 0.05.











# Check info: id, group, dele, pstm, wm, age, aoa_l2, months_abroad,
# l1_use_week, l2_use_week, years_word_int

# Set DELE as numeric

dem_all$DELE <- as.numeric(dem_all$DELE)


#### BY L1
# Create a table with mean + sd for: wm, pstm, dele, aoa_l2, months_abroad,
# l1_use_week, l2_use_week, years_word_int

dem_all %>%
  group_by(., l1) %>%
  filter(., l1 != "es") %>%
  summarise(., dele = round(mean(DELE),2),
            dele_sd = round(sd(DELE),2),
            n = length(unique(participant))) %>%
  knitr::kable()

# |l1 |  dele| dele_sd|  n|
# |:--|-----:|-------:|--:|
# |en | 38.48|    8.19| 65|
# |ma | 39.17|    7.56| 64|


dem_all %>%
  group_by(., l1) %>%
  summarise(aoa_l2 = round(mean(AoA_L2),2),
            aoa_l2_sd = round(sd(AoA_L2),2),
            abroad = round(mean(mo_ES_country),2),
            abroad_sd = round(sd(mo_ES_country),2),
            l1_use = round(mean(percent_l1_week),2),
            l1_use_sd = round(sd(percent_l1_week),2),
            l2_use = round(mean(percent_l2_week),2),
            l2_use_sd = round(sd(percent_l2_week),2),
            n = length(unique(participant))
  ) %>% knitr::kable()
# |l1 | aoa_l2| aoa_l2_sd| abroad| abroad_sd| l1_use| l1_use_sd| l2_use| l2_use_sd|  n|
# |:--|------:|---------:|------:|---------:|------:|---------:|------:|---------:|--:|
# |en |  16.29|      5.55|  38.08|     33.48|  66.62|     17.44|  33.31|     17.44| 65|
# |es |   7.03|      6.82|   0.30|      1.64|  10.73|      9.28|  88.93|      9.43| 30|
# |ma |  18.88|      3.59|  40.83|     45.46|  58.05|     22.14|  41.64|     21.66| 64|


## Homogeneity of variances tests
dem_all <- dem_all %>%
  filter(., l1 != "es")

bartlett.test(DELE ~ l1, data = dem_all) # K-squared = 0.39147, df = 1, p-value = 0.5315
bartlett.test(AoA_L2 ~ l1, data = dem_all) # K-squared = 11.629, df = 1, p-value = 0.0006495
bartlett.test(mo_ES_country ~ l1, data = dem_all) # K-squared = 5.8093, df = 1, p-value = 0.01594
bartlett.test(percent_l1_week ~ l1, data = dem_all) # K-squared = 3.5658, df = 1, p-value = 0.05898
bartlett.test(percent_l2_week ~ l1, data = dem_all) # K-squared = 2.9384, df = 1, p-value = 0.0865


# time abroad 
t.test(mo_ES_country ~ l1, data = dem_all, var.equal = TRUE)

TOSTtwo(m1 = 38.08, sd1 = 33.48, n1 = 65, # en
        m2 = 40.83, sd2 = 45.46, n2 = 64, # ma
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)
# The null hypothesis test was non-significant, t(115.76) = -0.391, p = 0.697, given an alpha of 0.05.


# dele 
t.test(DELE ~ l1, data = dem_all, var.equal = TRUE)
# t = -0.50059, df = 127, p-value = 0.6175

TOSTtwo(m1 = 38.48, sd1 = 8.19, n1 = 65, # EN
        m2 = 39.17, sd2 = 7.56, n2 = 64, # MA
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)
# t(126.48) = -0.497, p = 0.620


# L2 weekly use
# significant and outside area between dotted lines
TOSTtwo(m1 = 33.31, sd1 = 17.44, n1 = 65, # EN
        m2 = 41.64, sd2 = 21.66, n2 = 64, # MA
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)
# t(120.69) = -2.404, p = 0.0178

t.test(percent_l2_week ~ l1, data = dem_all, var.equal = TRUE)
# t = -2.4087, df = 127, p-value = 0.01745



arrange(dem_all, (desc(percent_l2_week))) %>% View()

t.test(percent_l2_week ~ l1, data = dem_all %>% filter(participant != 'IES04' & participant != 'IES17' & 
                                                         participant != 'AES32' & participant != 'IES28'), var.equal = TRUE)
# t = -1.9516, df = 123, p-value = 0.05326
bartlett.test(percent_l2_week ~ l1, data = dem_all %>% filter(participant != 'IES04' & participant != 'IES17' & 
                                                                participant != 'AES32' & participant != 'IES28'))
# K-squared = 3.6896, df = 1, p-value = 0.05475
t.test(DELE ~ l1, data = dem_all %>% filter(participant != 'IES04' & participant != 'IES17' & 
                                                         participant != 'AES32' & participant != 'IES28'), var.equal = TRUE)
# t = -0.14658, df = 123, p-value = 0.8837
bartlett.test(DELE ~ l1, data = dem_all %>% filter(participant != 'IES04' & participant != 'IES17' & 
                                                                participant != 'AES32' & participant != 'IES28'))

# Bartlett's K-squared = 0.23176, df = 1, p-value = 0.6302

dem_all %>%
  filter(., participant != 'IES04' & participant != 'IES17' & 
           participant != 'AES32' & participant != 'IES28') %>%
  group_by(., l1) %>%
  summarise(
            l2_use = round(mean(percent_l2_week),2),
            l2_use_sd = round(sd(percent_l2_week),2),
            n = length(unique(participant))
  ) %>% knitr::kable()
# |l1 | l2_use| l2_use_sd|  n|
#   |:--|------:|---------:|--:|
#   |en |  34.84|     16.91| 61|
#   |ma |  41.64|     21.66| 64|




## BY GROUP

# Create a table with mean + sd for: wm, pstm, dele, aoa_l2, months_abroad,
# l1_use_week, l2_use_week, years_word_int

dem_all %>%
  group_by(., group) %>%
  filter(., group != "mon") %>%
  summarise(., dele = round(mean(DELE),2),
               dele_sd = round(sd(DELE),2),
            n = length(unique(participant))) %>%
  knitr::kable()

# |group |  dele| dele_sd|  n|
# |:-----|-----:|-------:|--:|
# |aes   | 45.44|    4.26| 32|
# |ams   | 45.50|    3.97| 32|
# |ies   | 31.73|    4.58| 33|
# |ims   | 32.84|    4.23| 32|




dem_all %>%
  group_by(., group) %>%
  summarise(aoa_l2 = round(mean(AoA_L2),2),
            aoa_l2_sd = round(sd(AoA_L2),2),
            abroad = round(mean(mo_ES_country),2),
            abroad_sd = round(sd(mo_ES_country),2),
            l1_use = round(mean(percent_l1_week),2),
            l1_use_sd = round(sd(percent_l1_week),2),
            l2_use = round(mean(percent_l2_week),2),
            l2_use_sd = round(sd(percent_l2_week),2),
            n = length(unique(participant))
            ) %>% knitr::kable()
# |group | aoa_l2| aoa_l2_sd| abroad| abroad_sd| l1_use| l1_use_sd| l2_use| l2_use_sd|  n|
# |:-----|------:|---------:|------:|---------:|------:|---------:|------:|---------:|--:|
# |aes   |  15.06|      4.35|  51.28|     38.83|  61.25|     16.16|  38.59|     16.23| 32|
# |ams   |  17.88|      2.83|  45.12|     55.68|  52.81|     22.57|  46.56|     21.76| 32|
# |ies   |  17.48|      6.35|  25.29|     20.94|  71.82|     17.27|  28.18|     17.27| 33|
# |ims   |  19.88|      4.01|  36.53|     32.58|  63.28|     20.74|  36.72|     20.74| 32|
# |mon   |   7.03|      6.82|   0.30|      1.64|  10.73|      9.28|  88.93|      9.43| 30|



## Homogeneity of variances tests
dem_all <- dem_all %>%
  filter(., group != "mon")
## All seem ok
bartlett.test(DELE ~ l1, data = dem_all)
bartlett.test(AoA_L2 ~ l1, data = dem_all) # p-value = 0.0006495
bartlett.test(mo_ES_country ~ l1, data = dem_all) # p-value = 0.01594
bartlett.test(percent_l1_week ~ l1, data = dem_all)
bartlett.test(percent_l2_week ~ l1, data = dem_all)





# dele tost intermediate groups
# all good
TOSTtwo(m1 = 31.73, sd1 = 4.58, n1 = 33, # EN
        m2 = 32.84, sd2 = 4.23, n2 = 32, # MA
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)

# dele tost advanced groups
# all good
TOSTtwo(m1 = 45.44, sd1 = 4.26, n1 = 32, # EN
        m2 = 45.50, sd2 = 3.97, n2 = 32, # MA
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)




# age of acquistion tost intermediate
# not significant but outside the area between dotted lines 
TOSTtwo(m1 = 17.48, sd1 = 6.35, n1 = 33, # EN
        m2 = 19.88, sd2 = 4.01, n2 = 32, # MA
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)

# age of acquistion tost advanced
# signicant and outside area between dotted lines
TOSTtwo(m1 = 15.06, sd1 = 4.35, n1 = 32, # EN
        m2 = 17.88, sd2 = 2.83, n2 = 32, # MA
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)





# time abroad tost intermediate
# not significant but outside the area between dotted lines 
TOSTtwo(m1 = 25.29, sd1 = 20.94, n1 = 33, # en
        m2 = 36.53, sd2 = 32.58, n2 = 32, # ma
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)

# time abroad tost advanced
# good
TOSTtwo(m1 = 51.28, sd1 = 38.83, n1 = 32, # en
        m2 = 45.12, sd2 = 55.68, n2 = 32, # ma
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)

# |group | | l1_use| l1_use_sd| l2_use| l2_use_sd|  n|
# |:-----|-|------:|---------:|------:|---------:|--:|
# |aes   | |  61.25|     16.16|  38.59|     16.23| 32|
# |ams   | |  52.81|     22.57|  46.56|     21.76| 32|
# |ies   | |  71.82|     17.27|  28.18|     17.27| 33|
# |ims   | |  63.28|     20.74|  36.72|     20.74| 32|

# l2 use in a normal week tost intermediate 
# not significant but outside the area between dotted lines 
TOSTtwo(m1 = 28.18, sd1 = 17.27, n1 = 33, # en
        m2 = 36.72, sd2 = 20.74, n2 = 32, # ma
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)

# l2 use in a normal week tost advanced 
# not significant but outside the area between dotted lines 
TOSTtwo(m1 = 38.59, sd1 = 16.23, n1 = 32, # en
        m2 = 46.56, sd2 = 21.76, n2 = 32, # ma
        low_eqbound_d = -0.3,
        high_eqbound_d = 0.3,
        alpha = 0.05)
