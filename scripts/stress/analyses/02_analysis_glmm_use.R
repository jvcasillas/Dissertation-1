#
# Original script by Joseph
# Modified by Laura
# Last update: 2021/03/09
#
# GLMMs -----------------------------------------------------------------------
#
# - Target fixation as a function of group, stress (condition), and coda
#   at the offset of the first syllable (time_zero == 20)
# - This model builds on the t-test analyses by looking for between group
#   differences in target fixation at the offset of first syllable
#
# -----------------------------------------------------------------------------



# Onset c1

# Load data and models --------------------------------------------------------

source(here::here("scripts", "00_load_libs.R"))
# source(here::here("scripts", "02_load_data.R"))
stress10 <- read_csv(here::here("data", 'clean', 'stress_10ms_final_onset_c1.csv'))

mon_onset_c1_final   <- readRDS(here("mods", 'stress', "glmm", "onset_c1", "mon_onset_c1_final.rds"))
prop_0_mon_mod_cond  <- readRDS(here("mods", 'stress', "glmm", "onset_c1", "prop_0_mon_mod_cond.rds"))
prop_0_mon_mod_wm    <- readRDS(here("mods", 'stress', "glmm", "onset_c1", "prop_0_mon_mod_wm.rds"))
prop_0_mon_mod_int   <- readRDS(here("mods", 'stress', "glmm", "onset_c1", "prop_0_mon_mod_int.rds"))


# -----------------------------------------------------------------------------



# Data prep -------------------------------------------------------------------

# Filter time course to offset of 1st syllable (time_zero == 20, to account for 200 ms to launch saccade)
# Create sum coded fixed factors (condition)


df_stress <- stress10 %>%
  filter(., time_zero == 20) %>%
  mutate(condition_sum = if_else(cond == 1, "Present", "Preterit"),
         condition_sum = fct_relevel(condition_sum, "Present"))

stress_mon <- df_stress %>%
  filter(., l1 == 'es')



# -----------------------------------------------------------------------------




##########################     MONOLINGUALS    ##########################

# Random effects building -----------------------------------------------------

if(F) {
  
  prop_0_ranefA <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                           (1 | participant),
                         data = stress_mon, family = 'binomial',                 
                         control = glmerControl(optimizer = 'bobyqa'))
  
  prop_0_ranefB <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                           (1 | participant) +
                           (1 | target),
                         data = stress_mon, family = 'binomial', 
                         control = glmerControl(optimizer = 'bobyqa'))
  
  anova(prop_0_ranefA, prop_0_ranefB, refit = F) # keep intercept for target bc significant
  #                 Df    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
  # prop_0_ranefA    2 6377.4 6385.8 -3186.7   6373.4                     
  # prop_0_ranefB    3 5908.8 5921.3 -2951.4   5902.8 470.67  1  < 2.2e-16 ***
  
  prop_0_ranefC <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                           (1 + condition_sum | participant) +
                           (1 | target),
                         data = stress_mon, family = 'binomial', 
                         control = glmerControl(optimizer = 'bobyqa'))
  
  anova(prop_0_ranefB, prop_0_ranefC, refit = F) # keep slope for condition bc significant
  #                 Df    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
  # prop_0_ranefB    3 5908.8 5921.3 -2951.4   5902.8                     
  # prop_0_ranefC    5 5768.2 5789.1 -2879.1   5758.2 144.56  2  < 2.2e-16 ***
  
}

# -----------------------------------------------------------------------------





# Test fixed effects ----------------------------------------------------------

if(F) {
  prop_0_mon_mod_0 <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                              (1 + condition_sum | participant) +   
                              (1 | target),
                            data = stress_mon, family = 'binomial', 
                            control = glmerControl(optimizer = 'bobyqa'))
  
  prop_0_mon_mod_cond  <- update(prop_0_mon_mod_0, . ~ . + condition_sum)
  prop_0_mon_mod_wm    <- update(prop_0_mon_mod_cond,  . ~ . + ospan)
  prop_0_mon_mod_int  <- update(prop_0_mon_mod_wm,    . ~ . + condition_sum:ospan)
  
  anova(prop_0_mon_mod_0, prop_0_mon_mod_cond, test = "Chisq")  
  #                     npar    AIC    BIC  logLik deviance  Chisq  Df Pr(>Chisq)    
  # prop_0_mon_mod_0       5 5768.2 5789.1 -2879.1   5758.2                     
  # prop_0_mon_mod_cond    6 5769.1 5794.2 -2878.6   5757.1 1.0593  1     0.3034
  
  anova(prop_0_mon_mod_0, prop_0_mon_mod_wm, test = "Chisq") 
  #                   npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
  # prop_0_mon_mod_0     5 5768.2 5789.1 -2879.1   5758.2                     
  # prop_0_mon_mod_wm    7 5770.5 5799.7 -2878.3   5756.5 1.6891  2     0.4298
  
  anova(prop_0_mon_mod_0, prop_0_mon_mod_int, test = "Chisq") 
  #                    npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
  # prop_0_mon_mod_0      5768.2 5789.1 -2879.1   5758.2                     
  # prop_0_mon_mod_int    8 5772.5 5805.9 -2878.3   5756.5 1.6892  3     0.6393
  
  mon_onset_c1_final <- prop_0_mon_mod_0
  
  summary(mon_onset_c1_final) 
  # Fixed effects:
  #                       Estimate Std. Error z value Pr(>|z|)    
  # (Intercept)            -0.1378     0.1894  -0.727    0.467
  
}


# Save models

saveRDS(mon_onset_c1_final, here("mods", "stress", 
                               "glmm", "onset_c1", "mon_onset_c1_final.rds"))
saveRDS(prop_0_mon_mod_cond, here("mods", "stress", 
                                  "glmm", "onset_c1", "prop_0_mon_mod_cond.rds"))
saveRDS(prop_0_mon_mod_wm, here("mods", "stress", 
                                "glmm", "onset_c1", "prop_0_mon_mod_wm.rds"))
saveRDS(prop_0_mon_mod_int, here("mods", "stress", 
                                 "glmm", "onset_c1", "prop_0_mon_mod_int.rds"))

mon_onset_c1_final %>%
  ggplot(., aes(x = condition_sum, y = mean(target_prop))) + 
  geom_hline(yintercept = 0.5, color = "black", size = 0.75,
             lty = 3) + 
  stat_summary(fun.data = mean_cl_boot, geom = 'errorbar', 
               position = position_dodge(width = 0.5), 
               width = 0.35, color = 'grey40') + 
  stat_summary(fun.y = mean, geom = 'point', size = 4,
               position = position_dodge(width = 0.5)) + 
  coord_cartesian(ylim = c(0, 1)) + ylab('Target fixations') + xlab('Stress condition') + 
  ggtitle('Mean target fixations as a function of group\nand target type') +    
  theme_bw(base_size = 16, base_family = 'Times') -> stress_target_fix_mon_onset_c1

ggsave(paste0("./figs/stress/glmm/mon/stress_target_fix_mon_onset_c1.png"), stress_target_fix_mon_onset_c1, width = 150,
       height = 120, units = "mm", dpi = 600)

# -----------------------------------------------------------------------------




# Onset v1

# Load data and models --------------------------------------------------------

source(here::here("scripts", "00_load_libs.R"))
# source(here::here("scripts", "02_load_data.R"))
stress10 <- read_csv(here::here("data", 'clean', 'stress_10ms_final_onset_v1.csv'))

mon_onset_v1_final     <- readRDS(here("mods", 'stress', "glmm", "onset_v1", "mon_onset_v1_final.rds"))
prop_0_mon_mod_cond  <- readRDS(here("mods", 'stress', "glmm", "onset_v1", "prop_0_mon_mod_cond.rds"))
prop_0_mon_mod_wm    <- readRDS(here("mods", 'stress', "glmm", "onset_v1", "prop_0_mon_mod_wm.rds"))
prop_0_mon_mod_int   <- readRDS(here("mods", 'stress', "glmm", "onset_v1", "prop_0_mon_mod_int.rds"))

l2_onset_v1_delewm_final <- readRDS(here("mods", 'stress', "glmm", "onset_v1", "l2_onset_v1_delewm_final.rds"))
prop_0_l2_mod_l1 <- readRDS(here("mods", 'stress', "glmm", "onset_v1", "prop_0_l2_mod_l1.rds"))
prop_0_l2_mod_cond <- readRDS(here("mods", 'stress', "glmm", "onset_v1", "prop_0_l2_mod_cond.rds"))
prop_0_l2_mod_wm <- readRDS(here("mods", 'stress', "glmm", "onset_v1", "prop_0_l2_mod_wm.rds"))
prop_0_l2_mod_dele <- readRDS(here("mods", 'stress', "glmm", "onset_v1", "prop_0_l2_mod_dele.rds"))
l2_onset_v1_use_final <- readRDS(here("mods", 'stress', "glmm", "onset_v1", "l2_onset_v1_use_final.rds"))
prop_0_l2_mod_int  <- readRDS(here("mods", 'stress', "glmm", "onset_v1", "prop_0_l2_mod_int.rds"))
prop_0_l2_mod_int_wm  <- readRDS(here("mods", 'stress', "glmm", "onset_v1", "prop_0_l2_mod_int_wm.rds"))
prop_0_l2_mod_int_dele  <- readRDS(here("mods", 'stress', "glmm", "onset_v1", "prop_0_l2_mod_int_dele.rds"))
prop_0_l2_mod_int_use  <- readRDS(here("mods", 'stress', "glmm", "onset_v1", "prop_0_l2_mod_int_use.rds"))

# -----------------------------------------------------------------------------



# Data prep -------------------------------------------------------------------

# Filter time course to offset of 1st syllable (time_zero == 20, to account for 200 ms to launch saccade)
# Create sum coded fixed factors (condition)


df_stress <- stress10 %>%
  filter(., time_zero == 20) %>%
  mutate(condition_sum = if_else(cond == 1, 1, -1))

stress_mon <- df_stress %>%
  filter(., l1 == 'es')

stress_l2 <- df_stress %>%
  filter(., l1 != 'es') %>%
  mutate(., l1 = fct_relevel(l1, "en", "ma"))
        


# -----------------------------------------------------------------------------




##########################     MONOLINGUALS    ##########################

# Random effects building -----------------------------------------------------

if(F) {

prop_0_ranefA <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                    (1 | participant),
                    data = stress_mon, family = 'binomial',                 
                    control = glmerControl(optimizer = 'bobyqa'))

prop_0_ranefB <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                    (1 | participant) +
                    (1 | target),
                    data = stress_mon, family = 'binomial', 
                    control = glmerControl(optimizer = 'bobyqa'))

anova(prop_0_ranefA, prop_0_ranefB, refit = F) # keep intercept for target bc significant
#                 Df    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
# prop_0_ranefA    2 6262.7 6271.1 -3129.4   6258.7                         
# prop_0_ranefB    3 5804.4 5816.9 -2899.2   5798.4 460.38  1  < 2.2e-16 ***

prop_0_ranefC <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                    (1 + condition_sum | participant) +
                    (1 | target),
                    data = stress_mon, family = 'binomial', 
                    control = glmerControl(optimizer = 'bobyqa'))

anova(prop_0_ranefB, prop_0_ranefC, refit = F) # keep slope for condition bc significant
#                 Df    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
# prop_0_ranefB    3 5804.4 5816.9 -2899.2   5798.4                         
# prop_0_ranefC    5 5680.4 5701.2 -2835.2   5670.4 127.96  2  < 2.2e-16 ***

 }

# -----------------------------------------------------------------------------





# Test fixed effects ----------------------------------------------------------

if(F) {
prop_0_mon_mod_0 <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                    (1 + condition_sum | participant) +   
                    (1 | target),
                    data = stress_mon, family = 'binomial', 
                    control = glmerControl(optimizer = 'bobyqa'))

prop_0_mon_mod_cond  <- update(prop_0_mon_mod_0, . ~ . + condition_sum)
prop_0_mon_mod_wm    <- update(prop_0_mon_mod_cond,  . ~ . + ospan)
prop_0_mon_mod_int  <- update(prop_0_mon_mod_wm,    . ~ . + condition_sum:ospan)

anova(prop_0_mon_mod_0, prop_0_mon_mod_cond, test = "Chisq")   # no effect of condition bc not significant
#                     npar    AIC    BIC  logLik deviance  Chisq  Df Pr(>Chisq)    
# prop_0_mon_mod_0       5 5680.4 5701.2 -2835.2   5670.4                     
# prop_0_mon_mod_cond    6 5681.0 5706.1 -2834.5   5669.0 1.3438  1     0.2464

anova(prop_0_mon_mod_0, prop_0_mon_mod_wm, test = "Chisq") # no effect of wm bc not significant
#                   npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
# prop_0_mon_mod_0     5 5680.4 5701.2 -2835.2   5670.4                     
# prop_0_mon_mod_wm    7 5681.1 5710.3 -2833.5   5667.1 3.3389  2     0.1883

anova(prop_0_mon_mod_0, prop_0_mon_mod_int, test = "Chisq") # no interaction wm x condition
#                    npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
# prop_0_mon_mod_0      5 5680.4 5701.2 -2835.2   5670.4                     
# prop_0_mon_mod_int    8 5682.0 5715.4 -2833.0   5666.0 4.3618  3     0.2249

mon_onset_v1_final <- prop_0_mon_mod_0
  
summary(mon_onset_v1_final) 
# Fixed effects:
#                       Estimate Std. Error z value Pr(>|z|)    
# (Intercept)             0.4016     0.1865   2.154   0.0313 *
}
# -----------------------------------------------------------------------------






##########################     L2 SPEAKERS    ##########################

# Random effects building -----------------------------------------------------

if(F) {
  
  prop_0_ranefA <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                           (1 | participant),
                         data = stress_l2, family = 'binomial',                 
                         control = glmerControl(optimizer = 'bobyqa'))
  
  prop_0_ranefB <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                           (1 | participant) +
                           (1 | target),
                         data = stress_l2, family = 'binomial', 
                         control = glmerControl(optimizer = 'bobyqa'))
  
  anova(prop_0_ranefA, prop_0_ranefB, refit = F) # keep intercept for target bc significant
  #               npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
  # prop_0_ranefA    2 26951 26962 -13474    26947                         
  # prop_0_ranefB    3 26490 26507 -13242    26484 462.66  1  < 2.2e-16 ***
  
  prop_0_ranefC <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                           (1 + condition_sum | participant) +
                           (1 | target),
                         data = stress_l2, family = 'binomial', 
                         control = glmerControl(optimizer = 'bobyqa'))
  
  anova(prop_0_ranefB, prop_0_ranefC, refit = F) # keep slope for condition bc significant
  #               npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
  # prop_0_ranefB    3 26490 26507 -13242    26484                         
  # prop_0_ranefC    5 25570 25599 -12780    25560 923.87  2  < 2.2e-16 ***
  
}

# -----------------------------------------------------------------------------





# Test fixed effects ----------------------------------------------------------

if(F) {
  prop_0_l2_mod_0 <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                              (1 + condition_sum | participant) +   
                              (1 | target),
                            data = stress_l2, family = 'binomial', 
                            control = glmerControl(optimizer = 'bobyqa'))
  
  prop_0_l2_mod_l1    <- update(prop_0_l2_mod_0,     . ~ . + l1)
  prop_0_l2_mod_cond  <- update(prop_0_l2_mod_l1,    . ~ . + condition_sum)
  prop_0_l2_mod_int   <- update(prop_0_l2_mod_cond,  . ~ . + l1:condition_sum)
  
  anova(prop_0_l2_mod_0, prop_0_l2_mod_l1, test = "Chisq")
  #                  npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # prop_0_l2_mod_0     5 25570 25599 -12780    25560                     
  # prop_0_l2_mod_l1    6 25572 25606 -12780    25560 0.1414  1     0.7069
  
  anova(prop_0_l2_mod_0, prop_0_l2_mod_cond, test = "Chisq")   
  #                    npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # prop_0_l2_mod_0       5 25570 25599 -12780    25560                     
  # prop_0_l2_mod_cond    7 25574 25613 -12780    25560 0.7191  2      0.698
  
  anova(prop_0_l2_mod_0, prop_0_l2_mod_int, test = "Chisq") 
  # npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # prop_0_l2_mod_0      5 25570 25599 -12780    25560                     
  # prop_0_l2_mod_int    8 25576 25621 -12780    25560 0.8291  3     0.8425
  
  # ROUTE 1
  
  prop_0_l2_mod_wm       <- update(prop_0_l2_mod_0,  . ~ . + ospan)
  prop_0_l2_mod_int_wm   <- update(prop_0_l2_mod_wm,    . ~ . + l1:condition_sum:ospan)
  
  anova(prop_0_l2_mod_0, prop_0_l2_mod_wm, test = "Chisq") 
  # npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # prop_0_l2_mod_0     5 25570 25599 -12780    25560                     
  # prop_0_l2_mod_wm    6 25572 25606 -12780    25560 0.0438  1     0.8343
  
  anova(prop_0_l2_mod_0, prop_0_l2_mod_int_wm, test = "Chisq") # no interaction wm x condition
  #                      npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # prop_0_l2_mod_0         5 25570 25599 -12780    25560                    
  # prop_0_l2_mod_int_wm    8 25574 25619 -12779    25558 2.1518  3     0.5415
  
  
  # ROUTE 2
  
  prop_0_l2_mod_dele       <- update(prop_0_l2_mod_0,  . ~ . + DELE)
  prop_0_l2_mod_int_dele   <- update(prop_0_l2_mod_dele,  . ~ . + l1:condition_sum:DELE)
  
  anova(prop_0_l2_mod_0, prop_0_l2_mod_dele, test = "Chisq") 
  # npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # prop_0_l2_mod_0       5 25570 25599 -12780    25560                     
  # prop_0_l2_mod_dele    6 25572 25606 -12780    25560 0.1428  1     0.7055
  
  anova(prop_0_l2_mod_0, prop_0_l2_mod_int_dele, test = "Chisq") # no interaction wm x condition
  # npar   AIC   BIC logLik deviance Chisq Df Pr(>Chisq)
  # prop_0_l2_mod_0         5 25570 25599 -12780    25560                    
  # prop_0_l2_mod_int_wm   8 25575 25620 -12780    25559 1.0387  3    0.7919
  
  l2_onset_v1_delewm_final <- prop_0_l2_mod_0
  
  summary(l2_onset_v1_delewm_final) 
  # Fixed effects:
  #                       Estimate Std. Error z value Pr(>|z|)    
  # (Intercept)           -0.30163    0.09241  -3.264   0.0011 **
 
 
  # ROUTE 3
  
  prop_0_l2_mod_use       <- update(prop_0_l2_mod_0,  . ~ . + percent_l2_week)
  prop_0_l2_mod_int_use   <- update(prop_0_l2_mod_use,  . ~ . + l1:condition_sum:percent_l2_week)
  
  anova(prop_0_l2_mod_0, prop_0_l2_mod_use, test = "Chisq") 
  # npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
  # prop_0_l2_mod_0      5 25570 25599 -12780    25560                       
  # prop_0_l2_mod_use    6 25567 25601 -12778    25555 5.0449  1     0.0247 *
  
  anova(prop_0_l2_mod_use, prop_0_l2_mod_int_use, test = "Chisq") # no interaction wm x condition
  # npar   AIC   BIC logLik deviance Chisq Df Pr(>Chisq)
  # prop_0_l2_mod_use        6 25567 25601 -12778    25555                     
  # prop_0_l2_mod_int_use    8 25569 25614 -12776    25553 2.5376  2     0.2812
  
  l2_onset_v1_use_final <- prop_0_l2_mod_use
  
  summary(l2_onset_v1_use_final)
  # Fixed effects:
  #                  Estimate Std. Error z value Pr(>|z|)  
  # (Intercept)     -0.073183   0.136114  -0.538   0.5908  
  # percent_l2_week -0.006121   0.002700  -2.267   0.0234 *
  
  
}
  # -----------------------------------------------------------------------------
  



  
  
  # Save models -----------------------------------------------------------------
  
  if(F) {
    
    saveRDS(mon_onset_v1_final, here("mods", "stress", 
                               "glmm", "onset_v1", "mon_onset_v1_final.rds"))
    saveRDS(prop_0_mon_mod_cond, here("mods", "stress", 
                                   "glmm", "onset_v1", "prop_0_mon_mod_cond.rds"))
    saveRDS(prop_0_mon_mod_wm, here("mods", "stress", 
                                      "glmm", "onset_v1", "prop_0_mon_mod_wm.rds"))
    saveRDS(prop_0_mon_mod_int, here("mods", "stress", 
                                      "glmm", "onset_v1", "prop_0_mon_mod_int.rds"))
    
    
    saveRDS(l2_onset_v1_delewm_final, here("mods", "stress", 
                                   "glmm", "onset_v1", "l2_onset_v1_delewm_final.rds"))
    saveRDS(l2_onset_v1_use_final, here("mods", "stress", 
                                  "glmm", "onset_v1", "l2_onset_v1_use_final.rds"))
    saveRDS(prop_0_l2_mod_l1, here("mods", "stress", 
                                     "glmm", "onset_v1", "prop_0_l2_mod_l1.rds"))
    saveRDS(prop_0_l2_mod_cond, here("mods", "stress", 
                                     "glmm", "onset_v1", "prop_0_l2_mod_cond.rds"))
    saveRDS(prop_0_l2_mod_wm, here("mods", "stress", 
                                     "glmm", "onset_v1", "prop_0_l2_mod_wm.rds"))
    saveRDS(prop_0_l2_mod_dele, here("mods", "stress", 
                                     "glmm", "onset_v1", "prop_0_l2_mod_dele.rds"))
    saveRDS(prop_0_l2_mod_int, here("mods", "stress", 
                                       "glmm", "onset_v1", "prop_0_l2_mod_int.rds"))
    saveRDS(prop_0_l2_mod_int_wm, here("mods", "stress", 
                                     "glmm", "onset_v1", "prop_0_l2_mod_int_wm.rds"))
    saveRDS(prop_0_l2_mod_int_dele, here("mods", "stress", 
                                       "glmm", "onset_v1", "prop_0_l2_mod_int_dele.rds"))
    saveRDS(prop_0_l2_mod_int_use, here("mods", "stress", 
                                       "glmm", "onset_v1", "prop_0_l2_mod_int_use.rds"))
  
    }




mon_onset_v1_final %>%
  ggplot(., aes(x = condition_sum, y = mean(target_prop))) + 
  geom_hline(yintercept = 0.5, color = "black", size = 0.75,
             lty = 3) + 
  stat_summary(fun.data = mean_cl_boot, geom = 'errorbar', 
               position = position_dodge(width = 0.5), 
               width = 0.35, color = 'grey40') + 
  stat_summary(fun.y = mean, geom = 'point', size = 4,
               position = position_dodge(width = 0.5)) + 
  coord_cartesian(ylim = c(0, 1)) + ylab('Target fixations') + xlab('Stress condition') + 
  ggtitle('Mean target fixations as a function of group\nand target type') +    
  theme_bw(base_size = 16, base_family = 'Times') -> stress_target_fix_mon_onset_v1

ggsave(paste0("./figs/stress/glmm/mon/stress_target_fix_mon_onset_v1.png"), stress_target_fix_mon_onset_v1, width = 150,
       height = 120, units = "mm", dpi = 600)


l2_onset_v1_delewm_final %>%
  ggplot(., aes(x = l1, y = mean(target_prop), 
                dodge = condition_sum, color = condition_sum, 
                l1 = interaction(l1, condition_sum))) +
  geom_hline(yintercept = 0.5, color = "black", size = 0.75, 
             lty = 3) + 
  stat_summary(fun.data = mean_cl_boot, geom = 'errorbar', 
               position = position_dodge(width = 0.5), 
               width = 0.35, color = 'grey40') + 
  stat_summary(fun.y = mean, geom = 'point', size = 4,
               position = position_dodge(width = 0.5)) + 
  coord_cartesian(ylim = c(0, 1)) + ylab('Target fixations') + xlab('L1') + 
  scale_x_discrete(labels = c("English", "Mandarin")) +
  ggtitle('Mean target fixations at V1 onset') + 
  # scale_color_gradient2(limits = c(-1, 1)) +
  theme_bw(base_size = 16, base_family = 'Times') -> stress_target_fix_l2_delewm_onset_v1

  ggsave(paste0("./figs/stress/glmm/l2/stress_target_fix_l2_delewm_onset_v1.png"), stress_target_fix_l2_delewm_onset_v1, width = 150,
         height = 120, units = "mm", dpi = 600)

  
  
  
  l2_onset_v1_use_final %>%
   
    ggplot(., aes(x = l1, y = mean(target_prop), 
                  color = percent_l2_week, 
                  l1 = interaction(l1, percent_l2_week))) +
    facet_grid(. ~ factor(condition_sum, levels=c(1, -1))) +
    geom_hline(yintercept = 0.5, color = "black", size = 0.75, 
               lty = 3) + 
    # stat_summary(fun.data = mean_cl_boot, geom = 'errorbar', 
    #              position = position_dodge(width = 0.5), 
    #              width = 0.35, color = 'grey40') + 
    stat_summary(fun.y = mean, geom = 'point', size = 4,
                 position = position_dodge(width = 0.5)) + 
    coord_cartesian(ylim = c(0, 1)) + ylab('Target fixations') + xlab('L1') + 
    scale_x_discrete(labels = c("English", "Mandarin")) +
    scale_color_continuous(name = 'L2 weekly use (%)') +
    ggtitle('Mean target fixations at V1 onset as a function of weekly L2 use (%)') + 
    
    theme_bw(base_size = 16, base_family = 'Times') -> stress_target_fix_l2_use_onset_v1
  
  ggsave(paste0("./figs/stress/glmm/l2/stress_target_fix_l2_use_onset_v1.png"), stress_target_fix_l2_use_onset_v1, width = 150,
         height = 120, units = "mm", dpi = 600)
# -----------------------------------------------------------------------------







# Onset c2

# Load data and models --------------------------------------------------------

source(here::here("scripts", "00_load_libs.R"))
# source(here::here("scripts", "02_load_data.R"))
stress10 <- read_csv(here::here("data", 'clean', 'stress_10ms_final_onset_c2.csv'))

# mon_onset_c2_final     <- readRDS(here("mods", 'stress', "glmm", "onset_c2", "mon_onset_c2_final.rds"))
# prop_0_mon_mod_cond  <- readRDS(here("mods", 'stress', "glmm", "onset_c2", "prop_0_mon_mod_cond.rds"))
# prop_0_mon_mod_wm    <- readRDS(here("mods", 'stress', "glmm", "onset_c2", "prop_0_mon_mod_wm.rds"))
# prop_0_mon_mod_int   <- readRDS(here("mods", 'stress', "glmm", "onset_c2", "prop_0_mon_mod_int.rds"))
# 
# l2_onset_c2_all_final <- readRDS(here("mods", 'stress', "glmm", "onset_c2", "l2_onset_c2_all_final.rds"))
# prop_0_l2_mod_cond <- readRDS(here("mods", 'stress', "glmm", "onset_c2", "prop_0_l2_mod_cond.rds"))
# prop_0_l2_mod_wm <- readRDS(here("mods", 'stress', "glmm", "onset_c2", "prop_0_l2_mod_wm.rds"))
# prop_0_l2_mod_dele <- readRDS(here("mods", 'stress', "glmm", "onset_c2", "prop_0_l2_mod_dele.rds"))
# prop_0_l2_mod_use <- readRDS(here("mods", 'stress', "glmm", "onset_c2", "prop_0_l2_mod_use.rds"))
# prop_0_l2_mod_int_wm  <- readRDS(here("mods", 'stress', "glmm", "onset_c2", "prop_0_l2_mod_int_wm.rds"))
# prop_0_l2_mod_int_dele  <- readRDS(here("mods", 'stress', "glmm", "onset_c2", "prop_0_l2_mod_int_dele.rds"))
# prop_0_l2_mod_use <- readRDS(, here("mods", "stress", "glmm", "onset_c2", "prop_0_l2_mod_use.rds"))
# prop_0_l2_mod_int_use  <- readRDS(here("mods", 'stress', "glmm", "onset_c2", "prop_0_l2_mod_int_use.rds"))

# -----------------------------------------------------------------------------



# Data prep -------------------------------------------------------------------

# Filter time course to offset of 1st syllable (time_zero == 20, to account for 200 ms to launch saccade)
# Create sum coded fixed factors (condition)


df_stress <- stress10 %>%
  filter(., time_zero == 20) %>%
  mutate(., condition_sum = if_else(cond == "1", 1, -1)) #   1 = present

stress_mon <- df_stress %>%
  filter(., l1 == 'es')

stress_l2 <- df_stress %>%
  filter(., l1 != 'es') %>%
  filter(., participant != 'ies04' & participant != 'ies17' & participant != 'ies28' & participant != 'aes32') %>%
  mutate(., l1 = fct_relevel(l1, "en", "ma"))



# -----------------------------------------------------------------------------




##########################     MONOLINGUALS    ##########################

# Random effects building -----------------------------------------------------

if(F) {
  
  prop_0_ranefA <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                           (1 | participant),
                         data = stress_mon, family = 'binomial',                 
                         control = glmerControl(optimizer = 'bobyqa'))
  
  prop_0_ranefB <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                           (1 | participant) +
                           (1 | target),
                         data = stress_mon, family = 'binomial', 
                         control = glmerControl(optimizer = 'bobyqa'))
  
  anova(prop_0_ranefA, prop_0_ranefB, refit = F) # keep intercept for target bc significant
  #                 Df    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
  # prop_0_ranefA    2 5378.7 5387.0 -2687.3   5374.7                         
  # prop_0_ranefB    3 5092.2 5104.7 -2543.1   5086.2 288.47  1  < 2.2e-16 ***
    
  prop_0_ranefC <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                           (1 + condition_sum | participant) +
                           (1 | target),
                         data = stress_mon, family = 'binomial', 
                         control = glmerControl(optimizer = 'bobyqa'))
  
  anova(prop_0_ranefB, prop_0_ranefC, refit = F) # keep slope for condition bc significant
  #                 Df    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
  # prop_0_ranefB    3 5092.2 5104.7 -2543.1   5086.2                     
  # prop_0_ranefC    5 4884.9 4905.7 -2437.4   4874.9 211.34  2  < 2.2e-16 ***
  
}

# -----------------------------------------------------------------------------





# Test fixed effects ----------------------------------------------------------

if(F) {
  prop_0_mon_mod_0 <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                              (1 + condition_sum | participant) +   
                              (1 | target),
                            data = stress_mon, family = 'binomial', 
                            control = glmerControl(optimizer = 'bobyqa'))
  
  prop_0_mon_mod_cond  <- update(prop_0_mon_mod_0, . ~ . + condition_sum)
  prop_0_mon_mod_wm    <- update(prop_0_mon_mod_cond,  . ~ . + ospan)
  prop_0_mon_mod_int  <- update(prop_0_mon_mod_wm,    . ~ . + condition_sum:ospan)
  
  anova(prop_0_mon_mod_0, prop_0_mon_mod_cond, test = "Chisq")  
  #                     npar    AIC    BIC  logLik deviance  Chisq  Df Pr(>Chisq)    
  # prop_0_mon_mod_0       5 4884.9 4905.7 -2437.4   4874.9          
  # prop_0_mon_mod_cond    6 4886.4 4911.4 -2437.2   4874.4 0.4957  1     0.4814
  
  anova(prop_0_mon_mod_0, prop_0_mon_mod_wm, test = "Chisq") 
  #                   npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
  # prop_0_mon_mod_0     5 4884.9 4905.7 -2437.4   4874.9                     
  # prop_0_mon_mod_wm    7 4887.8 4917.0 -2436.9   4873.8 1.0795  2     0.5829
  
  anova(prop_0_mon_mod_0, prop_0_mon_mod_int, test = "Chisq") 
  #                    npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
  # prop_0_mon_mod_0      5 4884.9 4905.7 -2437.4   4874.9                     
  # prop_0_mon_mod_int    8 4889.5 4922.9 -2436.7   4873.5 1.3905  3     0.7078
  
  
  mon_onset_c2_final <- prop_0_mon_mod_0
  
  summary(mon_onset_c2_final) 
  # Fixed effects:
  #                       Estimate Std. Error z value Pr(>|z|)    
  # (Intercept)             1.1825     0.2271   5.207 1.91e-07 ***
  
}

# -----------------------------------------------------------------------------






##########################     L2 SPEAKERS    ##########################

# Random effects building -----------------------------------------------------

if(F) {
  
  prop_0_ranefA <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                           (1 | participant),
                         data = stress_l2, family = 'binomial',                 
                         control = glmerControl(optimizer = 'bobyqa'))
  
  prop_0_ranefB <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                           (1 | participant) +
                           (1 | target),
                         data = stress_l2, family = 'binomial', 
                         control = glmerControl(optimizer = 'bobyqa'))
  
  anova(prop_0_ranefA, prop_0_ranefB, refit = F) # intercept
  #               npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
  # prop_0_ranefA    2 26124 26135 -13060    26120                     
  # prop_0_ranefB    3 25618 25635 -12806    25612 507.75  1  < 2.2e-16 ***
  
  prop_0_ranefC <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                           (1 + condition_sum | participant) +
                           (1 | target),
                         data = stress_l2, family = 'binomial', 
                         control = glmerControl(optimizer = 'bobyqa'))
  
  anova(prop_0_ranefB, prop_0_ranefC, refit = F) # slope
  #               npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
  # prop_0_ranefB    3 25618 25635 -12806    25612          
  # prop_0_ranefC    5 24740 24768 -12365    24730 882.41  2  < 2.2e-16 ***
  
}

# -----------------------------------------------------------------------------





# Test fixed effects ----------------------------------------------------------

if(F) {
  prop_0_l2_mod_0 <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                             (1 + condition_sum | participant) +   
                             (1 | target),
                           data = stress_l2, family = 'binomial', 
                           control = glmerControl(optimizer = 'bobyqa'))
  
  prop_0_l2_mod_l1    <- update(prop_0_l2_mod_0,     . ~ . + l1)
  prop_0_l2_mod_cond  <- update(prop_0_l2_mod_l1,    . ~ . + condition_sum)
  prop_0_l2_mod_int   <- update(prop_0_l2_mod_cond,  . ~ . + l1:condition_sum)
  
  anova(prop_0_l2_mod_0, prop_0_l2_mod_l1, test = "Chisq")
  #                  npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # prop_0_l2_mod_0     5 24740 24768 -12365    24730                    
  # prop_0_l2_mod_l1    6 24742 24775 -12365    24730 0.083  1     0.7733
  
  anova(prop_0_l2_mod_0, prop_0_l2_mod_cond, test = "Chisq")   
  #                    npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # prop_0_l2_mod_0       5 24740 24768 -12365    24730                    
  # prop_0_l2_mod_cond    7 24743 24782 -12364    24729 0.968  2     0.6163
  
  anova(prop_0_l2_mod_0, prop_0_l2_mod_int, test = "Chisq") 
  #                   npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # prop_0_l2_mod_0      5 24740 24768 -12365    24730                     
  # prop_0_l2_mod_int    8 24744 24789 -12364    24728 1.0161  3     0.7973
  
  # ROUTE 1
  
  prop_0_l2_mod_wm       <- update(prop_0_l2_mod_0,  . ~ . + ospan)
  prop_0_l2_mod_int_wm   <- update(prop_0_l2_mod_wm,    . ~ . + l1:condition_sum:ospan)
  
  anova(prop_0_l2_mod_0, prop_0_l2_mod_wm, test = "Chisq") 
  #                  npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # prop_0_l2_mod_0     5 24740 24768 -12365    24730                     
  # prop_0_l2_mod_wm    6 24742 24775 -12365    24730 0.0164  1     0.8982
  
  anova(prop_0_l2_mod_0, prop_0_l2_mod_int_wm, test = "Chisq") # no interaction wm x condition
  #                      npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # prop_0_l2_mod_0         5 24740 24768 -12365    24730                     
  # prop_0_l2_mod_int_wm    8 24744 24788 -12364    24728 2.0041  3     0.5716
  
  l2_onset_c2_all_final <- prop_0_l2_mod_0 # for DELE and WM
  
  summary(l2_onset_c2_all_final) 
  # Fixed effects:
  #                       Estimate Std. Error z value Pr(>|z|)    
  # (Intercept)           -0.29616    0.09794  -3.024  0.00249 **
  
  
  # ROUTE 2
  
  prop_0_l2_mod_dele       <- update(prop_0_l2_mod_0,  . ~ . + DELE_z)
  prop_0_l2_mod_int_dele   <- update(prop_0_l2_mod_dele,  . ~ . + l1:condition_sum:DELE_z)
  
  anova(prop_0_l2_mod_0, prop_0_l2_mod_dele, test = "Chisq") 
  # npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # prop_0_l2_mod_0       5 24740 24768 -12365    24730                    
  # prop_0_l2_mod_dele    6 24742 24775 -12365    24730 0.034  1     0.8538
  
  anova(prop_0_l2_mod_0, prop_0_l2_mod_int_dele, test = "Chisq") # no interaction wm x condition
  #                      npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # prop_0_l2_mod_0         5 24740 24768 -12365    24730                     
  # prop_0_l2_mod_int_dele    8 24745 24790 -12365    24729 0.2945  3     0.9611
  
  
  # ROUTE 3
  
  prop_0_l2_mod_use  <- update(prop_0_l2_mod_0,  . ~ . + percent_l2_week)
  
  anova(prop_0_l2_mod_0, prop_0_l2_mod_use, test = "Chisq") 
  #                   npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
  # prop_0_l2_mod_0      5 24740 24768 -12365    24730                       
  # prop_0_l2_mod_use    6 24736 24769 -12362    24724 5.9657  1    0.01459 *
  
  prop_0_l2_mod_int_use   <- update(prop_0_l2_mod_use,  . ~ . + l1:condition_sum:percent_l2_week)
  
  anova(prop_0_l2_mod_use, prop_0_l2_mod_int_use, test = "Chisq")
  #                       npar   AIC   BIC logLik deviance Chisq Df Pr(>Chisq)
  # prop_0_l2_mod_use        6 24736 24769 -12362    24724                     
  # prop_0_l2_mod_int_use    8 24737 24782 -12360    24721 2.7048  2     0.2586
  
  l2_onset_c2_use_final <- prop_0_l2_mod_use
  
  summary(l2_onset_c2_use_final)
  # Estimate Std. Error z value Pr(>|z|)  
  # (Intercept)     -0.030531   0.144633  -0.211   0.8328  
  # percent_l2_week -0.006954   0.002813  -2.472   0.0134 *
}
# -----------------------------------------------------------------------------






# Save models -----------------------------------------------------------------

if(F) {
  
  saveRDS(mon_onset_c2_final, here("mods", "stress", 
                                 "glmm", "onset_c2", "mon_onset_c2_final.rds"))
  saveRDS(prop_0_mon_mod_cond, here("mods", "stress", 
                                    "glmm", "onset_c2", "prop_0_mon_mod_cond.rds"))
  saveRDS(prop_0_mon_mod_wm, here("mods", "stress", 
                                  "glmm", "onset_c2", "prop_0_mon_mod_wm.rds"))
  saveRDS(prop_0_mon_mod_int, here("mods", "stress", 
                                   "glmm", "onset_c2", "prop_0_mon_mod_int.rds"))
  
  
  saveRDS(l2_onset_c2_all_final, here("mods", "stress", 
                                "glmm", "onset_c2", "l2_onset_c2_final.rds"))
  saveRDS(prop_0_l2_mod_l1, here("mods", "stress", 
                                 "glmm", "onset_c2", "prop_0_l2_mod_l1.rds"))
  saveRDS(prop_0_l2_mod_cond, here("mods", "stress", 
                                   "glmm", "onset_c2", "prop_0_l2_mod_cond.rds"))
  saveRDS(prop_0_l2_mod_wm, here("mods", "stress", 
                                 "glmm", "onset_c2", "prop_0_l2_mod_wm.rds"))
  saveRDS(prop_0_l2_mod_dele, here("mods", "stress", 
                                   "glmm", "onset_c2", "prop_0_l2_mod_dele.rds"))
  saveRDS(prop_0_l2_mod_int, here("mods", "stress", 
                                  "glmm", "onset_c2", "prop_0_l2_mod_int.rds"))
  saveRDS(prop_0_l2_mod_int_wm, here("mods", "stress", 
                                     "glmm", "prop_0_l2_mod_int_wm.rds"))
  saveRDS(prop_0_l2_mod_int_dele, here("mods", "stress", 
                                       "glmm", "onset_c2", "prop_0_l2_mod_int_dele.rds"))
  saveRDS(prop_0_l2_mod_use, here("mods", "stress", 
                                  "glmm", "onset_c2", "prop_0_l2_mod_use.rds"))
  saveRDS(prop_0_l2_mod_int_use, here("mods", "stress", 
                                      "glmm", "onset_c2", "prop_0_l2_mod_int_use.rds"))
  saveRDS(l2_onset_c2_use_final, here("mods", "stress", 
                                       "glmm", "onset_c2", "l2_onset_c2_use_final.rds"))
  
}


mon_onset_c2_final %>%
  ggplot(., aes(x = condition_sum, y = mean(target_prop))) + 
  geom_hline(yintercept = 0.5, color = "black", size = 0.75,
             lty = 3) + 
  stat_summary(fun.data = mean_cl_boot, geom = 'errorbar', 
               position = position_dodge(width = 0.5), 
               width = 0.35, color = 'grey40') + 
  stat_summary(fun.y = mean, geom = 'point', size = 4,
               position = position_dodge(width = 0.5)) + 
  coord_cartesian(ylim = c(0, 1)) + ylab('Target fixations') + xlab('Stress condition') + 
  ggtitle('Mean target fixations as a function of group\nand target type') +    
  theme_bw(base_size = 16, base_family = 'Times') -> stress_target_fix_mon_onset_c2

ggsave(paste0("./figs/stress/glmm/mon/stress_target_fix_mon_onset_c2.png"), stress_target_fix_mon_onset_c2, width = 150,
       height = 120, units = "mm", dpi = 600)


l2_onset_c2_all_final %>%
  ggplot(., aes(x = l1, y = mean(target_prop), 
                dodge = condition_sum, color = condition_sum, 
                l1 = interaction(l1, condition_sum))) +
  geom_hline(yintercept = 0.5, color = "black", size = 0.75, 
             lty = 3) + 
  stat_summary(fun.data = mean_cl_boot, geom = 'errorbar', 
               position = position_dodge(width = 0.5), 
               width = 0.35, color = 'grey40') + 
  stat_summary(fun.y = mean, geom = 'point', size = 4,
               position = position_dodge(width = 0.5)) + 
  coord_cartesian(ylim = c(0, 1)) + ylab('Target fixations') + xlab('L1') + 
  scale_x_discrete(labels = c("English", "Mandarin")) +
  ggtitle('Mean target fixations at V1 onset') + 
  # scale_color_gradient2(limits = c(-1, 1)) +
  theme_bw(base_size = 16, base_family = 'Times') -> stress_target_fix_l2_all_onset_c2

ggsave(paste0("./figs/stress/glmm/l2/stress_target_fix_l2_all_onset_c2.png"), stress_target_fix_l2_all_onset_c2, width = 150,
       height = 120, units = "mm", dpi = 600)
# -----------------------------------------------------------------------------








# Onset c3

# Load data and models --------------------------------------------------------

source(here::here("scripts", "00_load_libs.R"))
# source(here::here("scripts", "02_load_data.R"))
stress10 <- read_csv(here::here("data", 'clean', 'stress_10ms_final_onset_c3.csv'))

# prop_0_mon_mod_0     <- readRDS(here("mods", 'stress', "glmm", "onset_c3", "mon_onset_v1_final.rds"))
# prop_0_mon_mod_cond  <- readRDS(here("mods", 'stress', "glmm", "onset_c3", "prop_0_mon_mod_cond.rds"))
# prop_0_mon_mod_wm    <- readRDS(here("mods", 'stress', "glmm", "onset_c3", "prop_0_mon_mod_wm.rds"))
# prop_0_mon_mod_int   <- readRDS(here("mods", 'stress', "glmm", "onset_c3", "prop_0_mon_mod_int.rds"))
# 
# prop_0_l2_mod_0 <- readRDS(here("mods", 'stress', "glmm", "onset_c3", "l2_onset_v1_delewm_final.rds"))
# prop_0_l2_mod_cond <- readRDS(here("mods", 'stress', "glmm", "onset_c3", "prop_0_l2_mod_cond.rds"))
# prop_0_l2_mod_wm <- readRDS(here("mods", 'stress', "glmm", "onset_c3", "prop_0_l2_mod_wm.rds"))
# prop_0_l2_mod_dele <- readRDS(here("mods", 'stress', "glmm", "onset_c3", "prop_0_l2_mod_dele.rds"))
# prop_0_l2_mod_use <- readRDS(here("mods", 'stress', "glmm", "onset_c3", "12_onset_v1_use_final.rds"))
# prop_0_l2_mod_int_wm  <- readRDS(here("mods", 'stress', "glmm", "onset_c3", "prop_0_l2_mod_int_wm.rds"))
# prop_0_l2_mod_int_dele  <- readRDS(here("mods", 'stress', "glmm", "onset_c3", "prop_0_l2_mod_int_dele.rds"))
# prop_0_l2_mod_use <- readRDS(here("mods", "stress", "glmm", "onset_c3", "prop_0_l2_mod_use.rds"))
# prop_0_l2_mod_int_use  <- readRDS(here("mods", 'stress', "glmm", "onset_c3", "prop_0_l2_mod_int_use.rds"))

# -----------------------------------------------------------------------------



# Data prep -------------------------------------------------------------------

# Filter time course to offset of 1st syllable (time_zero == 20, to account for 200 ms to launch saccade)
# Create sum coded fixed factors (condition)


df_stress <- stress10 %>%
  filter(., time_zero == 20) %>%
  mutate(., condition_sum = if_else(cond == "1", 1, -1)) #   1 = present

stress_mon <- df_stress %>%
  filter(., l1 == 'es')

stress_l2 <- df_stress %>%
  filter(., l1 != 'es') %>%
  filter(., participant != 'ies04' & participant != 'ies17' & participant != 'ies28' & participant != 'aes32') %>%
  mutate(., l1 = fct_relevel(l1, "en", "ma"))



# -----------------------------------------------------------------------------




##########################     MONOLINGUALS    ##########################

# Random effects building -----------------------------------------------------

if(F) {
  
  prop_0_ranefA <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                           (1 | participant),
                         data = stress_mon, family = 'binomial',                 
                         control = glmerControl(optimizer = 'bobyqa'))
  
  prop_0_ranefB <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                           (1 | participant) +
                           (1 | target),
                         data = stress_mon, family = 'binomial', 
                         control = glmerControl(optimizer = 'bobyqa'))
  
  anova(prop_0_ranefA, prop_0_ranefB, refit = F) 
  #                 Df    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
  # prop_0_ranefA    2 3700.7 3709.0 -1848.3   3696.7                         
  # prop_0_ranefB    3 3409.1 3421.7 -1701.6   3403.1 293.53  1  < 2.2e-16 ***
  
  prop_0_ranefC <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                           (1 + condition_sum | participant) +
                           (1 | target),
                         data = stress_mon, family = 'binomial', 
                         control = glmerControl(optimizer = 'bobyqa'))
  
  anova(prop_0_ranefB, prop_0_ranefC, refit = F) 
  #                 Df    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
  # prop_0_ranefB    3 3409.1 3421.7 -1701.6   3403.1                         
  # prop_0_ranefC    5 3280.8 3301.6 -1635.4   3270.8 132.38  2  < 2.2e-16 ***
  
}

# -----------------------------------------------------------------------------





# Test fixed effects ----------------------------------------------------------

if(F) {
  prop_0_mon_mod_0 <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                              (1 + condition_sum | participant) +   
                              (1 | target),
                            data = stress_mon, family = 'binomial', 
                            control = glmerControl(optimizer = 'bobyqa'))
  
  prop_0_mon_mod_cond  <- update(prop_0_mon_mod_0, . ~ . + condition_sum)
  prop_0_mon_mod_wm    <- update(prop_0_mon_mod_cond,  . ~ . + ospan)
  prop_0_mon_mod_int  <- update(prop_0_mon_mod_wm,    . ~ . + condition_sum:ospan)
  
  anova(prop_0_mon_mod_0, prop_0_mon_mod_cond, test = "Chisq")  
  #                     npar    AIC    BIC  logLik deviance  Chisq  Df Pr(>Chisq)    
  # prop_0_mon_mod_0       5 3280.8 3301.6 -1635.4   3270.8                     
  # prop_0_mon_mod_cond    6 3282.5 3307.5 -1635.2   3270.5 0.2897  1     0.5904
  
  anova(prop_0_mon_mod_0, prop_0_mon_mod_wm, test = "Chisq") 
  #                   npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
  # prop_0_mon_mod_0     5 3280.8 3301.6 -1635.4   3270.8                    
  # prop_0_mon_mod_wm    7 3284.4 3313.6 -1635.2   3270.4 0.356  2      0.837
  
  anova(prop_0_mon_mod_0, prop_0_mon_mod_int, test = "Chisq") 
  #                    npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)
  # prop_0_mon_mod_0      5 3280.8 3301.6 -1635.4   3270.8                    
  # prop_0_mon_mod_int    8 3285.1 3318.5 -1634.6   3269.1 1.643  3     0.6497
  
  mon_onset_c3_final <- prop_0_mon_mod_0
  
  summary(mon_onset_c3_final) 
  # Fixed effects:
  #                       Estimate Std. Error z value Pr(>|z|)    
  # (Intercept)             2.7658     0.3631   7.617  2.6e-14 ***
  
}

# -----------------------------------------------------------------------------






##########################     L2 SPEAKERS    ##########################

# Random effects building -----------------------------------------------------

if(F) {
  
  prop_0_ranefA <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                           (1 | participant),
                         data = stress_l2, family = 'binomial',                 
                         control = glmerControl(optimizer = 'bobyqa'))
  
  prop_0_ranefB <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                           (1 | participant) +
                           (1 | target),
                         data = stress_l2, family = 'binomial', 
                         control = glmerControl(optimizer = 'bobyqa'))
  
  anova(prop_0_ranefA, prop_0_ranefB, refit = F) # intercept
  #               npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
  # prop_0_ranefA    2 23424 23435 -11710    23420                     
  # prop_0_ranefB    3 23002 23018 -11498    22996 424.46  1  < 2.2e-16 ***
  
  prop_0_ranefC <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                           (1 + condition_sum | participant) +
                           (1 | target),
                         data = stress_l2, family = 'binomial', 
                         control = glmerControl(optimizer = 'bobyqa'))
  
  anova(prop_0_ranefB, prop_0_ranefC, refit = F) # slope
  #               npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)    
  # prop_0_ranefB    3 23002 23018 -11498    22996                     
  # prop_0_ranefC    5 22126 22154 -11058    22116 879.51  2  < 2.2e-16 ***
  
}

# -----------------------------------------------------------------------------





# Test fixed effects ----------------------------------------------------------

if(F) {
  prop_0_l2_mod_0 <- glmer(cbind(target_count, 10 - target_count) ~ 1 +
                             (1 + condition_sum | participant) +   
                             (1 | target),
                           data = stress_l2, family = 'binomial', 
                           control = glmerControl(optimizer = 'bobyqa'))
  
  prop_0_l2_mod_l1    <- update(prop_0_l2_mod_0,     . ~ . + l1)
  prop_0_l2_mod_cond  <- update(prop_0_l2_mod_l1,    . ~ . + condition_sum)
  prop_0_l2_mod_int   <- update(prop_0_l2_mod_cond,  . ~ . + l1:condition_sum)
  
  anova(prop_0_l2_mod_0, prop_0_l2_mod_l1, test = "Chisq")
  #                  npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # prop_0_l2_mod_0     5 22126 22154 -11058    22116                       
  # prop_0_l2_mod_l1    6 22125 22158 -11056    22113 3.3256  1    0.06821 .
  
  anova(prop_0_l2_mod_0, prop_0_l2_mod_cond, test = "Chisq")   
  #                    npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # prop_0_l2_mod_0       5 22126 22154 -11058    22116                       
  # prop_0_l2_mod_cond    7 22123 22162 -11055    22109 6.9681  2    0.03068 *
  
  anova(prop_0_l2_mod_cond, prop_0_l2_mod_int, test = "Chisq") 
  #                    npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # prop_0_l2_mod_cond    7 22123 22162 -11055    22109                    
  # prop_0_l2_mod_int     8 22125 22170 -11054    22109  0.38  1     0.5376
  
  # ROUTE 1
  
  prop_0_l2_mod_wm       <- update(prop_0_l2_mod_cond,  . ~ . + ospan)
  prop_0_l2_mod_int_wm   <- update(prop_0_l2_mod_wm,    . ~ . + l1:condition_sum:ospan)
  
  anova(prop_0_l2_mod_cond, prop_0_l2_mod_wm, test = "Chisq") 
  #                  npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # prop_0_l2_mod_cond    7 22123 22162 -11055    22109                     
  # prop_0_l2_mod_wm      8 22124 22169 -11054    22108 0.9879  1     0.3203
  
  anova(prop_0_l2_mod_cond, prop_0_l2_mod_int_wm, test = "Chisq") # no interaction wm x condition
  #                      npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # prop_0_l2_mod_cond      7 22123 22162 -11055    22109                     
  # prop_0_l2_mod_int_wm   10 22125 22181 -11053    22105 3.8324  3     0.2801
  
  l2_onset_c3_wm_final <- prop_0_l2_mod_cond
  
  summary(l2_onset_c3_wm_final) 
  # Fixed effects:
  #               Estimate Std. Error z value Pr(>|z|)    
  # (Intercept)    0.99976    0.12059   8.290   <2e-16 ***
  # l1ma          -0.25833    0.13990  -1.846   0.0648 .  
  # condition_sum -0.15232    0.07945  -1.917   0.0552 .
  
  
  # ROUTE 2
  
  prop_0_l2_mod_dele       <- update(prop_0_l2_mod_cond,  . ~ . + DELE_z)
  prop_0_l2_mod_int_dele   <- update(prop_0_l2_mod_dele,  . ~ . + l1:condition_sum:DELE_z)
  
  anova(prop_0_l2_mod_cond, prop_0_l2_mod_dele, test = "Chisq") 
  #                    npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # prop_0_l2_mod_cond    7 22123 22162 -11055    22109                        
  # prop_0_l2_mod_dele    8 22112 22157 -11048    22096 13.44  1  0.0002463 ***
  
  anova(prop_0_l2_mod_dele, prop_0_l2_mod_int_dele, test = "Chisq") # no interaction wm x condition
  #                        npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # prop_0_l2_mod_dele        8 22112 22157 -11048    22096                      
  # prop_0_l2_mod_int_dele   10 22109 22165 -11044    22089 7.057  2    0.02935 *
  
  l2_onset_c3_dele_final <- prop_0_l2_mod_int_dele
  
  summary(l2_onset_c3_dele_final) 
  #                           Estimate Std. Error z value Pr(>|z|)    
  # (Intercept)                1.10415    0.11957   9.235  < 2e-16 ***
  # l1ma                      -0.26933    0.13337  -2.019 0.043441 *  
  # condition_sum             -0.18779    0.08133  -2.309 0.020934 *  
  # DELE_z                     0.31939    0.08338   3.830 0.000128 ***
  # l1en:condition_sum:DELE_z -0.25396    0.09595  -2.647 0.008129 ** 
  # l1ma:condition_sum:DELE_z  0.03215    0.09962   0.323 0.746920 
  
  
  stress_l2 %>%
    group_by(condition_sum) %>%
    summarize(gaze_fix = mean(target_count)) %>% print()
  #   condition_sum gaze_fix
  # 1            -1     6.94
  # 2             1     6.47
  
  # ROUTE 3
  
  prop_0_l2_mod_use  <- update(prop_0_l2_mod_cond,  . ~ . + use_z)
  prop_0_l2_mod_int_use   <- update(prop_0_l2_mod_use,  . ~ . + l1:condition_sum:use_z)
  
  anova(prop_0_l2_mod_cond, prop_0_l2_mod_use, test = "Chisq") 
  #                    npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)  
  # prop_0_l2_mod_cond    7 22123 22162 -11055    22109                       
  # prop_0_l2_mod_use     8 22118 22163 -11051    22102 7.449  1   0.006347 **
  
  anova(prop_0_l2_mod_use, prop_0_l2_mod_int_use, test = "Chisq")
  #                       npar   AIC   BIC logLik deviance  Chisq Df Pr(>Chisq)
  # prop_0_l2_mod_use        8 22118 22163 -11051    22102                     
  # prop_0_l2_mod_int_use   10 22121 22177 -11051    22101 0.5688  2     0.7525
  
  l2_onset_c3_use_final <- prop_0_l2_mod_use
  
  summary(l2_onset_c3_use_final)
  #               Estimate Std. Error z value Pr(>|z|)    
  # (Intercept)    1.11964    0.12606   8.882  < 2e-16 ***
  # l1ma          -0.32312    0.13791  -2.343  0.01913 *  
  # condition_sum -0.15396    0.07939  -1.939  0.05248 .  
  # use_z          0.26537    0.09587   2.768  0.00564 **  
    
}

saveRDS(l2_onset_c3_dele_final, here("mods", "stress", 
                                   "glmm", "onset_c3", "l2_onset_c3_dele_final_z.rds"))
saveRDS(l2_onset_c3_use_final, here("mods", "stress", 
                                     "glmm", "onset_c3", "l2_onset_c3_use_final_z.rds"))
                                     
# -----------------------------------------------------------------------------






# Save models -----------------------------------------------------------------

if(F) {
  
  saveRDS(mon_onset_c3_final, here("mods", "stress", 
                                 "glmm", "onset_c3", "mon_onset_c3_final.rds"))
  saveRDS(prop_0_mon_mod_cond, here("mods", "stress", 
                                    "glmm", "onset_c3", "prop_0_mon_mod_cond.rds"))
  saveRDS(prop_0_mon_mod_wm, here("mods", "stress", 
                                  "glmm", "onset_c3", "prop_0_mon_mod_wm.rds"))
  saveRDS(prop_0_mon_mod_int, here("mods", "stress", 
                                   "glmm", "onset_c3", "prop_0_mon_mod_int.rds"))
  
  
  saveRDS(l2_onset_c3_wm_final, here("mods", "stress", 
                                "glmm", "onset_c3", "l2_onset_c3_wm_final.rds"))
  saveRDS(l2_onset_c3_dele_final, here("mods", "stress", 
                                     "glmm", "onset_c3", "l2_onset_c3_dele_final.rds"))
  saveRDS(l2_onset_c3_use_final, here("mods", "stress", 
                                     "glmm", "onset_c3", "l2_onset_c3_use_final.rds"))
  saveRDS(prop_0_l2_mod_l1, here("mods", "stress", 
                                 "glmm", "onset_c3", "prop_0_l2_mod_l1.rds"))
  saveRDS(prop_0_l2_mod_cond, here("mods", "stress", 
                                   "glmm", "onset_c3", "prop_0_l2_mod_cond.rds"))
  saveRDS(prop_0_l2_mod_int, here("mods", "stress", 
                                  "glmm", "onset_c3", "prop_0_l2_mod_int.rds"))
  
  saveRDS(prop_0_l2_mod_wm, here("mods", "stress", 
                                 "glmm", "onset_c3", "prop_0_l2_mod_wm.rds"))
  saveRDS(prop_0_l2_mod_int_wm, here("mods", "stress", 
                                     "glmm", "prop_0_l2_mod_int_wm.rds"))
  saveRDS(prop_0_l2_mod_dele, here("mods", "stress", 
                                   "glmm", "onset_c3", "prop_0_l2_mod_dele.rds"))
  saveRDS(prop_0_l2_mod_int_dele, here("mods", "stress", 
                                       "glmm", "onset_c3", "prop_0_l2_mod_int_dele.rds"))
  saveRDS(prop_0_l2_mod_use, here("mods", "stress", 
                                  "glmm", "onset_c3", "prop_0_l2_mod_use.rds"))
  saveRDS(prop_0_l2_mod_int_use, here("mods", "stress", 
                                      "glmm", "onset_c3", "prop_0_l2_mod_int_use.rds"))
  
  
}


mon_onset_c3_final %>%
  ggplot(., aes(x = condition_sum, y = mean(target_prop))) + 
  geom_hline(yintercept = 0.5, color = "black", size = 0.75,
             lty = 3) + 
  stat_summary(fun.data = mean_cl_boot, geom = 'errorbar', 
               position = position_dodge(width = 0.5), 
               width = 0.35, color = 'grey40') + 
  stat_summary(fun.y = mean, geom = 'point', size = 4,
               position = position_dodge(width = 0.5)) + 
  coord_cartesian(ylim = c(0, 1)) + ylab('Target fixations') + xlab('Stress condition') + 
  ggtitle('Mean target fixations as a function of group\nand target type') +    
  theme_bw(base_size = 16, base_family = 'Times') -> stress_target_fix_mon_onset_c3

ggsave(paste0("./figs/stress/glmm/mon/stress_target_fix_mon_onset_c3.png"), stress_target_fix_mon_onset_c3, width = 150,
       height = 120, units = "mm", dpi = 600)



l2_onset_c3_use_final %>%
  ggplot(., aes(x = l1, y = mean(target_prop), 
                color = percent_l2_week, 
                l1 = interaction(l1, percent_l2_week))) +
  facet_grid(. ~ factor(condition_sum, levels=c(1, -1))) +
  geom_hline(yintercept = 0.5, color = "black", size = 0.75, 
             lty = 3) + 
  # stat_summary(fun.data = mean_cl_boot, geom = 'errorbar', 
  #              position = position_dodge(width = 0.5), 
  #              width = 0.35, color = 'grey40') + 
  stat_summary(fun.y = mean, geom = 'point', size = 4,
               position = position_dodge(width = 0.5)) + 
  coord_cartesian(ylim = c(0, 1)) + ylab('Target fixations') + xlab('L1') + 
  scale_x_discrete(labels = c("English", "Mandarin")) +
  scale_color_continuous(name = 'L2 weekly use (%)') +
  ggtitle('Mean target fixations at V1 onset as a function of weekly L2 use (%)') + 
  theme_bw(base_size = 16, base_family = 'Times') -> stress_target_fix_l2_use_onset_c3

ggsave(paste0("./figs/stress/glmm/l2/stress_target_fix_l2_use_onset_c3.png"), stress_target_fix_l2_use_onset_c3, width = 150,
       height = 120, units = "mm", dpi = 600)



l2_onset_c3_wm_final %>%
  ggplot(., aes(x = l1, y = mean(target_prop), 
                color = ospan, 
                l1 = interaction(l1, ospan))) +
  facet_grid(. ~ factor(condition_sum, levels=c(1, -1))) +
  geom_hline(yintercept = 0.5, color = "black", size = 0.75, 
             lty = 3) + 
  # stat_summary(fun.data = mean_cl_boot, geom = 'errorbar', 
  #              position = position_dodge(width = 0.5), 
  #              width = 0.35, color = 'grey40') + 
  stat_summary(fun.y = mean, geom = 'point', size = 4,
               position = position_dodge(width = 0.5)) + 
  coord_cartesian(ylim = c(0, 1)) + ylab('Target fixations') + xlab('L1') + 
  scale_x_discrete(labels = c("English", "Mandarin")) +
  scale_color_continuous(name = 'Verbal WM') +
  ggtitle('Mean target fixations at V1 onset as a function of WM') + 
  theme_bw(base_size = 16, base_family = 'Times') -> stress_target_fix_l2_wm_onset_c3

ggsave(paste0("./figs/stress/glmm/l2/stress_target_fix_l2_wm_onset_c3.png"), stress_target_fix_l2_wm_onset_c3, width = 150,
       height = 120, units = "mm", dpi = 600)



l2_onset_c3_dele_final %>%
  ggplot(., aes(x = l1, y = mean(target_prop), 
                color = DELE, 
                l1 = interaction(l1, DELE))) +
  facet_grid(. ~ factor(condition_sum, levels=c(1, -1))) +
  geom_hline(yintercept = 0.5, color = "black", size = 0.75, 
             lty = 3) + 
  # stat_summary(fun.data = mean_cl_boot, geom = 'errorbar', 
  #              position = position_dodge(width = 0.5), 
  #              width = 0.35, color = 'grey40') + 
  stat_summary(fun.y = mean, geom = 'point', size = 4,
               position = position_dodge(width = 0.5)) + 
  coord_cartesian(ylim = c(0, 1)) + ylab('Target fixations') + xlab('L1') + 
  scale_x_discrete(labels = c("English", "Mandarin")) +
  scale_color_continuous(name = 'Proficiency') +
  ggtitle('Mean target fixations at V1 onset as a function of proficiency') + 
  theme_bw(base_size = 16, base_family = 'Times') -> stress_target_fix_l2_dele_onset_c3

ggsave(paste0("./figs/stress/glmm/l2/stress_target_fix_l2_dele_onset_c3.png"), stress_target_fix_l2_dele_onset_c3, width = 150,
       height = 120, units = "mm", dpi = 600)

# -----------------------------------------------------------------------------




# Model descriptives ----------------------------------------------------------
# 
MuMIn::r.squaredGLMM(prop_0_mod_final) 
#                   R2m       R2c
# theoretical 0.1644337 0.8341553
# delta       0.1528810 0.7755495

summary(prop_0_mod_final)
#                         Estimate Std. Error z value Pr(>|z|)    
#   (Intercept)           2.2390     0.1803  12.420  < 2e-16 ***
#   groupaes             -1.0022     0.2271  -4.413 1.02e-05 ***
#   groupies             -1.4069     0.2242  -6.275 3.49e-10 ***
#   groupams             -1.3068     0.2259  -5.786 7.21e-09 ***
#   groupims             -1.6783     0.2256  -7.438 1.02e-13 ***
  
confint(prop_0_mod_final, method = "Wald")
#                            2.5 %       97.5 %
# (Intercept)           1.885679  2.5923575
# groupaes             -1.447381 -0.5571072
# groupies             -1.846331 -0.9674771
# groupams             -1.749521 -0.8641394
# groupims             -2.120526 -1.2360267



# Calculate mean target fixation as a function of group, condition, 
# for each participant. We will plot the mean and calculate the 
# bootstrapped 95% confidence interval and plot it all. 
prop_0_mod_final %>%
  group_by(., group, condition_sum, participant) %>%
  #  summarise(., meanFix = mean(target_prop)) %>%
  ggplot(., aes(x = group, y = mean(target_prop), 
                dodge = condition_sum, color = condition_sum, 
                group = interaction(group, condition_sum))) +
  geom_hline(yintercept = 0.5, color = "black", size = 0.75, 
             lty = 3) + 
  stat_summary(fun.data = mean_cl_boot, geom = 'errorbar', 
               position = position_dodge(width = 0.5), 
               width = 0.35, color = 'grey40') + 
  stat_summary(fun.y = mean, geom = 'point', size = 4,
               position = position_dodge(width = 0.5)) + 
  coord_cartesian(ylim = c(0, 1)) + ylab('Target fixations') + xlab('Group') + 
  #  scale_x_discrete(labels = c("mon", "aes", "ies", "ams", "ims")) +         
  ggtitle('Mean target fixations as a function of group\nand target type') +    
  #  scale_color_brewer(palette = "Set1", name = '', labels = c('Present', 'Preterit')) +       
  theme_bw(base_size = 16, base_family = 'Times') -> stress_rel_target_fix

# Graph to check the effects of WM

prop_0_mod_final %>%
#  group_by(., group, condition_sum, participant) %>%
  ggplot(., aes(x = group, y = mean(target_prop), 
                dodge = WM_set, color = WM_set,
                group = interaction(group, WM_set))) +
  geom_hline(yintercept = 0.5, color = "black", size = 0.75, 
             lty = 3) + 
  stat_summary(fun.data = mean_cl_boot, geom = 'errorbar', 
               position = position_dodge(width = 0.5), 
               width = 0.35, color = 'grey40') + 
  stat_summary(fun.y = mean, geom = 'point', size = 4,
               position = position_dodge(width = 0.5)) + 
  coord_cartesian(ylim = c(0, 1)) + ylab('Target fixations') + xlab('Group') + 
  #  scale_x_discrete(labels = c("mon", "aes", "ies", "ams", "ims")) +       
  ggtitle('Mean target fixations as a function of group and verbal WM') + 
  # scale_color_brewer(palette = "Set1", name = '', labels = c('Present', 'Preterit')) + 
  theme_bw(base_size = 16, base_family = 'Times') -> stress_rel_coda_fix


