# Morphosyntactic predictability: tidy stress data ----------------------------
#
# This script will load and tidy the raw eye tracking data
# with 10 ms bins and save the output to data/clean
#
# Last update: 09/03/2020
#
# -----------------------------------------------------------------------------




# Source libs -----------------------------------------------------------------

source(here::here("scripts", "00_load_libs.R"))



# -----------------------------------------------------------------------------

# load data
stress10 <- read.delim(here("data", "stress_10bin.txt"))

# Check gaze fixation columns have different values
unique(stress10$AVERAGE_IA_1_SAMPLE_COUNT)  # looking at target according to IA_#_ID
unique(stress10$AVERAGE_IA_2_SAMPLE_COUNT)  # looking at distractor
unique(stress10$AVERAGE_IA_0_SAMPLE_COUNT)  # elsewhere


# Tidy data -------------------------------------------------------------------

# Read data
stress10 <- stress10 %>%
  
  # create variable group
  separate(., col = RECORDING_SESSION_LABEL,
           into = c("group", "group_member"),
           sep = 3,
           remove = FALSE) %>%
  
  #select and rename variables of interest
  select(., RECORDING_SESSION_LABEL, TRIAL_INDEX, BIN_INDEX,
         AVERAGE_IA_0_SAMPLE_COUNT, AVERAGE_IA_0_SAMPLE_COUNT_.,
         AVERAGE_IA_1_SAMPLE_COUNT, AVERAGE_IA_1_SAMPLE_COUNT_.,
         AVERAGE_IA_2_SAMPLE_COUNT, AVERAGE_IA_2_SAMPLE_COUNT_.,
         ACCURACY, RT, block, cond, 
         id, lex_freq, phonot_freq, 
         t01, t02, t03, t04, t05, t06, t07, target, version, group) %>%
  dplyr::rename(., participant = RECORDING_SESSION_LABEL,
                trial = TRIAL_INDEX, 
                bin = BIN_INDEX,
                target_count = AVERAGE_IA_1_SAMPLE_COUNT, 
                target_prop = AVERAGE_IA_1_SAMPLE_COUNT_.,
                offset_prev_word = t01,
                onset_v1 = t02,
                onset_c2 = t03,
                onset_c3 = t04,
                onset_v2 = t05,
                offset_target = t06,
                endSentence = t07,
                sentence_id = id) %>%
  
  # remove incorrect
  filter(., ACCURACY == 1) %>%
  
  # drop unused levels of factors
  droplevels(.) %>%
  
  # Create eLog variable and respective wts
  mutate(.,eLog = log((target_count + 0.5) / (10 - target_count + 0.5)),
         wts = 1 / (target_count + 0.5) + 1 / (10 - target_count + 0.5)) %>%
  
  # CHANGE onset_c3 DEPENDING ON TRIGGER TO ANALYZE
  dplyr::select(participant, group, target, cond, target, bin,
                target_count, target_prop, eLog, wts, onset_v1) %>% 
  gather(., landmark, lm_bin, -c(participant:wts)) %>%
  mutate(., lm_bin = (lm_bin / 10) %>% ceiling(.),
         t_onset = if_else(bin == lm_bin, TRUE, FALSE)) %>%

  group_by(., participant, target) %>%
  mutate(., time_zero = onset_pupil(bin, t_onset, event = c("TRUE"))) %>%
  ungroup(.)

# Load verbal WM
dem <- read_csv(here("data", "pupurri_analysis.csv"))
dem <- dem %>%
  select(., participant, DELE, percent_l2_week, WM_set)

dem$participant <- tolower(dem$participant)
dem$DELE <- as.numeric(dem$DELE)

# Add verbal wm score to eyetracking data frame
stress10 <- merge(x = stress10, y = dem, by = "participant", all.x=TRUE)

# Create L1 column
stress10 <- separate(data = stress10,
                     col = group,
                     into = c("prof", "l1"),
                     sep = 1,
                     remove = FALSE)

stress10$l1 <- str_replace(stress10$l1, "es", "en")
stress10$l1 <- str_replace(stress10$l1, "ms", "ma")
stress10$l1 <- str_replace(stress10$l1, "on", "es")

# This needs to be updated for all df at triggers other than C3
stress10$DELE[is.na(stress10$DELE) & stress10$l1 == 'es'] <- 56
stress10$percent_l2_week[is.na(stress10$percent_l2_week) & stress10$l1 == 'es'] <- 0

stress10 <- stress10 %>%
  # mutate(DELE = DELE + runif(n(), min = -0.15, max = 0.15) * (n() > 1)) %>%
  mutate(., ospan = (WM_set - mean(WM_set))/sd(WM_set),
         use_z = (percent_l2_week - mean(percent_l2_week))/sd(percent_l2_week),
         DELE_z = (DELE - mean(DELE))/sd(DELE)
  )
# 
# 
# wm <- read_csv("./data/clean/ospan_set_z_scores.csv")

# stress10 <- merge(x = stress10, y = wm, by = "participant", all.x=TRUE)

write_csv(stress10, here("data", "clean", "stress_10ms_final_onset_v1.csv"))


# -----------------------------------------------------------------------------

stress10 <- read_csv(here("data", "clean", "stress_10ms_final.csv"))

# Test plot for stress_unrelated
stress10$cond <- factor(stress10$cond, levels = c("1", "2"), 
                  labels = c("Present", "Preterit"))

timecourse10 <- stress10 %>%
  filter(time_zero > -20, time_zero < 80) %>% 
  ggplot(., aes(x = time_zero, y = target_prop, color = group)) +
  facet_grid(. ~ cond) +
  geom_vline(xintercept = 20, lty = 3) +
  geom_hline(yintercept = 0.5, color = "white", size = 3) +
  stat_summary(fun.y = mean, geom = "line") +
  ggtitle("Time course per verbal tense") +
  xlab("Time in 10 ms bins (0 = 1st syllable offset; 200 ms processing not accounted)") +
  ylab("Proportion of fixations on target") +
  scale_color_discrete(name="Group",
                      breaks = c("aes", 'ams', 'ies', 'ims', 'mon'),
                      labels = c("Adv EN", 'Adv MA', "Int EN", 'Int MA', 'ES Sp'))

ggsave('timecourse10.png',
       plot = timecourse10, dpi = 600, device = "png",
       path = './figs/stress/',
       height = 3.5, width = 8.5, units = 'in')
