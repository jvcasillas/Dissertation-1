# Time course plots -----------------------------------------------------------
#
# - This script plots the raw eye tracking time course data and the
#   model estimates from a growth curve analysis
#
# -----------------------------------------------------------------------------



# Source scripts and load models ----------------------------------------------

source(here::here("scripts", "00_load_libs.R"))
source(here::here("scripts", "02_load_data.R"))

# Get path to saved models 
gca_mods_path <- here("mods", "stress", "gca") 

# Load models as list and store full mod to global env
load(paste0(gca_mods_path, "/gca_mon_mods.Rdata"))
load(paste0(gca_mods_path, "/gca_l2_mods.Rdata"))
load(paste0(gca_mods_path, "/model_preds.Rdata"))
list2env(gca_mon_mods, globalenv())
list2env(gca_l2_mods, globalenv())
list2env(model_preds, globalenv())

# Set path for saving figs
figs_path <- here("figs", "stress", "gca")

# -----------------------------------------------------------------------------






# Plot raw data ---------------------------------------------------------------

stress50$cond <- as.factor(as.character(stress50$cond))

condition_names <- c(
  `1` = 'Present',
  `2` = 'Preterit'
)


stress_p1 <- stress50 %>%        
    #na.omit(.) %>%
    filter(., time_zero >= -10, time_zero <= 20) %>%
    mutate(., l1 = fct_relevel(l1, "es", "en", "ma")) %>%
    ggplot(., aes(x = time_zero, y = target_prop, fill = l1)) +
    facet_grid(. ~ cond, labeller = as_labeller(condition_names)) +
    geom_hline(yintercept = 0.5, color = 'white', size = 3) +
    geom_vline(xintercept = 0, color = 'grey40', lty = 3) +
    geom_vline(xintercept = 4, color = 'grey40', lty = 3) +
    stat_summary(fun.y = "mean", geom = "line", size = 1) +  
    stat_summary(fun.data = mean_cl_boot, geom = 'pointrange', size = 0.5,
                 stroke = 0.5, pch = 21) +
    scale_fill_brewer(palette = 'Set1', name = "L1",
                       labels = c("ES", "EN", "MA")) +
    scale_x_continuous(breaks = c(-10, 0, 10, 20),
                       labels = c("-500", "0", "500", "1000")) +
    labs(y = 'Proportion of target fixations',
         x = 'Time relative to target syllable offset (ms)',
         caption = "Mean +/- 95% CI") +
    annotate("text", x = 3.3, y = 0.02, label = '200ms',
             angle = 90, size = 3, hjust = 0) +
    theme_grey(base_size = 12, base_family = "Times")

l1_names <- c(
  `es` = 'Spanish speakers',
  `en` = 'English speakers',
  `ma` = 'Mandarin speakers'
)

stress_p2 <- stress50 %>%        
  #na.omit(.) %>%
  filter(., time_zero >= -10, time_zero <= 20) %>%
  mutate(., l1 = fct_relevel(l1, "es", "en", "ma")) %>%
  ggplot(., aes(x = time_zero, y = target_prop, fill = cond)) +
  facet_grid(. ~ l1, labeller = as_labeller(l1_names)) +
  geom_hline(yintercept = 0.5, color = 'white', size = 3) +
  geom_vline(xintercept = 0, color = 'grey40', lty = 3) +
  geom_vline(xintercept = 4, color = 'grey40', lty = 3) +
  stat_summary(fun.y = "mean", geom = "line", size = 1) +  
  stat_summary(fun.data = mean_cl_boot, geom = 'pointrange', size = 0.5,
               stroke = 0.5, pch = 21) +
  scale_fill_brewer(palette = 'Set1', name = "Tense",
                    labels = c("Present", "Preterit")) +
  scale_x_continuous(breaks = c(-10, 0, 10, 20),
                     labels = c("-500", "0", "500", "1000")) +
  labs(y = 'Proportion of target fixations',
       x = 'Time relative to target syllable offset (ms)',
       caption = "Mean +/- 95% CI") +
  annotate("text", x = 3.3, y = 0.02, label = '200ms',
           angle = 90, size = 3, hjust = 0) +
  theme_grey(base_size = 12, base_family = "Times")


ggsave('stress_l1.png',
       plot = stress_p1, dpi = 600, device = "png",
       path = figs_path,
       height = 3.5, width = 8.5, units = 'in')
ggsave('stress_tense.png',
       plot = stress_p1, dpi = 600, device = "png",
       path = figs_path,
       height = 3.5, width = 8.5, units = 'in')



# -----------------------------------------------------------------------------




# Plot GCA --------------------------------------------------------------------

# Monolingual Spanish speakers

stress_mon <- model_preds$fits_all_mon %>%
  mutate(condition = if_else(condition_sum == 1, "Present", "Preterit"),
         condition = fct_relevel(condition, "Present"))%>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                fill = condition)) + #, color = condition
  geom_hline(yintercept = 0, lty = 3, size = 0.4) +
  geom_vline(xintercept = 4, lty = 3, size = 0.4) +
  stat_summary(fun.y = "mean", geom = "line", size = 1) + 
  # geom_ribbon(alpha = 0.2, color = "grey", show.legend = F) +
  stat_summary(fun.data = mean_cl_boot, geom = 'ribbon',fun.args=list(conf.int=0.95),
               alpha = 0.5) +
  geom_point(aes(color = condition), size = 1.3, show.legend = F) +
  geom_point(aes(color = condition), size = 0.85, show.legend = F) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target") +
  theme_big + legend_adj #+ labs(color = "Condition")


# L2 speakers

# Proficiency




# Within group differences
stress_dele_l1_tog <- fits_all_l2_dele %>% #model_preds$
  mutate(condition = if_else(condition_sum == 1, "Present", "Preterit"),
         condition = fct_relevel(condition, "Present"),
         l1 = if_else(l1 == 'EN', 'English', 'Mandarin'),
         l1 = fct_relevel(l1, "English", "Mandarin")) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                fill = condition)) +
  facet_wrap(l1 ~ .) +
  geom_hline(yintercept = 0, lty = 3, size = 0.4) +
  geom_vline(xintercept = 4, lty = 3, size = 0.4) +
  stat_summary(fun.y = "mean", geom = "line", size = 1) + 
  # geom_ribbon(alpha = 0.2, color = "grey", show.legend = F) +
  stat_summary(fun.data = mean_cl_boot, geom = 'ribbon',fun.args=list(conf.int=0.95),
               alpha = 0.5) +
  geom_point(aes(color = DELE_z), alpha = .5, pch = 21, size = 0.85, show.legend = T) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target",
       color = 'L2 proficiency') +
  scale_fill_brewer(palette = 'Set1', name = "Tense",
                    labels = c("Present", "Preterit")) +
  theme_big + labs(color = "Proficiency") +
  theme(
    legend.position = c(0.06, 0.7),
    legend.key = element_blank(),
    legend.background = element_blank(),
    strip.background = element_blank(),
    axis.title.y = element_text(size = rel(.9), hjust = 0.95),
    axis.title.x = element_text(size = rel(.9)),
    legend.key.size = unit(0.75, 'lines'),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    plot.margin = unit(rep(2, 4), "mm"),
    panel.grid.major = element_line(colour = 'grey90', size = 0.15),
    panel.grid.minor = element_line(colour = 'grey90', size = 0.15)
  )


# stress_dele_l1_sep <- fits_all_l2_dele %>% #model_preds$
#   mutate(condition = if_else(condition_sum == 1, "Present", "Preterit"),
#          condition = fct_relevel(condition, "Present"),
#          l1 = if_else(l1 == 'EN', 'English', 'Mandarin'),
#          l1 = fct_relevel(l1, "English", "Mandarin")) %>%
#   ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin)) +
#   facet_grid(l1 ~ condition) +
#   geom_hline(yintercept = 0, lty = 3, size = 0.4) +
#   geom_vline(xintercept = 4, lty = 3, size = 0.4) +
#   stat_summary(fun.y = "mean", geom = "line", size = 1) + 
#   # geom_ribbon(alpha = 0.2, color = "grey", show.legend = F) +
#   stat_summary(fun.data = mean_cl_boot, geom = 'ribbon',fun.args=list(conf.int=0.95),
#                alpha = 0.5) +
#   geom_point(aes(color = DELE_z), alpha = .5, pch = 21, size = 0.85, show.legend = T) +
#   scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
#                      labels = c("-200", "0", "200", "400", "600")) +
#   labs(x = "Time (ms) relative to target syllable offset",
#        y = "Empirical logit of looks to target",
#        color = 'L2 proficiency') +
#   scale_fill_brewer(palette = 'Set1', name = "Tense",
#                     labels = c("Present", "Preterit")) +
#   theme_big + labs(color = "Proficiency") +
#   theme(
#     legend.position = c(0.06, 0.7),
#     legend.key = element_blank(),
#     legend.background = element_blank(),
#     strip.background = element_blank(),
#     axis.title.y = element_text(size = rel(.9), hjust = 0.95),
#     axis.title.x = element_text(size = rel(.9)),
#     legend.key.size = unit(0.75, 'lines'),
#     legend.text = element_text(size = 6),
#     legend.title = element_text(size = 7),
#     plot.margin = unit(rep(2, 4), "mm"),
#     panel.grid.major = element_line(colour = 'grey90', size = 0.15),
#     panel.grid.minor = element_line(colour = 'grey90', size = 0.15)
#   )


# Within condition differences
stress_dele_cond_tog <- fits_all_l2_dele %>% #model_preds$
  mutate(condition = if_else(condition_sum == 1, "Present", "Preterit"),
         condition = fct_relevel(condition, "Present"),
         l1 = if_else(l1 == 'EN', 'English', 'Mandarin'),
         l1 = fct_relevel(l1, "English", "Mandarin"))%>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                fill = l1, color = DELE_z)) +
  facet_wrap(condition ~ .) +
  geom_hline(yintercept = 0, lty = 3, size = 0.4) +
  geom_vline(xintercept = 4, lty = 3, size = 0.4) +
  stat_summary(fun.y = "mean", geom = "line", size = 1) +
  # geom_ribbon(alpha = 0.2, color = "grey", show.legend = F) +
  stat_summary(fun.data = mean_cl_boot, geom = 'ribbon',fun.args=list(conf.int=0.95),
               alpha = 0.5) +
  geom_point(aes(color = DELE_z), alpha = .6, pch = 21, 
             size = 0.85, show.legend = T) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target") +
  scale_fill_brewer(palette = 'Set1', name = "L1",
                    labels = c("English", "Mandarin")) +
  theme_big + labs(color = "Proficiency") + #legend_adj + 
  theme(
    legend.position = c(0.065, 0.65),
    legend.key = element_blank(),
    legend.background = element_blank(),
    strip.background = element_blank(),
    axis.title.y = element_text(size = rel(.9), hjust = 0.95),
    axis.title.x = element_text(size = rel(.9)),
    legend.key.size = unit(0.75, 'lines'),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6),
    plot.margin = unit(rep(2, 4), "mm"),
    panel.grid.major = element_line(colour = 'grey90', size = 0.15),
    panel.grid.minor = element_line(colour = 'grey90', size = 0.15)
  )

# Condition and L1 split
stress_dele_split <- fits_all_l2_dele %>% #model_preds$
  mutate(condition = if_else(condition_sum == 1, "Present", "Preterit"),
         condition = fct_relevel(condition, "Present"),
         l1 = if_else(l1 == 'EN', 'English', 'Mandarin'),
         l1 = fct_relevel(l1, "English", "Mandarin"))%>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                color = DELE_z)) +
  facet_grid(condition ~ l1) +
  geom_hline(yintercept = 0, lty = 3, size = 0.4) +
  geom_vline(xintercept = 4, lty = 3, size = 0.4) +
  stat_summary(fun.y = "mean", geom = "line", size = 1) + 
  # geom_ribbon(alpha = 0.2, color = "grey", show.legend = F) +
  stat_summary(fun.data = mean_cl_boot, geom = 'ribbon',fun.args=list(conf.int=0.95),
               alpha = 0.5) +
  geom_point(aes(color = DELE_z), size = 0.85, show.legend = F) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target") +
  theme_big + labs(color = "Proficiency") +
  theme(
    legend.position = c(0.1, 0.85),
    legend.key = element_blank(),
    legend.background = element_blank(),
    strip.background = element_blank(),
    axis.title.y = element_text(size = rel(.9), hjust = 0.95),
    axis.title.x = element_text(size = rel(.9)),
    legend.key.size = unit(0.75, 'lines'),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    plot.margin = unit(rep(2, 4), "mm"),
    panel.grid.major = element_line(colour = 'grey90', size = 0.15),
    panel.grid.minor = element_line(colour = 'grey90', size = 0.15)
  )




# L2 use

# Within group differences
stress_use_l1_tog <- fits_all_l2_use %>% #model_preds$
  mutate(condition = if_else(condition_sum == 1, "Present", "Preterit"),
         condition = fct_relevel(condition, "Present"),
         l1 = if_else(l1 == 'EN', 'English', 'Mandarin'),
         l1 = fct_relevel(l1, "English", "Mandarin")) %>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                fill = condition, color = use_z)) +
  facet_wrap(l1 ~ .) +
  geom_hline(yintercept = 0, lty = 3, size = 0.4) +
  geom_vline(xintercept = 4, lty = 3, size = 0.4) +
  stat_summary(fun.y = "mean", geom = "line", size = 1) + 
  # geom_ribbon(alpha = 0.2, color = "grey", show.legend = F) +
  stat_summary(fun.data = mean_cl_boot, geom = 'ribbon',fun.args=list(conf.int=0.95),
               alpha = 0.5) +
  geom_point(aes(color = use_z), alpha = .6, pch = 21, 
             size = 0.85, show.legend = T) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target") +
  scale_fill_brewer(palette = 'Set1', name = "Tense",
                    labels = c("Present", "Preterit")) +
  theme_big + labs(color = "Weekly L2 % use") +
  theme(
    legend.position = c(0.1, 0.68),
    legend.key = element_blank(),
    legend.background = element_blank(),
    strip.background = element_blank(),
    axis.title.y = element_text(size = rel(.9), hjust = 0.95),
    axis.title.x = element_text(size = rel(.9)),
    legend.key.size = unit(0.75, 'lines'),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    plot.margin = unit(rep(2, 4), "mm"),
    panel.grid.major = element_line(colour = 'grey90', size = 0.15),
    panel.grid.minor = element_line(colour = 'grey90', size = 0.15)
  )

# Within condition differences
stress_use_cond_tog <- fits_all_l2_use %>% #model_preds$
  mutate(condition = if_else(condition_sum == 1, "Present", "Preterit"),
         condition = fct_relevel(condition, "Present"),
         l1 = if_else(l1 == 'EN', 'English', 'Mandarin'),
         l1 = fct_relevel(l1, "English", "Mandarin"))%>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                fill = l1, color = use_z)) +
  facet_wrap(condition ~ .) +
  geom_hline(yintercept = 0, lty = 3, size = 0.4) +
  geom_vline(xintercept = 4, lty = 3, size = 0.4) +
  stat_summary(fun.y = "mean", geom = "line", size = 1) + 
  # geom_ribbon(alpha = 0.2, color = "grey", show.legend = F) +
  stat_summary(fun.data = mean_cl_boot, geom = 'ribbon',fun.args=list(conf.int=0.95),
               alpha = 0.5) +
  geom_point(aes(color = use_z), alpha = .6, pch = 21, 
             size = 0.85, show.legend = T) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target") +
  scale_fill_brewer(palette = 'Set1', name = "L1",
                    labels = c("EN", "MA")) +
  theme_big + labs(color = "Weekly L2 % use") +
  theme(
    legend.position = c(0.1, 0.68),
    legend.key = element_blank(),
    legend.background = element_blank(),
    strip.background = element_blank(),
    axis.title.y = element_text(size = rel(.9), hjust = 0.95),
    axis.title.x = element_text(size = rel(.9)),
    legend.key.size = unit(0.75, 'lines'),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    plot.margin = unit(rep(2, 4), "mm"),
    panel.grid.major = element_line(colour = 'grey90', size = 0.15),
    panel.grid.minor = element_line(colour = 'grey90', size = 0.15)
  )

# Condition and L1 split
stress_use_split <- fits_all_l2_dele %>% #model_preds$
  mutate(condition = if_else(condition_sum == 1, "Present", "Preterit"),
         condition = fct_relevel(condition, "Present"),
         l1 = if_else(l1 == 'EN', 'English', 'Mandarin'),
         l1 = fct_relevel(l1, "English", "Mandarin"))%>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                color = use_z)) +
  facet_grid(condition ~ l1) +
  geom_hline(yintercept = 0, lty = 3, size = 0.4) +
  geom_vline(xintercept = 4, lty = 3, size = 0.4) +
  stat_summary(fun.y = "mean", geom = "line", size = 1) + 
  # geom_ribbon(alpha = 0.2, color = "grey", show.legend = F) +
  stat_summary(fun.data = mean_cl_boot, geom = 'ribbon',fun.args=list(conf.int=0.95),
               alpha = 0.5) +
  geom_point(aes(color = use_z), size = 1.3, show.legend = F) +
  geom_point(aes(color = use_z), size = 0.85, show.legend = F) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target") +
  theme_big + labs(color = "Weekly L2 % use") +
  theme(
    legend.position = c(0.11, 0.85),
    legend.key = element_blank(),
    legend.background = element_blank(),
    strip.background = element_blank(),
    axis.title.y = element_text(size = rel(.9), hjust = 0.95),
    axis.title.x = element_text(size = rel(.9)),
    legend.key.size = unit(0.75, 'lines'),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    plot.margin = unit(rep(2, 4), "mm"),
    panel.grid.major = element_line(colour = 'grey90', size = 0.15),
    panel.grid.minor = element_line(colour = 'grey90', size = 0.15)
  )





# wm

# Within group differences
stress_wm_l1 <- model_preds$fits_all_l2_wm %>%
  mutate(condition = if_else(condition_sum == 1, "Present", "Preterit"),
         condition = fct_relevel(condition, "Present"),
         l1 = if_else(l1 == 'en', 'English', 'Mandarin'),
         l1 = fct_relevel(l1, "English", "Mandarin"))%>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                fill = condition, color = ospan)) +
  facet_wrap(l1 ~ .) +
  geom_hline(yintercept = 0, lty = 3, size = 0.4) +
  geom_vline(xintercept = 4, lty = 3, size = 0.4) +
  stat_summary(fun.y = "mean", geom = "line", size = 1) + 
  # geom_ribbon(alpha = 0.2, color = "grey", show.legend = F) +
  stat_summary(fun.data = mean_cl_boot, geom = 'ribbon',fun.args=list(conf.int=0.95),
               alpha = 0.5) +
  geom_point(aes(color = ospan), size = 1.3, show.legend = F) +
  geom_point(aes(color = ospan), size = 0.85, show.legend = F) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target") +
  scale_fill_brewer(palette = 'Set1', name = "Tense",
                    labels = c("Present", "Preterit")) +
  theme_big + labs(color = "Verbal WM") +
  theme(
    legend.position = c(0.1, 0.75),
    legend.key = element_blank(),
    legend.background = element_blank(),
    strip.background = element_blank(),
    axis.title.y = element_text(size = rel(.9), hjust = 0.95),
    axis.title.x = element_text(size = rel(.9)),
    legend.key.size = unit(0.75, 'lines'),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    plot.margin = unit(rep(2, 4), "mm"),
    panel.grid.major = element_line(colour = 'grey90', size = 0.15),
    panel.grid.minor = element_line(colour = 'grey90', size = 0.15)
  )

# Within condition differences
stress_wm_cond <- model_preds$fits_all_l2_wm %>%
  mutate(condition = if_else(condition_sum == 1, "Present", "Preterit"),
         condition = fct_relevel(condition, "Present"),
         l1 = if_else(l1 == 'en', 'English', 'Mandarin'),
         l1 = fct_relevel(l1, "English", "Mandarin"))%>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                fill = l1, color = ospan)) +
  facet_wrap(condition ~ .) +
  geom_hline(yintercept = 0, lty = 3, size = 0.4) +
  geom_vline(xintercept = 4, lty = 3, size = 0.4) +
  stat_summary(fun.y = "mean", geom = "line", size = 1) + 
  # geom_ribbon(alpha = 0.2, color = "grey", show.legend = F) +
  stat_summary(fun.data = mean_cl_boot, geom = 'ribbon',fun.args=list(conf.int=0.95),
               alpha = 0.5) +
  geom_point(aes(color = ospan), size = 1.3, show.legend = F) +
  geom_point(aes(color = ospan), size = 0.85, show.legend = F) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target") +
  scale_fill_brewer(palette = 'Set1', name = "L1",
                    labels = c("EN", "MA")) +
  theme_big + labs(color = "Verbal WM") +
  theme(
    legend.position = c(0.1, 0.75),
    legend.key = element_blank(),
    legend.background = element_blank(),
    strip.background = element_blank(),
    axis.title.y = element_text(size = rel(.9), hjust = 0.95),
    axis.title.x = element_text(size = rel(.9)),
    legend.key.size = unit(0.75, 'lines'),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    plot.margin = unit(rep(2, 4), "mm"),
    panel.grid.major = element_line(colour = 'grey90', size = 0.15),
    panel.grid.minor = element_line(colour = 'grey90', size = 0.15)
  )

# Condition and L1 split
stress_wm_split <- model_preds$fits_all_l2_wm %>%
  mutate(condition = if_else(condition_sum == 1, "Present", "Preterit"),
         condition = fct_relevel(condition, "Present"),
         l1 = if_else(l1 == 'en', 'English', 'Mandarin'),
         l1 = fct_relevel(l1, "English", "Mandarin"))%>%
  ggplot(., aes(x = time_zero, y = fit, ymax = ymax, ymin = ymin,
                color = ospan)) +
  facet_grid(condition ~ l1) +
  geom_hline(yintercept = 0, lty = 3, size = 0.4) +
  geom_vline(xintercept = 4, lty = 3, size = 0.4) +
  stat_summary(fun.y = "mean", geom = "line", size = 1) + 
  # geom_ribbon(alpha = 0.2, color = "grey", show.legend = F) +
  stat_summary(fun.data = mean_cl_boot, geom = 'ribbon',fun.args=list(conf.int=0.95),
               alpha = 0.5) +
  geom_point(aes(color = ospan), size = 1.3, show.legend = F) +
  geom_point(aes(color = ospan), size = 0.85, show.legend = F) +
  scale_x_continuous(breaks = c(-4, 0, 4, 8, 12),
                     labels = c("-200", "0", "200", "400", "600")) +
  labs(x = "Time (ms) relative to target syllable offset",
       y = "Empirical logit of looks to target") +
  theme_big + labs(color = "Verbal WM") +
  theme(
    legend.position = c(0.1, 0.85),
    legend.key = element_blank(),
    legend.background = element_blank(),
    strip.background = element_blank(),
    axis.title.y = element_text(size = rel(.9), hjust = 0.95),
    axis.title.x = element_text(size = rel(.9)),
    legend.key.size = unit(0.75, 'lines'),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 7),
    plot.margin = unit(rep(2, 4), "mm"),
    panel.grid.major = element_line(colour = 'grey90', size = 0.15),
    panel.grid.minor = element_line(colour = 'grey90', size = 0.15)
  )





ggsave(paste0(figs_path, "/stress_mon.png"), stress_mon, width = 150,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/stress_dele_l1_tog.png"), stress_dele_l1_tog, width = 150,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/stress_dele_cond_tog.png"), stress_dele_cond_tog, width = 150,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/stress_dele_split.png"), stress_dele_split, width = 150,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/stress_use_l1_tog.png"), stress_use_l1_tog, width = 150,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/stress_use_cond_tog.png"), stress_use_cond_tog, width = 150,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/stress_use_split.png"), stress_use_split, width = 150,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/stress_wm_l1.png"), stress_wm_l1, width = 150,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/stress_wm_cond.png"), stress_wm_cond, width = 150,
       height = 120, units = "mm", dpi = 600)
ggsave(paste0(figs_path, "/stress_wm_split.png"), stress_wm_split, width = 150,
       height = 120, units = "mm", dpi = 600)

# -----------------------------------------------------------------------------
