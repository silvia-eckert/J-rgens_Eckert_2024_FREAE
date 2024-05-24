# INFO ####
# This R script cleans the raw photometer data and adds the metadata
# for mock-treated plants of both Aketo and BThu T. vulgare chemotypes
# Author: Silvia Eckert
# Generated: 2024-04-25

# =========================================================================== #

# This R script is sourced by R_Figure_1.R

# =========================================================================== #

# AKETO ####
## Mock treatment ####
# F10, F18, F100: Read data 
Aketo_0uM5AZA <- tidy_plates(input = here("Working_data/Mono_plates"),
                       pattern = "BMA_exp1_.._grp1",
                       how_many = "multiple",
                       direction = "horizontal",
                       skip_lines = 5,
                       group_IDs = "Aketo",
                       experiment_names = "Mock-treated",
                       validity_method = "threshold",
                       threshold = 1,
                       treatment_labels = rep(c("F10", "F18", "F100", "Control"), each=2),
                       concentration_levels = rep(c(100,200),4)) %>% 
  mutate(Validity = if_else(Value < 0.06, "invalid", Validity)); Aketo_0uM5AZA
# F30: Read data
Mock_0uM5AZA_F30 <- tidy_plates(input = here("Working_data/Mixed_plates"),
                           pattern = "BMA_exp1_.._grp3",
                           how_many = "multiple",
                           direction = "horizontal",
                           skip_lines = 5,
                           group_IDs = "Mixed",
                           experiment_names = "Mock-treated",
                           validity_method = "threshold",
                           threshold = 1,
                           treatment_labels = rep(c("Aketo_F30", "BThu_F30",
                                                    "Aketo_Control", "BThu_Control"),
                                                  each = 2),
                           concentration_levels = rep(c(100,200),4)) %>% 
  mutate(Validity = if_else(Value < 0.06, "invalid", Validity)) %>% 
  dplyr::select(!Group) %>% 
  tidyr::separate(Treatment, c("Group", "Treatment")) %>% 
  dplyr::relocate(Position, Value, Validity, Treatment, Concentration,
                  Timepoint, File, Group, Experiment) %>% 
  split(., f = .$Group); Mock_0uM5AZA_F30

# F10, F18, F100: Subtract background absorption at T0
Aketo_0uM5AZA_subtracted <- subtract_T0(Aketo_0uM5AZA)
# F30: Subtract background absorption at T0
Aketo_0uM5AZA_F30_subtracted <- subtract_T0(Mock_0uM5AZA_F30$Aketo)

# F10, F18, F100: Calculate growth performance
Aketo_0uM5AZA_Fx_growth <- calculate_growth_performance(Aketo_0uM5AZA_subtracted,
                                           concentration_grouping = TRUE)
# F30: Calculate growth performance
Aketo_0uM5AZA_F30_growth <- calculate_growth_performance(Aketo_0uM5AZA_F30_subtracted,
                                                     concentration_grouping = TRUE)

# F10, F18, F100: Summarize growth performance
Aketo_0uM5AZA_Fx_sum <- summarize_growth_performance(Aketo_0uM5AZA_Fx_growth,
                                                  grouping = c("Group", "Experiment",
                                                               "Treatment", "Concentration",
                                                               "Timepoint"))
# F30: Summarize growth performance
Aketo_0uM5AZA_F30_sum <- summarize_growth_performance(Aketo_0uM5AZA_F30_growth,
                                                  grouping = c("Group", "Experiment",
                                                               "Treatment", "Concentration",
                                                               "Timepoint"))

# F10, F18, F100: Apply sign test
Aketo_0um5AZA_Fx_growth_T3 <- Aketo_0uM5AZA_Fx_growth[Aketo_0uM5AZA_Fx_growth$Timepoint == "T3",]
Aketo_0um5AZA_Fx_sum_T3 <- Aketo_0uM5AZA_Fx_sum[Aketo_0uM5AZA_Fx_sum$Timepoint == "T3",]
Aketo_0uM5AZA_Fx_stats <- apply_sign_test(stats_data = Aketo_0um5AZA_Fx_growth_T3,
                                          summarized_data = Aketo_0um5AZA_Fx_sum_T3,
                                          grouping = c("Group", "Experiment",
                                                       "Treatment", "Concentration"))
# F30: Apply sign test
Aketo_0um5AZA_F30_growth_T3 <- Aketo_0uM5AZA_F30_growth[Aketo_0uM5AZA_F30_growth$Timepoint == "T3",]
Aketo_0um5AZA_F30_sum_T3 <- Aketo_0uM5AZA_F30_sum[Aketo_0uM5AZA_F30_sum$Timepoint == "T3",]
Aketo_0uM5AZA_F30_stats <- apply_sign_test(stats_data = Aketo_0um5AZA_F30_growth_T3,
                                       summarized_data = Aketo_0um5AZA_F30_sum_T3,
                                       grouping = c("Group", "Experiment","Treatment",
                                                    "Concentration"))

# Combine results
Aketo_0uM5AZA_growth <- rbind(Aketo_0um5AZA_Fx_growth_T3,
                              Aketo_0um5AZA_F30_growth_T3)
Aketo_0uM5AZA_sum <- rbind(Aketo_0um5AZA_Fx_sum_T3,
                              Aketo_0um5AZA_F30_sum_T3)
Aketo_0uM5AZA_stats <- rbind(Aketo_0uM5AZA_Fx_stats,
                             Aketo_0uM5AZA_F30_stats)

# Plot results
Figure_1a <- plot_growth_performance(input_data = Aketo_0uM5AZA_sum,
                                     level_unit = "ppm",
                                     stats_data = Aketo_0uM5AZA_growth,
                                     treatment_order = c("F10","F18", "F30", "F100"),
                                     apply_sign_test = T,
                                     na = "",
                                     y_lab = "Growth performance [%]",
                                     col_facets = "Group",
                                     row_facets = "Experiment",
                                     p_values = "p.signif",
                                     grouping = c("Group", "Experiment",
                                                  "Treatment", "Concentration"),
                                     level_colors = c("white","gray40")) +
  scale_y_continuous(limits = c(-100, 150)) +
  labs(x = NULL) +
  theme(strip.text.y.right = element_blank(),
        axis.title = element_text(size=14)) +
  annotate("text", x = 0.5, y = 140, hjust = 0,
           label = bquote("n ="~.(min(Aketo_0uM5AZA_stats$n))*"-"*.(max(Aketo_0uM5AZA_stats$n)))); Figure_1a

# =========================================================================== #

# BTHU ####
## Mock treatment ####
# F10, F18, F100: Read data 
BThu_0uM5AZA <- tidy_plates(input = here("Working_data/Mono_plates"),
                             pattern = "BMA_exp1_.._grp2",
                             how_many = "multiple",
                             direction = "horizontal",
                             skip_lines = 5,
                             group_IDs = "BThu",
                             experiment_names = "Mock-treated",
                             validity_method = "threshold",
                             threshold = 1,
                             treatment_labels = rep(c("F10", "F18", "F100", "Control"), each=2),
                             concentration_levels = rep(c(100,200),4)) %>% 
  mutate(Validity = if_else(Value < 0.06, "invalid", Validity)); BThu_0uM5AZA


# F10, F18, F100: Subtract background absorption at T0
BThu_0uM5AZA_subtracted <- subtract_T0(BThu_0uM5AZA)
# F30: Subtract background absorption at T0
BThu_0uM5AZA_F30_subtracted <- subtract_T0(Mock_0uM5AZA_F30$BThu)

# F10, F18, F100: Calculate growth performance
BThu_0uM5AZA_Fx_growth <- calculate_growth_performance(BThu_0uM5AZA_subtracted,
                                                        concentration_grouping = TRUE)
# F30: Calculate growth performance
BThu_0uM5AZA_F30_growth <- calculate_growth_performance(BThu_0uM5AZA_F30_subtracted,
                                                         concentration_grouping = TRUE)

# F10, F18, F100: Summarize growth performance
BThu_0uM5AZA_Fx_sum <- summarize_growth_performance(BThu_0uM5AZA_Fx_growth,
                                                     grouping = c("Group", "Experiment",
                                                                  "Treatment", "Concentration",
                                                                  "Timepoint"))
# F30: Summarize growth performance
BThu_0uM5AZA_F30_sum <- summarize_growth_performance(BThu_0uM5AZA_F30_growth,
                                                      grouping = c("Group", "Experiment",
                                                                   "Treatment", "Concentration",
                                                                   "Timepoint"))

# F10, F18, F100: Apply sign test
BThu_0um5AZA_Fx_growth_T3 <- BThu_0uM5AZA_Fx_growth[BThu_0uM5AZA_Fx_growth$Timepoint == "T3",]
BThu_0um5AZA_Fx_sum_T3 <- BThu_0uM5AZA_Fx_sum[BThu_0uM5AZA_Fx_sum$Timepoint == "T3",]
BThu_0uM5AZA_Fx_stats <- apply_sign_test(stats_data = BThu_0um5AZA_Fx_growth_T3,
                                          summarized_data = BThu_0um5AZA_Fx_sum_T3,
                                          grouping = c("Group", "Experiment",
                                                       "Treatment", "Concentration"))
# F30: Apply sign test
BThu_0um5AZA_F30_growth_T3 <- BThu_0uM5AZA_F30_growth[BThu_0uM5AZA_F30_growth$Timepoint == "T3",]
BThu_0um5AZA_F30_sum_T3 <- BThu_0uM5AZA_F30_sum[BThu_0uM5AZA_F30_sum$Timepoint == "T3",]
BThu_0uM5AZA_F30_stats <- apply_sign_test(stats_data = BThu_0um5AZA_F30_growth_T3,
                                           summarized_data = BThu_0um5AZA_F30_sum_T3,
                                           grouping = c("Group", "Experiment","Treatment",
                                                        "Concentration"))

# Combine results
BThu_0uM5AZA_growth <- rbind(BThu_0um5AZA_Fx_growth_T3,
                              BThu_0um5AZA_F30_growth_T3)
BThu_0uM5AZA_sum <- rbind(BThu_0um5AZA_Fx_sum_T3,
                           BThu_0um5AZA_F30_sum_T3)
BThu_0uM5AZA_stats <- rbind(BThu_0uM5AZA_Fx_stats,
                             BThu_0uM5AZA_F30_stats)

# Plot results
Figure_1b <- plot_growth_performance(input_data = BThu_0uM5AZA_sum,
                                     level_unit = "ppm",
                                     stats_data = BThu_0uM5AZA_growth,
                                     treatment_order = c("F10","F18", "F30", "F100"),
                                     apply_sign_test = T,
                                     na = "",
                                     y_lab = "Growth performance [%]",
                                     col_facets = "Group",
                                     row_facets = "Experiment",
                                     p_values = "p.signif",
                                     grouping = c("Group", "Experiment","Treatment", "Concentration"),
                                     level_colors = c("white","gray40")) +
  scale_y_continuous(limits = c(-100, 150)) +
  labs(x = NULL, y = NULL) +
  theme(legend.position = "none",
        axis.title = element_text(size=14)) +
  annotate("text", x = 0.5, y = 140, hjust = 0,
           label = bquote("n ="~.(min(BThu_0uM5AZA_stats$n))*"-"*.(max(BThu_0uM5AZA_stats$n)))); Figure_1b

# =========================T=H=E==E=N=D====================================== #
