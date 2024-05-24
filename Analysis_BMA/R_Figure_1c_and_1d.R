# INFO ####
# This R script cleans the raw photometer data and adds the metadata
# for low-dose treated plants of both Aketo and BThu T. vulgare chemotypes
# Author: Silvia Eckert
# Generated: 2024-04-25

# =========================================================================== #

# This R script is sourced by R_Figure_1.R

# =========================================================================== #

# AKETO ####
## Low-dose treatment ####
# F10, F18, F100: Read data 
Aketo_50uM5AZA <- tidy_plates(input = here("Working_data/Mono_plates"),
                             pattern = "BMA_exp2_.._grp1",
                             how_many = "multiple",
                             direction = "horizontal",
                             skip_lines = 5,
                             group_IDs = "Aketo",
                             experiment_names = "Low-dose",
                             validity_method = "threshold",
                             threshold = 1,
                             treatment_labels = rep(c("F10", "F18", "F100", "Control"), each=2),
                             concentration_levels = rep(c(100,200),4)) %>% 
  mutate(Validity = if_else(Value < 0.06, "invalid", Validity)); Aketo_50uM5AZA
# F30: Read data
Mock_50uM5AZA_F30 <- tidy_plates(input = here("Working_data/Mixed_plates"),
                                pattern = "BMA_exp2_.._grp3",
                                how_many = "multiple",
                                direction = "horizontal",
                                skip_lines = 5,
                                group_IDs = "Mixed",
                                experiment_names = "Low-dose",
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
  split(., f = .$Group); Mock_50uM5AZA_F30

# F10, F18, F100: Subtract background absorption at T0
Aketo_50uM5AZA_subtracted <- subtract_T0(Aketo_50uM5AZA)
# F30: Subtract background absorption at T0
Aketo_50uM5AZA_F30_subtracted <- subtract_T0(Mock_50uM5AZA_F30$Aketo)

# F10, F18, F100: Calculate growth performance
Aketo_50uM5AZA_Fx_growth <- calculate_growth_performance(Aketo_50uM5AZA_subtracted,
                                                        concentration_grouping = TRUE)
# F30: Calculate growth performance
Aketo_50uM5AZA_F30_growth <- calculate_growth_performance(Aketo_50uM5AZA_F30_subtracted,
                                                         concentration_grouping = TRUE)

# F10, F18, F100: Summarize growth performance
Aketo_50uM5AZA_Fx_sum <- summarize_growth_performance(Aketo_50uM5AZA_Fx_growth,
                                                     grouping = c("Group", "Experiment",
                                                                  "Treatment", "Concentration",
                                                                  "Timepoint"))
# F30: Summarize growth performance
Aketo_50uM5AZA_F30_sum <- summarize_growth_performance(Aketo_50uM5AZA_F30_growth,
                                                      grouping = c("Group", "Experiment",
                                                                   "Treatment", "Concentration",
                                                                   "Timepoint"))

# F10, F18, F100: Apply sign test
Aketo_50um5AZA_Fx_growth_T3 <- Aketo_50uM5AZA_Fx_growth[Aketo_50uM5AZA_Fx_growth$Timepoint == "T3",]
Aketo_50um5AZA_Fx_sum_T3 <- Aketo_50uM5AZA_Fx_sum[Aketo_50uM5AZA_Fx_sum$Timepoint == "T3",]
Aketo_50uM5AZA_Fx_stats <- apply_sign_test(stats_data = Aketo_50um5AZA_Fx_growth_T3,
                                          summarized_data = Aketo_50um5AZA_Fx_sum_T3,
                                          grouping = c("Group", "Experiment",
                                                       "Treatment", "Concentration"))
# F30: Apply sign test
Aketo_50um5AZA_F30_growth_T3 <- Aketo_50uM5AZA_F30_growth[Aketo_50uM5AZA_F30_growth$Timepoint == "T3",]
Aketo_50um5AZA_F30_sum_T3 <- Aketo_50uM5AZA_F30_sum[Aketo_50uM5AZA_F30_sum$Timepoint == "T3",]
Aketo_50uM5AZA_F30_stats <- apply_sign_test(stats_data = Aketo_50um5AZA_F30_growth_T3,
                                           summarized_data = Aketo_50um5AZA_F30_sum_T3,
                                           grouping = c("Group", "Experiment","Treatment",
                                                        "Concentration"))

# Combine results
Aketo_50uM5AZA_growth <- rbind(Aketo_50um5AZA_Fx_growth_T3,
                              Aketo_50um5AZA_F30_growth_T3)
Aketo_50uM5AZA_sum <- rbind(Aketo_50um5AZA_Fx_sum_T3,
                           Aketo_50um5AZA_F30_sum_T3)
Aketo_50uM5AZA_stats <- rbind(Aketo_50uM5AZA_Fx_stats,
                             Aketo_50uM5AZA_F30_stats)

# Plot results
Figure_1c <- plot_growth_performance(input_data = Aketo_50uM5AZA_sum,
                                     level_unit = "ppm",
                                     stats_data = Aketo_50uM5AZA_growth,
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
        legend.position = "none",
        axis.title = element_text(size=14)) +
  annotate("text", x = 0.5, y = 140, hjust = 0,
           label = bquote("n ="~.(min(Aketo_50uM5AZA_stats$n))*"-"*.(max(Aketo_50uM5AZA_stats$n)))); Figure_1c

# =========================================================================== #

# BTHU ####
## Low-dose treatment ####
# F10, F18, F100: Read data 
BThu_50uM5AZA <- tidy_plates(input = here("Working_data/Mono_plates"),
                            pattern = "BMA_exp2_.._grp2",
                            how_many = "multiple",
                            direction = "horizontal",
                            skip_lines = 5,
                            group_IDs = "BThu",
                            experiment_names = "Low-dose",
                            validity_method = "threshold",
                            threshold = 1,
                            treatment_labels = rep(c("F10", "F18", "F100", "Control"), each=2),
                            concentration_levels = rep(c(100,200),4)) %>% 
  mutate(Validity = if_else(Value < 0.06, "invalid", Validity)); BThu_50uM5AZA


# F10, F18, F100: Subtract background absorption at T0
BThu_50uM5AZA_subtracted <- subtract_T0(BThu_50uM5AZA)
# F30: Subtract background absorption at T0
BThu_50uM5AZA_F30_subtracted <- subtract_T0(Mock_50uM5AZA_F30$BThu)

# F10, F18, F100: Calculate growth performance
BThu_50uM5AZA_Fx_growth <- calculate_growth_performance(BThu_50uM5AZA_subtracted,
                                                       concentration_grouping = TRUE)
# F30: Calculate growth performance
BThu_50uM5AZA_F30_growth <- calculate_growth_performance(BThu_50uM5AZA_F30_subtracted,
                                                        concentration_grouping = TRUE)

# F10, F18, F100: Summarize growth performance
BThu_50uM5AZA_Fx_sum <- summarize_growth_performance(BThu_50uM5AZA_Fx_growth,
                                                    grouping = c("Group", "Experiment",
                                                                 "Treatment", "Concentration",
                                                                 "Timepoint"))
# F30: Summarize growth performance
BThu_50uM5AZA_F30_sum <- summarize_growth_performance(BThu_50uM5AZA_F30_growth,
                                                     grouping = c("Group", "Experiment",
                                                                  "Treatment", "Concentration",
                                                                  "Timepoint"))

# F10, F18, F100: Apply sign test
BThu_50um5AZA_Fx_growth_T3 <- BThu_50uM5AZA_Fx_growth[BThu_50uM5AZA_Fx_growth$Timepoint == "T3",]
BThu_50um5AZA_Fx_sum_T3 <- BThu_50uM5AZA_Fx_sum[BThu_50uM5AZA_Fx_sum$Timepoint == "T3",]
BThu_50uM5AZA_Fx_stats <- apply_sign_test(stats_data = BThu_50um5AZA_Fx_growth_T3,
                                         summarized_data = BThu_50um5AZA_Fx_sum_T3,
                                         grouping = c("Group", "Experiment",
                                                      "Treatment", "Concentration"))
# F30: Apply sign test
BThu_50um5AZA_F30_growth_T3 <- BThu_50uM5AZA_F30_growth[BThu_50uM5AZA_F30_growth$Timepoint == "T3",]
BThu_50um5AZA_F30_sum_T3 <- BThu_50uM5AZA_F30_sum[BThu_50uM5AZA_F30_sum$Timepoint == "T3",]
BThu_50uM5AZA_F30_stats <- apply_sign_test(stats_data = BThu_50um5AZA_F30_growth_T3,
                                          summarized_data = BThu_50um5AZA_F30_sum_T3,
                                          grouping = c("Group", "Experiment","Treatment",
                                                       "Concentration"))

# Combine results
BThu_50uM5AZA_growth <- rbind(BThu_50um5AZA_Fx_growth_T3,
                             BThu_50um5AZA_F30_growth_T3)
BThu_50uM5AZA_sum <- rbind(BThu_50um5AZA_Fx_sum_T3,
                          BThu_50um5AZA_F30_sum_T3)
BThu_50uM5AZA_stats <- rbind(BThu_50uM5AZA_Fx_stats,
                            BThu_50uM5AZA_F30_stats)

# Plot results
Figure_1d <- plot_growth_performance(input_data = BThu_50uM5AZA_sum,
                                     level_unit = "ppm",
                                     stats_data = BThu_50uM5AZA_growth,
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
  labs(x = NULL, y = NULL) +
  theme(axis.title = element_text(size=14)) +
  annotate("text", x = 0.5, y = 140, hjust = 0,
           label = bquote("n ="~.(min(BThu_50uM5AZA_stats$n))*"-"*.(max(BThu_50uM5AZA_stats$n)))); Figure_1d

# =========================T=H=E==E=N=D====================================== #
