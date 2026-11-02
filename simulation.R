# Load necessary libraries
library(tidyverse)
library(survival)
library(survminer)
library(caret)
library(ggrepel)
library(gridExtra)

# ==========================================
# 3. SIMULATE DETECTION PROCESS
# ==========================================
simulate_detection <- function(risk_group, event_time, has_event, sens, spec) {

  schedule <- if (risk_group == "High Risk") {
    schedule_high
  } else if (risk_group == "Mid Risk") {
    schedule_mid
  } else {
    schedule_low
  }

  # Initialize counters for this patient
  n_tp <- 0; n_fp <- 0; n_tn <- 0; n_fn <- 0
  detect_time <- max_follow_up
  status <- 0

  for (test_day in schedule) {
    # 1. Determine Disease Status at this specific day
    is_sick <- has_event && (test_day >= event_time)

    # 2. Simulate Test Result
    if (is_sick) {
      test_result <- rbinom(1, 1, sens)
      if (test_result == 1) {
        n_tp <- n_tp + 1
        detect_time <- test_day
        status <- 1
        break # Stop testing after detection (biopsy & treat)
      } else {
        n_fn <- n_fn + 1 # False Negative
      }
    } else {
      test_result <- rbinom(1, 1, 1 - spec)
      if (test_result == 1) {
        n_fp <- n_fp + 1 # False Positive
      } else {
        n_tn <- n_tn + 1 # True Negative
      }
    }
  }

  return(list(
    detect_time = detect_time,
    status = status,
    n_tp = n_tp, n_fp = n_fp, n_tn = n_tn, n_fn = n_fn
  ))
}

# ==========================================
# 1. PARAMETERS & ASSUMPTIONS
# ==========================================

dpi <- 300
fig_width <- 16
fig_height <- 12
decimals <- 1
pirche_color_two <- c("#7CC89C", "#9D4180")
pirche_color_four <- c("#7cc89c", "#a49c53", "#bb6652", "#9d4180")
pirche_gradient_3 <- c("#7ec99d", "#b5834b", "#9e4281")
pirche_gradient_11 <- c("#7ec99d","#88bd82","#95b06b","#a1a259","#ac934e","#b5834b","#bb724e","#bd6157","#b95264","#af4772","#9e4281")
pirche_grey <- "#e3e1ec"
manuscript_root <- "/Users/pir-mn/Library/CloudStorage/GoogleDrive-matthias.niemann@pirche.com/Geteilte Ablagen/PIR R&D/Projects/cfDNA Target Trial Emulation/manuscript_figures/"
supplement <- "supplement/"

# Cohort settings
iterations <- 100  # number of repeats
n_patients <- 200   # how many patients to simulate
max_follow_up <- 3 * 360 # follow-up in days

# ----------------------------------------------------
# Define Set of (Sensitivity, Specificity) Tuples Here
# ----------------------------------------------------
test_scenarios <- tibble::tribble(
  ~scenario_id,             ~test_sensitivity, ~test_specificity,
  "Bloom_2017",             0.59,              0.85,
  "Sigdel_2019",            0.89,              0.73,
  "Whitlam_2019_1",         0.85,              0.75,
  "Whitlam_2019_2",         0.85,              0.79,
  "Oellerich_2019_1",       0.73,              0.69,
  "Oellerich_2019_2",       0.73,              0.73,
  "Huang_2019",             0.79,              0.72,
  "Dauber_2020",            0.88,              0.81,
  "Zhang_2020",             0.89,              0.74,
  "De_Vlaminck_2014",       0.58,              0.93,
  "Richmond_2020",          0.60,              0.92,
  "Khush_2019",             0.44,              0.80,
  "Agbor_Enoh_2021",        0.81,              0.85,
  "Knuettgen_2022",         0.76,              0.83,
  "Kim_2022",               0.79,              0.77,
  "Sorbini_2021",           0.64,              0.71,
  "De_Vlaminck_2015",       1.00,              0.73,
  "Sayah_2020",             0.73,              0.53,
  "Khush_2021",             0.56,              0.76,
  "Jang_2021",              0.77,              0.84,
  "Keller_2022",            0.76,              0.70,
  "Rosenheck_2022",         0.89,              0.83,
  "Schuetz_2017",           0.90,              0.93,
  "Levitsky_2021",          1.00,              0.80,
  "Baumann_2022",           0.81,              0.91,
  "Ventura_Aguiar_2022_1",  0.29,              1.00,
  "Ventura_Aguiar_2022_2",  0.86,              0.94
)

# Study parameters (sensitivity and specificity removed from here)
utah <- data.frame(
  name = "Utah",
  prop_high_risk = 0.57,
  hr_high_risk = 1.4,
  baseline_incidence_1yr = 0.323
)

leuven <- data.frame(
  name = "Leuven",
  prop_high_risk = 0.33,
  hr_high_risk = 1.28,
  baseline_incidence_1yr = 0.097
)

basel <- data.frame(
  name = "Basel",
  prop_high_risk = 0.27,
  hr_high_risk = 2.19,
  baseline_incidence_1yr = 0.107
)

simulation_parameters <- list(utah, leuven, basel)

cost_pirche <- 1500
cost_cfdna <- 2800 + 200 # cfDNA test + dsa test
mode <- "low"

schedule_high <- 30 * c(1, 2, 3, 4, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36)
schedule_low  <- 30 * c(3, 6, 9, 12, 24, 36)

print(paste("Standard Schedule:", paste(schedule_low, collapse = ", ")))
print(paste("High Schedule:", paste(schedule_high, collapse = ", ")))

### Plot frequency
plot_data <- data.frame(
  x_value = c(schedule_high, schedule_low),
  y_index = c(rep(1, length(schedule_high)), rep(2, length(schedule_low))),
  group_id = c(rep("EF", length(schedule_high)), rep("SF", length(schedule_low)))
)
plot_data$group_id <- factor(plot_data$group_id, levels = c("SF", "EF"))

schedule <- ggplot(plot_data, aes(x = x_value, y = y_index, color = group_id, shape = group_id)) +
  geom_point(size = 4, alpha = 0.8) +
  scale_color_manual(values = c("EF" = pirche_color_four[3], "SF" = pirche_color_four[1])) +
  scale_shape_manual(values = c("EF" = 17, "SF" = 15)) +
  ylim(0, 3) +
  scale_x_continuous(breaks = seq(0, max(plot_data$x_value), by = 180)) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  labs(
    title = "Frequency of immune monitoring",
    x = "Days post Tx",
    y = "Index",
    color = "Schedule",
    shape = "Schedule"
  )

ggsave(paste0(manuscript_root, "schedule.png"), plot = schedule, dpi = dpi, width = fig_width, height = 0.4 * fig_height, units = "cm", limitsize = FALSE)

# ====================================================
# 2. RUN SIMULATIONS ACROSS STUDIES & TEST TUPLES
# ====================================================
global_result <- data.frame()

for (index in seq_along(simulation_parameters)) {
  parameters <- simulation_parameters[[index]]

  for (scen_idx in seq_len(nrow(test_scenarios))) {
    current_scenario <- test_scenarios[scen_idx, ]
    sens_val <- current_scenario$test_sensitivity
    spec_val <- current_scenario$test_specificity
    scen_id  <- current_scenario$scenario_id

    cat(sprintf("Simulating Study: %s | Scenario: %s (Sens: %.2f, Spec: %.2f)\n",
                parameters$name, scen_id, sens_val, spec_val))

    joint_result <- data.frame()

    for (iteration in seq(iterations)) {
      set.seed(iteration * 42)

      cohort <- tibble(
        id = 1:n_patients,
        risk_group = sample(c("Low Risk", "High Risk"), n_patients, replace = TRUE,
                            prob = c(1 - parameters$prop_high_risk, parameters$prop_high_risk))
      )

      shape_param <- 0.9
      t_1yr <- 365
      scale_low <- t_1yr / (-log(1 - parameters$baseline_incidence_1yr))^(1 / shape_param)

      cohort <- cohort %>%
        mutate(
          scale = ifelse(risk_group == "High Risk",
                         scale_low / (parameters$hr_high_risk^(1 / shape_param)),
                         scale_low),
          true_event_time = rweibull(n(), shape = shape_param, scale = scale),
          has_event = true_event_time <= max_follow_up,
          true_event_time = pmin(true_event_time, max_follow_up)
        )

      # Simulate using current tuple's sens_val & spec_val
      results <- cohort %>%
        rowwise() %>%
        mutate(
          sim_res_mixed = list(simulate_detection(risk_group, true_event_time, has_event, sens_val, spec_val)),
          detect_time_mixed = sim_res_mixed$detect_time,
          status_mixed = sim_res_mixed$status,
          n_tp_mixed = sim_res_mixed$n_tp, n_fp_mixed = sim_res_mixed$n_fp, n_tn_mixed = sim_res_mixed$n_tn, n_fn_mixed = sim_res_mixed$n_fn,
          detection_delay_mixed = ifelse(status_mixed == 1, detect_time_mixed - true_event_time, NA)
        ) %>%
        ungroup() %>%
        rowwise() %>%
        mutate(
          sim_res_low = list(simulate_detection("Low Risk", true_event_time, has_event, sens_val, spec_val)),
          detect_time_low = sim_res_low$detect_time,
          status_low = sim_res_low$status,
          n_tp_low = sim_res_low$n_tp, n_fp_low = sim_res_low$n_fp, n_tn_low = sim_res_low$n_tn, n_fn_low = sim_res_low$n_fn,
          detection_delay_low = ifelse(status_low == 1, detect_time_low - true_event_time, NA)
        ) %>%
        ungroup() %>%
        rowwise() %>%
        mutate(
          sim_res_high = list(simulate_detection("High Risk", true_event_time, has_event, sens_val, spec_val)),
          detect_time_high = sim_res_high$detect_time,
          status_high = sim_res_high$status,
          n_tp_high = sim_res_high$n_tp, n_fp_high = sim_res_high$n_fp, n_tn_high = sim_res_high$n_tn, n_fn_high = sim_res_high$n_fn,
          detection_delay_high = ifelse(status_high == 1, detect_time_high - true_event_time, NA)
        ) %>%
        ungroup()

      delay_mixed <- results %>%
        filter(status_mixed == 1) %>%
        summarise(median_delay = median(detection_delay_mixed, na.rm = TRUE))

      delay_high <- results %>%
        filter(status_high == 1) %>%
        summarise(median_delay = median(detection_delay_high, na.rm = TRUE))

      delay_low <- results %>%
        filter(status_low == 1) %>%
        summarise(median_delay = median(detection_delay_low, na.rm = TRUE))

      tests_performed <- results %>%
        rowwise() %>%
        mutate(
          schedule_mixed = list(if (risk_group == "High Risk") schedule_high else schedule_low),
          n_tests_mixed = sum(unlist(schedule_mixed) > 0),
          schedule_high = list(schedule_high),
          n_tests_high = sum(unlist(schedule_high) > 0),
          schedule_low = list(schedule_low),
          n_tests_low = sum(unlist(schedule_low) > 0)
        )

      tests_performed_strategy_summary <- tests_performed %>%
        ungroup() %>%
        summarise(
          total_tests_mixed = sum(n_tests_mixed),
          total_tests_high  = sum(n_tests_high),
          total_tests_low   = sum(n_tests_low)
        )

      result <- data.frame(
        strategy = c("EF", "SF", "PF"),
        costs = c(
          tests_performed_strategy_summary$total_tests_high * cost_cfdna / 1e6,
          tests_performed_strategy_summary$total_tests_low  * cost_cfdna / 1e6,
          (tests_performed_strategy_summary$total_tests_mixed * cost_cfdna + nrow(cohort) * cost_pirche) / 1e6
        ),
        tests = c(
          tests_performed_strategy_summary$total_tests_high,
          tests_performed_strategy_summary$total_tests_low,
          tests_performed_strategy_summary$total_tests_mixed
        ),
        median_delay = c(
          delay_high$median_delay,
          delay_low$median_delay,
          delay_mixed$median_delay
        ),
        iteration = iteration,
        simulation = parameters$name,
        scenario = scen_id,
        sensitivity = sens_val,
        specificity = spec_val
      )
      result$delay_per_cost <- result$median_delay / result$costs

      joint_result <- bind_rows(joint_result, result)
      global_result <- bind_rows(global_result, result)
    }

    # Cost-effectiveness plot for each study & scenario
    labeled_points <- joint_result %>%
      group_by(strategy) %>%
      mutate(
        target_x = median(median_delay, na.rm = TRUE),
        target_y = median(costs, na.rm = TRUE),
        dist_to_center = sqrt((median_delay - target_x)^2 + (costs - target_y)^2)
      ) %>%
      filter(dist_to_center == min(dist_to_center)) %>%
      slice(1)

    cost_effectiveness <- ggplot(joint_result, aes(x = median_delay, y = costs, color = strategy)) +
      geom_point(size = 1, alpha = 0.25) +
      geom_text(
        data = labeled_points,
        aes(label = paste0(strategy, "\n(", sprintf(paste0("%.", decimals, "f"), target_y), "M $, ", sprintf(paste0("%.", decimals, "f"), target_x), "d)"),
            x = target_x, y = target_y),
        vjust = -0.5, fontface = "bold", lineheight = 0.8, size = 4, show.legend = FALSE
      ) +
      geom_point(data = labeled_points, size = 4, shape = 3, color = "black") +
      scale_color_manual(values = c("EF" = pirche_color_four[3], "SF" = pirche_color_four[1], "PF" = pirche_color_four[4])) +
      theme_light() +
      labs(
        title = paste0(parameters$name, " (", scen_id, ")\nCost effectiveness (n = ", n_patients, ")"),
        x = "median days to detection (lower is better)",
        y = "costs in millions (lower is better)",
        color = "Schedule"
      )

    ggsave(paste0(manuscript_root, "sim_effectiveness/cost_effectiveness_", parameters$name, "_", scen_id, ".png"),
           plot = cost_effectiveness, dpi = dpi, width = fig_width, height = fig_height, units = "cm", limitsize = FALSE)
  }
}

global_result$strategy <- factor(global_result$strategy, levels = c("SF", "EF", "PF"))
global_result$simulation <- as.factor(global_result$simulation)
global_result$scenario <- as.factor(global_result$scenario)

# 2. Create the plot
global_plot <- ggplot(global_result, aes(x = simulation, y = median_delay, fill = strategy)) +
  geom_boxplot() +
  scale_fill_manual(values=c("SF" = pirche_gradient_11[11], "EF" = pirche_gradient_11[1], "PF" = pirche_gradient_11[6])) +
  labs(
    title = "Time-to-detection by Surveillance Schedule",
    x = "Study",
    y = "time-to-detection in days",
    fill = "Schedule"
  ) +
  theme_light()

ggsave(paste0(manuscript_root, "delay.png"), plot = global_plot, dpi = dpi, width = fig_width, height = fig_height, units = "cm", limitsize = FALSE)

# plot the publication-specific efficiency

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

# 1. Aggregate data to find medians and quartiles
summary_data <- global_result %>%
  group_by(strategy, simulation) %>%
  summarise(
    x_median = median(median_delay, na.rm = TRUE),
    x_q1     = quantile(median_delay, 0.25, na.rm = TRUE),
    x_q3     = quantile(median_delay, 0.75, na.rm = TRUE),
    y_median = median(tests, na.rm = TRUE),
    y_q1     = quantile(tests, 0.25, na.rm = TRUE),
    y_q3     = quantile(tests, 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# 2. Create the segment dataframe and drop any incomplete pairs
connection_lines <- summary_data %>%
  filter(strategy %in% c("EF", "SF")) %>%
  select(simulation, strategy, x_median, y_median) %>%
  pivot_wider(
    names_from = strategy,
    values_from = c(x_median, y_median)
  ) %>%
  # CRITICAL: Ensures rows with missing EF or SF coordinates are handled safely
  drop_na(x_median_EF, y_median_EF, x_median_SF, y_median_SF)

# 3. Build the plot
line_width <- 0.5
wisker_width_x <- 6   # Controls the cap width on the horizontal axis (geom_errorbar)
wisker_height_y <- 40   # Controls the cap height on the vertical axis (geom_errorbarh)

test_effectiveness <- ggplot(summary_data, aes(x = x_median, y = y_median)) +

  # Moving geom_segment above error bars so it draws first
  geom_segment(
    data = connection_lines,
    aes(
      x = x_median_EF, y = y_median_EF,
      xend = x_median_SF, yend = y_median_SF
    ),
    color = "darkgray",       # Darker grey for better visibility
    linetype = "dashed",
    linewidth = 0.8,
    inherit.aes = FALSE     # Completely isolates aesthetic mappings
  ) +

  # Vertical error bars
  geom_errorbar(
    aes(ymin = y_q1, ymax = y_q3, color = strategy),
    width = wisker_width_x,   # Uses X-axis unit scaling
    linewidth = line_width,
    alpha = 0.8
  ) +

  # Horizontal error bars
  geom_errorbarh(
    aes(xmin = x_q1, xmax = x_q3, color = strategy),
    height = wisker_height_y, # Uses Y-axis unit scaling
    linewidth = line_width,
    alpha = 0.8
  ) +

  # Center point representing the median
  geom_point(aes(color = strategy), size = 2, shape = 16) +

  facet_grid(. ~ simulation) +

  # Dynamic Text Labels
  geom_text(
    aes(
      label = paste0(
        strategy, "\n(",
        sprintf(paste0("%.", 0, "f"), y_median), "n, ",
        sprintf(paste0("%.", decimals, "f"), x_median), "d)"
      ),
      color = strategy
    ),
    vjust = -0.2,
    hjust = 0.5,
    lineheight = 0.8,
    size = 3,
    show.legend = FALSE
  ) +

  scale_color_manual(values = c("EF" = pirche_gradient_3[1], "SF" = pirche_gradient_3[3], "PF" = pirche_gradient_3[2])) +
  theme_light() +

  # CRITICAL FIX: Use coord_cartesian instead of xlim/ylim to prevent data clipping
  coord_cartesian(
    ylim = c(0.9 * min(global_result$tests, na.rm = TRUE), 1.2 * max(global_result$tests, na.rm = TRUE)),
    xlim = c(min(global_result$median_delay, na.rm = TRUE) - 10, max(global_result$median_delay, na.rm = TRUE) - 290)
  ) +

  labs(
    title = paste0("Simulated strategy efficiency (n = ", n_patients, ")"),
    x = "time-to-detection in days",
    y = "number of tests",
    color = "schedule"
  ) +
  theme(
    legend.position = "none"
  )

# 4. Save the plot
p_layout <- rbind(c(1))
panel <- arrangeGrob(test_effectiveness, layout_matrix = p_layout, ncol = 1)
ggsave(paste0(manuscript_root, "test_effectiveness.png"), plot = panel, dpi = dpi, width = fig_width, height = fig_height, units = "cm", limitsize = FALSE)

# extract numeric data for the manuscript

delay_comparison <- summary_data %>%
  # Filter for only the two strategies you want to compare
  filter(strategy %in% c("SF", "EF", "PF")) %>%
  # Keep only the necessary columns
  select(simulation, strategy, x_median) %>%
  # Move strategies into columns to have them side-by-side
  pivot_wider(names_from = strategy, values_from = x_median) %>%
  # Calculate the percentage (Extended as a % of Standard)
  mutate(
    percentage_ES = (EF / SF) * 100,
    reduction_pct_ES = 100 - percentage_ES,
    percentage_PS = (PF / SF) * 100,
    reduction_pct_PS = 100 - percentage_PS
  )

print(delay_comparison, width = Inf)

test_comparison <- summary_data %>%
  # Filter for the relevant strategies
  filter(strategy %in% c("SF", "EF", "PF")) %>%
  # Select identifying columns and the 'tests' column
  select(simulation, strategy, y_median) %>%
  # Pivot to wide format
  pivot_wider(names_from = strategy, values_from = y_median) %>%
  # Calculate the percentage
  mutate(
    percentage_ES = (EF / SF) * 100,
    reduction_pct_ES = 100 - percentage_ES,
    percentage_PS = (PF / SF) * 100,
    reduction_pct_PS = 100 - percentage_PS
  )

print(test_comparison, width = Inf)

my_comparisons <- list( c("SF", "PF") )

global_result_norm <- global_result %>%
  group_by(iteration, simulation, scenario) %>%
  mutate(
    # Min-Max Scaling calculated strictly within each strategy's group
    Delay_Scaled = (median_delay - min(median_delay, na.rm = TRUE)) /
      (max(median_delay, na.rm = TRUE) - min(median_delay, na.rm = TRUE)),

    Tests_Scaled = (tests - min(tests, na.rm = TRUE)) /
      (max(tests, na.rm = TRUE) - min(tests, na.rm = TRUE)),

    # Calculate composite score on locally weighted inputs
    nttd = 1 - (Delay_Scaled + Tests_Scaled) / 2
  ) %>%
  ungroup() %>% # CRITICAL: Ungroup so the subsequent sorting and Wilcox tests run across strategies

  # Reorder the Strategy factor by the new localized median descending
  mutate(strategy = fct_reorder(strategy, nttd, .fun = median, .desc = FALSE))

nttd_boxplot <- ggplot(global_result_norm, aes(x = strategy, y = nttd, color=strategy, fill = strategy)) +

  # Create the boxplots
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +

  # STATISTICAL SIGNIFICANCE (Now functions perfectly because strategy is on the x-axis)
  stat_compare_means(
    method = "wilcox.test",
    comparisons = my_comparisons,
    label = "p.signif",
    step.increase = 0.12,      # Slightly increased to avoid overlapping bars
    size = 4,
    tip.length = 0.005
  ) +

  # NEW: Separate into simulation panels side-by-side
  facet_grid(. ~ simulation) +

  # Custom colors matching your setup
  scale_fill_manual(values = c("EF" = pirche_gradient_3[1], "SF" = pirche_gradient_3[3], "PF" = pirche_gradient_3[2])) +
  scale_color_manual(values = c("EF" = pirche_gradient_3[1], "SF" = pirche_gradient_3[3], "PF" = pirche_gradient_3[2])) +

  # Styling
  theme_light() +
  labs(
    x = "strategy",
    y = "efficiency\nMNT = (1 - mean(n, ttd))"
  ) +
  scale_y_continuous(limits = c(min(global_result_norm$nttd), max(global_result_norm$nttd) + length(my_comparisons) * 0.10)) +

  theme(
    legend.position = "none"
  )
print(nttd_boxplot)
ggsave(paste0(manuscript_root, supplement, "sim_mnt.png"), plot = nttd_boxplot, dpi = dpi,  width = fig_width, height = fig_height, units = "cm", limitsize = FALSE)

panel <- arrangeGrob(test_effectiveness, nttd_boxplot, layout_matrix = rbind(c(1), c(2)), ncol = 1)
ggsave(paste0(manuscript_root, "sim_effectiveness_multipanel.png"), plot = panel, dpi = dpi, width = fig_width, height = fig_height, units = "cm", limitsize = FALSE)

# Calculate and print the median nxttd footprint per strategy
strategy_medians <- global_result_norm %>%
  group_by(strategy, simulation) %>%
  summarise(
    Median_nxttd = median(nttd, na.rm = TRUE),
    Min_nxttd    = min(nttd, na.rm = TRUE),
    Max_nxttd    = max(nttd, na.rm = TRUE),
    Sample_Size  = n(),
    .groups = "drop"
  ) %>%
  # Arrange from lowest loss (most efficient) to highest loss
  arrange(Median_nxttd)

# Print the structured summary table to the console
print(strategy_medians)