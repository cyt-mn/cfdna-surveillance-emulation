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
  } else {
    if (risk_group == "Mid Risk") {
      schedule_mid
    } else {
      schedule_low
    }
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
      # Patient has rejection. Test is Positive with prob = Sensitivity
      test_result <- rbinom(1, 1, sens)
      if (test_result == 1) {
        n_tp <- n_tp + 1
        detect_time <- test_day
        status <- 1
        break # Stop testing after detection (biopsy & treat)
      } else {
        n_fn <- n_fn + 1 # False Negative (Missed it!)
      }
    } else {
      # Patient is stable. Test is Positive with prob = (1 - Specificity)
      # This is a False Positive (False Alarm)
      test_result <- rbinom(1, 1, 1 - spec)
      if (test_result == 1) {
        n_fp <- n_fp + 1
        # Assumption: False Positive leads to biopsy, which is negative.
        # Monitoring continues? Or stops?
        # Usually monitoring continues, so we DO NOT break here.
      } else {
        n_tn <- n_tn + 1 # True Negative (Reassured)
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
iterations <- 1000  # number of repeats
n_patients <- 200    # how many patients to simulate
max_follow_up <- 3*360 # follow-up in days

first_attempt <- data.frame(
  prop_high_risk=0.5, # 50% of patients are High PIRCHE
  hr_high_risk=1.5,
  baseline_incidence_1yr=0.015, # 1.5% ABMR rate in Low Risk at 1 year
  test_sensitivity=0.83, # dd-cfDNA + dnDSA Test Performance
  test_specificity=0.79 # dd-cfDNA + dnDSA Test Performance
)

utah <- data.frame(
  prop_high_risk=0.57,
  hr_high_risk=1.4,
  baseline_incidence_1yr=0.323,
  test_sensitivity=0.80, # dd-cfDNA + dnDSA Test Performance
  test_specificity=0.81 # dd-cfDNA + dnDSA Test Performance
)

leuven <- data.frame(
  prop_high_risk=0.33,
  hr_high_risk=1.28,
  baseline_incidence_1yr=0.097,
  test_sensitivity=0.80, # dd-cfDNA + dnDSA Test Performance
  test_specificity=0.81 # dd-cfDNA + dnDSA Test Performance
)

basel <- data.frame(
  prop_high_risk=0.27,
  hr_high_risk=2.19,
  baseline_incidence_1yr=0.107,
  test_sensitivity=0.80, # dd-cfDNA + dnDSA Test Performance
  test_specificity=0.81 # dd-cfDNA + dnDSA Test Performance
)

# Table with Parameters:
# Utah
# High-risk percentage: 57% (Table 1)
# Hazard: 1.40 (Table 1, any rejection, 3y)
# Baseline risk: 32.3% (Table 1, any rejection with PIRCHE-T2 low) @ 1 year
# Leuven
# High-risk percentage: 33% (Figure 3)
# Hazard: 1.28 (37% /  ((27%+31%) / 2), Figure 3, high risk / low+int, @ 3 years)
# Baseline risk: 23.5% (21%+26%) / 2 (Figure 3, low+int 3y incidence) @ 1 year
# Does not take non-linearity into account, so attempt with 1/3rd of 3 years values
# Baseline risk: 9.7% (27%+31%) / 2 / 3 (Figure 3, low+int 3y incidence) @ 1 year
# Basel
# High-risk percentage: 27% = 94/(262+83) (Table 1)
# Hazard: 2.19 (23.4% / 10.7%, Table 2)
# Baseline risk: 10.7% = 6% x (83/(83+262)) + 12.2% * (262/(83+262))(Table 2) @ 1 year


simulation_parameters <- list(utah, leuven, basel)


cost_pirche <- 1500
cost_cfdna <- 2800 + 200 # cfDNA test + dsa test

mode <- "low"

########################## less frequent simulation
schedule_high <- 30* c(1, 2, 3, 4, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36)
# schedule_high <- 30* c(1, 2, 3, 4, 6, 12, 18, 24, 30, 36)
# schedule_high <- 30 * c(seq(from = 1, to = 35, by = 2), 36)

# Low Risk: Standard (e.g., Months 1, 3, 6, 12, then yearly)
schedule_low <- 30 * c(3, 6, 9, 12, 24, 36)

print(paste("Standard Schedule:", paste(schedule_low, collapse=", ")))
print(paste("High Schedule:", paste(schedule_high, collapse=", ")))

### plot frequency
plot_data <- data.frame(
  x_value = c(schedule_high, schedule_low),
  y_index = c(rep(1, length(schedule_high)), rep(2, length(schedule_low))),
  group_id = c(rep("EF", length(schedule_high)), rep("SF", length(schedule_low)))
)

plot_data$group_id <- factor(plot_data$group_id, levels = c("SF", "EF"))

# 3. Create the ggplot
schedule <- ggplot(plot_data, aes(x = x_value, y = y_index, color = group_id, shape = group_id)) +

  # Add points (increase size for visibility)
  geom_point(size = 4, alpha = 0.8) +

  # Optional: Manually choose your colors and symbols
  scale_color_manual(values = c("EF" = pirche_color_four[3], "SF" = pirche_color_four[1])) +
  scale_shape_manual(values = c("EF" = 17, "SF" = 15)) +
  ylim(0,3) +
  scale_x_continuous(breaks = seq(0, max(plot_data$x_value), by = 180)) +
  # Clean up the theme
  theme_minimal() +
  theme(
    axis.title.y = element_blank(), # Removes the label "y"
    axis.text.y  = element_blank(), # Removes the numbers (1, 2, 3...)
    axis.ticks.y = element_blank()  # Removes the little tick marks
  ) +
  labs(
    title = "Frequency of immune monitoring",
    x = "Days post Tx",
    y = "Index",
    color = "Schedule",
    shape = "Schedule"
  )

ggsave(paste0(manuscript_root, "schedule.png"), plot = schedule, dpi = dpi, width = fig_width, height = 0.4*fig_height, units = "cm", limitsize = FALSE)


global_result <- data.frame(strategy = character(), costs=numeric(), median_delay=numeric(), iteration=numeric(), simulation=numeric())

for (index in seq(length(simulation_parameters))) {
  # ==========================================
  # 2. GENERATE PATIENT COHORT & EVENTS
  # ==========================================
  parameters <- simulation_parameters[[index]]
  print(parameters)
  joint_result <- data.frame(strategy = character(), costs=numeric(), median_delay=numeric(), iteration=numeric(), simulation=numeric())

  for (iteration in seq(iterations)) {
    set.seed(iteration * 42)
    # set.seed(42) # For reproducibility
    # Create Dataframe
    cohort <- tibble(
      id = 1:n_patients,
      risk_group = sample(c("Low Risk", "High Risk"), n_patients, replace = TRUE,
                          prob = c(1 - parameters$prop_high_risk, parameters$prop_high_risk))
    )

    # Simulate (unknown) "True" Time to ABMR (Weibull Distribution)
    # We solve for the scale parameter (lambda) so that the 1-year incidence
    # matches our assumption for the Low Risk group.
    shape_param <- 0.9 # Shape < 1 implies decreasing hazard (early risk is highest)
    t_1yr <- 365
    scale_low <- t_1yr / (-log(1 - parameters$baseline_incidence_1yr))^(1/shape_param)

    # Generate event times
    cohort <- cohort %>%
      mutate(
        # Apply Hazard Ratio: Scale parameter decreases as Risk increases
        scale = ifelse(risk_group == "High Risk",
                       scale_low / (parameters$hr_high_risk^(1/shape_param)),
                       scale_low),

        # Random draw from Weibull
        true_event_time = rweibull(n(), shape = shape_param, scale = scale),

        # Censor events happening after follow-up time
        has_event = true_event_time <= max_follow_up,
        true_event_time = pmin(true_event_time, max_follow_up)
      )

    # Apply simulation row-wise
    results <- cohort %>%
      rowwise() %>%
      mutate(
        sim_res_mixed = list(simulate_detection(risk_group, true_event_time, has_event, parameters$test_sensitivity, parameters$test_specificity)),
        detect_time_mixed = sim_res_mixed$detect_time,
        status_mixed = sim_res_mixed$status, # 1 = Detected, 0 = Censored/Missed
        n_tp_mixed = sim_res_mixed$n_tp, n_fp_mixed = sim_res_mixed$n_fp, n_tn_mixed = sim_res_mixed$n_tn, n_fn_mixed = sim_res_mixed$n_fn,
        # Calculate "Delay" (Time detected - Time onset)
        # Only relevant for detected patients
        detection_delay_mixed = ifelse(status_mixed == 1, detect_time_mixed - true_event_time, NA)
      ) %>%
      ungroup() %>%
      rowwise() %>%
      mutate(
        sim_res_low = list(simulate_detection("Low Risk", true_event_time, has_event, parameters$test_sensitivity, parameters$test_specificity)),
        detect_time_low = sim_res_low$detect_time,
        status_low = sim_res_low$status, # 1 = Detected, 0 = Censored/Missed
        n_tp_low = sim_res_low$n_tp, n_fp_low = sim_res_low$n_fp, n_tn_low = sim_res_low$n_tn, n_fn_low = sim_res_low$n_fn,
        # Calculate "Delay" (Time detected - Time onset)
        # Only relevant for detected patients
        detection_delay_low = ifelse(status_low == 1, detect_time_low - true_event_time, NA)
      ) %>%
      ungroup() %>%
      rowwise() %>%
      mutate(
        sim_res_high = list(simulate_detection("High Risk", true_event_time, has_event, parameters$test_sensitivity, parameters$test_specificity)),
        detect_time_high = sim_res_high$detect_time,
        status_high = sim_res_high$status, # 1 = Detected, 0 = Censored/Missed
        n_tp_high = sim_res_high$n_tp, n_fp_high = sim_res_high$n_fp, n_tn_high = sim_res_high$n_tn, n_fn_high = sim_res_high$n_fn,
        # Calculate "Delay" (Time detected - Time onset)
        # Only relevant for detected patients
        detection_delay_high = ifelse(status_high == 1, detect_time_high - true_event_time, NA)
      ) %>%
      ungroup()

    # ==========================================
    # 4. ANALYSIS & VISUALIZATION
    # ==========================================
    # A. Cumulative Incidence of Detection (Kaplan-Meier)
    fit <- survfit(Surv(true_event_time, status_mixed) ~ risk_group, data = results)

    survplot <- ggsurvplot(
      fit,
      data = results,
      fun = "event", # Plot cumulative events (detections) rather than survival
      risk.table = TRUE,
      pval = TRUE,
      title = "Time to ABMR stratified by PIRCHE",
      xlab = "Days Post-Transplant",
      ylab = "ABMR incidence",
      palette = c("red", "blue"),
      legend.labs = c("High Risk", "Low Risk")
    )

    # B. Analysis of Detection Delay (for those detected)
    delay_summary_mixed <- results %>%
      filter(status_mixed == 1) %>%
      group_by(risk_group) %>%
      summarise(
        n_detected = n(),
        median_delay_days_mixed = median(detection_delay_mixed),
        mean_delay_days_mixed = mean(detection_delay_mixed),
        sd_delay_mixed = sd(detection_delay_mixed)
      )
    delay_mixed <- results %>%
      filter(status_mixed == 1) %>%
      summarise(
        n_detected = n(),
        median_delay_days_mixed = median(detection_delay_mixed),
        mean_delay_days_mixed = mean(detection_delay_mixed),
        sd_delay_mixed = sd(detection_delay_mixed)
      )

    delay_summary_high <- results %>%
      filter(status_high == 1) %>%
      group_by(risk_group) %>%
      summarise(
        n_detected = n(),
        median_delay_days_high = median(detection_delay_high),
        mean_delay_days_high = mean(detection_delay_high),
        sd_delay_high = sd(detection_delay_high)
      )
    delay_high <- results %>%
      filter(status_high == 1) %>%
      summarise(
        n_detected = n(),
        median_delay_days_high = median(detection_delay_high),
        mean_delay_days_high = mean(detection_delay_high),
        sd_delay_high = sd(detection_delay_high)
      )

    delay_summary_low <- results %>%
      filter(status_low == 1) %>%
      group_by(risk_group) %>%
      summarise(
        n_detected = n(),
        median_delay_days_low = median(detection_delay_low),
        mean_delay_days_low = mean(detection_delay_low),
        sd_delay_low = sd(detection_delay_low)
      )
    delay_low <- results %>%
      filter(status_low == 1) %>%
      summarise(
        n_detected = n(),
        median_delay_days_low = median(detection_delay_low),
        mean_delay_days_low = mean(detection_delay_low),
        sd_delay_low = sd(detection_delay_low)
      )

    # print("Summary of Diagnostic Delay (Days between Onset and Test Positive) mixed:")
    # print(delay_summary_mixed)
    # print(delay_mixed)
    # print("Summary of Diagnostic Delay (Days between Onset and Test Positive) high:")
    # print(delay_summary_high)
    # print(delay_high)
    # print("Summary of Diagnostic Delay (Days between Onset and Test Positive) low:")
    # print(delay_summary_low)
    # print(delay_low)
    # print("Summary of Diagnostic Delay (Days between Onset and Test Positive) mid:")
    # print(delay_summary_mid)
    # print(delay_mid)

    # C. Efficiency Metric (Tests per Detection)
    # Calculate total tests performed
    tests_performed <- results %>%
      rowwise() %>%
      mutate(
        schedule_mixed = list(if (risk_group == "High Risk") schedule_high else schedule_low),
        # Count tests performed up to detection or censoring
        # n_tests_mixed = sum(unlist(schedule_mixed) <= detect_time_mixed), # use this if only count tests until positivity
        n_tests_mixed = sum(unlist(schedule_mixed) > 0),
        schedule_high = list(schedule_high),
        # n_tests_high = sum(unlist(schedule_high) <= detect_time_high),
        n_tests_high = sum(unlist(schedule_high) > 0),  # use this if only count tests until positivity
        schedule_low = list(schedule_low),
        # n_tests_low = sum(unlist(schedule_low) <= detect_time_low)  # use this if only count tests until positivity
        n_tests_low = sum(unlist(schedule_low) > 0)
      )

    tests_performed_strategy_summary <- tests_performed %>%
      group_by(1) %>%
      summarise(
        total_tests_mixed = sum(n_tests_mixed),
        total_detections_mixed = sum(status_mixed),
        tests_per_detection_mixed = total_tests_mixed / total_detections_mixed,
        total_tests_high = sum(n_tests_high),
        total_detections_high = sum(status_high),
        tests_per_detection_high = total_tests_high / total_detections_high,
        total_tests_low = sum(n_tests_low),
        total_detections_low = sum(status_low),
        tests_per_detection_low = total_tests_low / total_detections_low
      )

    tests_performed_mixed_summary <- tests_performed  %>%
      group_by(risk_group) %>%
      summarise(
        total_tests_mixed = sum(n_tests_mixed),
        total_detections_mixed = sum(status_mixed),
        tests_per_detection = total_tests_mixed / total_detections_mixed
      )

    # print("Efficiency Analysis:")
    # print(tests_performed_mixed_summary)
    # print(tests_performed_strategy_summary)

    result <- data.frame(strategy = c("EF", "SF", "PF"),
                         costs=c(tests_performed_strategy_summary$total_tests_high * cost_cfdna / 1000000, tests_performed_strategy_summary$total_tests_low * cost_cfdna / 1000000, (tests_performed_strategy_summary$total_tests_mixed * cost_cfdna + nrow(cohort) * cost_pirche) / 1000000),
                         tests=c(tests_performed_strategy_summary$total_tests_high, tests_performed_strategy_summary$total_tests_low, tests_performed_strategy_summary$total_tests_mixed),
                          median_delay=c(delay_high$median_delay_days_high, delay_low$median_delay_days_low, delay_mixed$median_delay_days_mixed))
    result$delay_per_cost <- result$median_delay / result$costs
    result$iteration <- iteration
    result$simulation <- index

    joint_result <- rbind(joint_result, result)
    global_result <- rbind(global_result, result)
  }


  # 3. Create the ggplot

  labeled_points <- joint_result %>%
    group_by(strategy) %>%
    mutate(
      # Calculate the target center for this group
      target_x = median(median_delay),
      target_y = median(costs),

      # Calculate Euclidean distance from every point to that center
      dist_to_center = sqrt((median_delay - target_x)^2 + (costs - target_y)^2)
    ) %>%
    # Keep only the point with the smallest distance
    filter(dist_to_center == min(dist_to_center)) %>%
    slice(1) # Handle ties

  cost_effectiveness <- ggplot(joint_result, aes(x = median_delay, y = costs, color = strategy)) +
    # Add points (increase size for visibility)
    geom_point(size = 1, alpha = 0.25) +
    geom_text(
      data = labeled_points,
      aes(label = paste0(strategy, "\n(", sprintf(paste0("%.", decimals, "f"), target_y), "M $, ", sprintf(paste0("%.", decimals, "f"), target_x), "d)"), x = target_x, y = target_y + 0.0),
      vjust = -0.5,
      fontface = "bold",
      lineheight = 0.8,
      size=4,
      show.legend = FALSE
    ) +
    geom_point(data = labeled_points, size = 4, shape = 3, color="black") +
    # Optional: Manually choose your colors and symbols
    scale_color_manual(values = c("EF" = pirche_color_four[3], "SF" = pirche_color_four[1], "PF" = pirche_color_four[4])) +
    # Clean up the theme
    theme_light() +
    ylim(0.9*min(joint_result$costs), 1.1*max(joint_result$costs)) +
    xlim(0.8*min(joint_result$median_delay), max(joint_result$median_delay)) +
    labs(
      title = paste0("Cost effectiveness of schedules\nin detecting ABMR (n = ", n_patients, ")"),
      x = "median days to detection (lower is better)",
      y = "costs in millions (lower is better)",
      color = "Schedule"
    )
  ggsave(paste0(manuscript_root, "cost_effectiveness_", index ,".png"), plot = cost_effectiveness, dpi = dpi, width = fig_width, height = height, units = "cm", limitsize = FALSE)


}

global_result$simulation <- as.factor(global_result$simulation)
mapping <- c("1" = "Utah", "2" = "Leuven", "3" = "Basel")

# Apply the mapping
global_result <- global_result %>%
  mutate(simulation = recode(as.character(simulation), !!!mapping))
global_result$simulation <- as.factor(global_result$simulation)

global_result$strategy <- factor(global_result$strategy, levels = c("SF", "EF", "PF"))

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
    xlim = c(min(global_result$median_delay, na.rm = TRUE) - 55, max(global_result$median_delay, na.rm = TRUE) - 35)
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

my_comparisons <- list( c("EF", "SF"), c("SF", "PF") )

global_result_norm <- global_result %>%
  group_by(iteration, simulation) %>%
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
  scale_y_continuous(limits = c(min(global_result_norm$nttd), max(global_result_norm$nttd) + length(my_comparisons) * 0.08)) +

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