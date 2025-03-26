###############################################################################
# 1) Libraries & Reading Data
###############################################################################
library(dplyr)
library(stringr)
library(tidyr)
library(purrr)
library(ggplot2)
#install.packages("ineq")  # if not already installed
library(ineq)

# Source your custom SLiM parser, if needed (remove if not actually needed here)
source("/mnt/data2/avalecha/slim/SLiM/My_scripts/Andro_4/functions_Andro_4_lab_pc.R")

# Directory where .rds files are stored
results_dir <- "/mnt/data2/avalecha/slim/SLiM/Simulation_Results_Andro_4/beta_2_20_reps_2/"

# List all RDS files that match sim_s_..._r_..._rep_...
result_files <- list.files(path = results_dir, pattern = "^sim_s_.*_r_.*_rep_.*\\.rds$", full.names = TRUE)

# Read & parse replicate from each filename
simulation_list <- lapply(result_files, function(file_path) {
  df <- readRDS(file_path)
  
  # e.g. "sim_s_0.1_r_0.5_rep_3.rds" -> "3"
  rep_val <- str_extract(basename(file_path), "(?<=_rep_)\\d+")
  df$replicate <- as.integer(rep_val)
  
  df
})

# Combine into one data frame
combined_results <- bind_rows(simulation_list)

# (Optional) Save combined results
out_rds <- file.path(results_dir, "combined_results.rds")
if (!file.exists(out_rds)) {
  saveRDS(combined_results, file = out_rds)
} else {
  message("File 'combined_results.rds' already exists. Skipping save.")
}

# Create a base directory for plots
plots_dir <- "/mnt/data2/avalecha/slim/SLiM/Simulation_Results_Andro_4/beta_2_20_reps_2/Plots_second_set"
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}


##########################################################
# Plot 3 (3_1, 3_2,.... one for each s,r combo) :  Scatter Plot of male and herma sel coeffs
##########################################################


# For example:
sel_coeff_plots_dir <- file.path(plots_dir, "sel_coeff_plots")

# Create the sub‐folder if it doesn't exist
if (!dir.exists(sel_coeff_plots_dir)) {
  dir.create(sel_coeff_plots_dir, recursive = TRUE)
}



gens_to_keep <- c(1, seq(20, 200, by = 20))

# Build a single long data frame containing polymorphic + fixed
all_sel_df <- combined_results %>%
  filter(generation %in% gens_to_keep) %>%
  rowwise() %>%
  do({
    # We'll extract three tibbles (polymorphic, freq=1, freq=0) and combine them
    tib_poly <- tibble(
      locus_position = .$position_poly_loci,
      male_sel       = .$poly_pos_male_sels,
      herma_sel      = .$poly_pos_herma_sels,
      freq_locus     = .$freq_poly_loci,
      locus_type     = rep("polymorphic", length(.$position_poly_loci))
    )
    
    tib_fixed1 <- tibble(
      locus_position = .$freq_1_fixed_loci,
      male_sel       = .$freq_1_male_sels,
      herma_sel      = .$freq_1_herma_sels,
      freq_locus     = rep(1.0, length(.$freq_1_fixed_loci)),
      locus_type     = rep("fixed_1", length(.$freq_1_fixed_loci))
    )
    
    tib_fixed0 <- tibble(
      locus_position = .$freq_0_fixed_loci,
      male_sel       = .$freq_0_male_sels,
      herma_sel      = .$freq_0_herma_sels,
      freq_locus     = rep(0.0, length(.$freq_0_fixed_loci)),
      locus_type     = rep("fixed_0", length(.$freq_0_fixed_loci))
    )
    
    # Combine them into one data frame for this generation
    df_all <- bind_rows(tib_poly, tib_fixed1, tib_fixed0)
    
    # Add columns for s, r, generation
    df_all$s          <- .$s
    df_all$r          <- .$r
    df_all$generation <- .$generation
    df_all$replicate <- .$replicate
    
    df_all
  }) %>%
  ungroup()


# Assume you have 'all_sel_df' with columns:
#  s, r, generation, locus_position, male_sel, herma_sel, freq_locus, locus_type

# Get the unique (s, r) pairs
param_pairs <- all_sel_df %>%
  distinct(s, r,replicate) %>%
  arrange(s, r,replicate)


### Comment for the moment #####

# # For each (s, r) combination, filter and plot
# for (idx in seq_len(nrow(param_pairs))) {
#   
#   curr_s   <- param_pairs$s[idx]
#   curr_r   <- param_pairs$r[idx]
#   curr_rep <- param_pairs$replicate[idx]
#   
#   # Filter data for this one replicate
#   df_sub <- all_sel_df %>%
#     filter(s == curr_s, r == curr_r, replicate == curr_rep)
#   
#   # Determine symmetrical axis limits
#   max_val <- max(abs(df_sub$male_sel), abs(df_sub$herma_sel), na.rm = TRUE)
#   
#   # Build the plot
#   scatter_sel_plot <- ggplot(df_sub, aes(
#     x     = male_sel,
#     y     = herma_sel,
#     color = freq_locus,
#     shape = locus_type
#   )) +
#     geom_point(alpha = 0.7) +
#     scale_color_viridis_c(option = "inferno") +
#     facet_wrap(~ generation, labeller = label_both) +
#     geom_abline(slope = -1, intercept = 0, color="gray50", linetype="dashed") +
#     coord_equal(xlim = c(-max_val, max_val), ylim = c(-max_val, max_val)) +
#     labs(
#       x = "Male Selection Coefficient",
#       y = "Hermaphrodite Selection Coefficient",
#       color = "Frequency",
#       shape = "Locus Type",
#       title = paste(
#         "Scatter of (Male, Herma) Sel Coeffs, s=",
#         curr_s, "r=", curr_r, "rep=", curr_rep
#       )
#     ) +
#     theme_bw(base_size = 14)
#   
#   ggsave(
#     filename = paste0("ScatterSel_s", curr_s, "_r", curr_r,
#                       "_rep_", curr_rep, ".png"),
#     plot = scatter_sel_plot,
#     path = sel_coeff_plots_dir,
#     width = 10,
#     height = 6,
#     dpi = 300
#   )
# }



# ##########################################################
# # Plot 4 Sliding Window analysis - Chromosome Organisation
# ##########################################################
#

########################################## with the sliding window analysis

#
# If you want the total set of possible positions (0..L-1),
# you can define it from your maximum L or from union of all locus_position.
all_positions <- sort(unique(all_sel_df$locus_position))



# We'll build a new data frame, "sw_results", that has columns:
#  s, r, generation, window_center, avg_male, avg_herma
sw_results <- all_sel_df %>%
  filter(generation %in% gens_to_keep) %>%
  group_by(s, r, generation, replicate) %>%
  do({
    # For the current (s, r, generation) subset:
    df_grp <- .
    
    # We'll unify with "all_positions" so that any missing locus is freq=0
    df_unified <- tibble(
      locus_position = all_positions
    ) %>%
      left_join(df_grp, by = "locus_position") %>%
      # Replace NAs appropriately:
      mutate(
        # freq_locus: if not in df_grp => 0
        freq_locus = ifelse(is.na(freq_locus), 0, freq_locus),
        male_sel   = ifelse(is.na(male_sel),   0, male_sel),
        herma_sel  = ifelse(is.na(herma_sel),  0, herma_sel)
      )
    
    # Now we have vectors of positions, freq, male_sel, herma_sel for this generation
    pos_vec   <- df_unified$locus_position
    freq_vec  <- df_unified$freq_locus
    male_vec  <- df_unified$male_sel
    herma_vec <- df_unified$herma_sel
    
    # Do the sliding window
    # (Adjust window_size, step to taste)
    sw_df <- sliding_window_sel_freq(
      positions = pos_vec,
      freq      = freq_vec,
      male_sel  = male_vec,
      herma_sel = herma_vec,
      window_size = 20,
      step       = 1
    )
    
    sw_df
  }) %>%
  ungroup()



################################################################################
## 1) Sliding window for "male-benefit" loci
################################################################################
sw_results_male_benefit <- all_sel_df %>%
  filter(generation %in% gens_to_keep) %>%
  group_by(s, r, generation, replicate) %>%
  do({
    df_grp <- .
    
    # For each generation subset, unify with all_positions
    df_unified <- tibble(locus_position = all_positions) %>%
      left_join(df_grp, by = "locus_position") %>%
      mutate(
        # freq_locus_male_benefit = 1 if male_sel>0, else 0
        freq_locus_male_benefit = ifelse(male_sel > 0, 1, 0)
      )
    
    pos_vec  <- df_unified$locus_position
    freq_vec <- df_unified$freq_locus_male_benefit
    male_vec  <- df_unified$male_sel
    herma_vec <- df_unified$herma_sel
    
    
    # In the sliding window function, the 'male_sel' and 'herma_sel' arguments 
    # are not actually used for the frequency-based result, so you can pass 
    # dummy zeros if you want. Or, if your 'sliding_window_sel_freq()' 
    # specifically needs them, you can keep them as is but they won't matter 
    # for the final freq-based output.
    sw_df <- sliding_window_sel_freq(
      positions  = pos_vec,
      freq       = freq_vec,
      male_sel   = male_vec,
      herma_sel  = herma_vec,
      window_size = 20,
      step       = 1
    )
    
    sw_df
  }) %>%
  ungroup()


################################################################################
## 2) Sliding window for "herma-benefit" loci
################################################################################
## This is just the complement of the above (where male_sel <= 0 => herma-benefit).
## We'll define freq_locus_herma_benefit = 1 - freq_locus_male_benefit.
sw_results_herma_benefit <- all_sel_df %>%
  filter(generation %in% gens_to_keep) %>%
  group_by(s, r, generation, replicate) %>%
  do({
    df_grp <- .
    
    df_unified <- tibble(locus_position = all_positions) %>%
      left_join(df_grp, by = "locus_position") %>%
      mutate(
        freq_locus_male_benefit  = ifelse(male_sel > 0, 1, 0),
        freq_locus_herma_benefit = 1 - freq_locus_male_benefit
      )
    
    pos_vec  <- df_unified$locus_position
    freq_vec <- df_unified$freq_locus_herma_benefit
    male_vec  <- df_unified$male_sel
    herma_vec <- df_unified$herma_sel
    
    sw_df <- sliding_window_sel_freq(
      positions  = pos_vec,
      freq       = freq_vec,
      male_sel   = male_vec, 
      herma_sel  = herma_vec, 
      window_size = 20,
      step       = 1
    )
    
    sw_df
  }) %>%
  ungroup()




# A subdirectory for saving plots
sw_plots_dir <- file.path(plots_dir, "sliding_window_plots_new")
if (!dir.exists(sw_plots_dir)) {
  dir.create(sw_plots_dir, recursive = TRUE)
}

sw_paramrep <- sw_results %>%
  distinct(s, r, replicate) %>%
  arrange(s, r, replicate)


### Comment for the moment #####

# 
# for (i in seq_len(nrow(sw_paramrep))) {
#   
#   curr_s   <- sw_paramrep$s[i]
#   curr_r   <- sw_paramrep$r[i]
#   curr_rep <- sw_paramrep$replicate[i]
#   
#   # Filter for this single replicate
#   df_sub <- sw_results %>%
#     filter(s == curr_s, r == curr_r, replicate == curr_rep)
#   
#   # Basic line plot of (avg_male, avg_herma) vs window_center, faceted by generation
#   p <- ggplot(df_sub, aes(x = window_center)) +
#     geom_line(aes(y = avg_male, color = "Male"), size = 1) +
#     geom_line(aes(y = avg_herma, color = "Hermaphrodite"), size = 1) +
#     facet_wrap(~ generation, labeller = label_both) +
#     scale_color_manual(values = c("Male" = "blue", "Hermaphrodite" = "red")) +
#     labs(
#       x = "Chromosome Position (Window Center)",
#       y = "Avg(Selection Coeff × Frequency)",
#       color = "Type",
#       title = paste0(
#         "Sliding-Window Weighted Selection, s=", curr_s,
#         ", r=", curr_r, ", rep=", curr_rep
#       )
#     ) +
#     theme_bw(base_size = 14)
#   
#   # Optionally add index text or anything else
#   # ...
#   
#   ggsave(
#     filename = paste0("SW_SelFreq_s", curr_s, "_r", curr_r,
#                       "_rep_", curr_rep, ".png"),
#     plot = p,
#     path = sw_plots_dir,
#     width = 10,
#     height = 6,
#     dpi = 300
#   )
# }
# 
#  
# 


for (i in seq_len(nrow(sw_paramrep))) {
  
  curr_s   <- sw_paramrep$s[i]
  curr_r   <- sw_paramrep$r[i]
  curr_rep <- sw_paramrep$replicate[i]
  
  # Filter for this single replicate
  df_sub <- sw_results %>%
    filter(s == curr_s, r == curr_r, replicate == curr_rep)
  #df_sub_male_benefit <- rep1_male_benefit %>%
  #  filter(s == curr_s, r == curr_r, replicate == curr_rep)
  df_sub_herma_benefit <- sw_results_herma_benefit %>%
    filter(s == curr_s, r == curr_r, replicate == curr_rep)
  
  p <- ggplot() +
    # Raw measurement curves
    geom_line(data = df_sub, aes(x = window_center, y = avg_male, color = "Male", linetype = "Actual"), size = 1) +
    geom_line(data = df_sub, aes(x = window_center, y = avg_herma, color = "Hermaphrodite", linetype = "Actual"), size = 1) +
    # Male benefit curves
    # geom_line(data = df_sub_male_benefit, aes(x = window_center, y = avg_male, color = "Male", linetype = "Male benefit"), size = 1) +
    #  geom_line(data = df_sub_male_benefit, aes(x = window_center, y = avg_herma, color = "Hermaphrodite", linetype = "Male benefit"), size = 1) +
    # Hermaphrodite benefit curves
    geom_line(data = df_sub_herma_benefit, aes(x = window_center, y = avg_male, color = "Male", linetype = "Hermaphrodite benefit"), size = 1) +
    geom_line(data = df_sub_herma_benefit, aes(x = window_center, y = avg_herma, color = "Hermaphrodite", linetype = "Hermaphrodite benefit"), size = 1) +
    facet_wrap(~ generation, labeller = label_both) +
    scale_color_manual(values = c("Male" = "blue", "Hermaphrodite" = "red")) +
    scale_linetype_manual(values = c("Actual" = 1, "Hermaphrodite benefit" = 3)) +
    labs(
      x = "Chromosome Position (Window Center)",
      y = "Avg(Selection Coeff × Frequency)",
      color = "Sex",
      linetype = "Treatment",
      title = paste0(
        "Sliding-Window Weighted Selection, s=", curr_s,
        ", r=", curr_r, ", rep=", curr_rep
      )
    ) +
    theme_bw(base_size = 14)
  
  
  # Optionally add index text or anything else
  # ...

  ggsave(
    filename = paste0("SW_SelFreq_s", curr_s, "_r", curr_r,
                      "_rep_", curr_rep, ".png"),
    plot = p,
    path = sw_plots_dir,
    width = 10,
    height = 6,
    dpi = 300
  )
}




# 
# 
# #############################################################
# #Gini coefficients with time (Absolute difference not centered curves)
# #############################################################
# 
# Compute the Gini coefficient for each (s, r, replicate, generation) using sw_results
gini_time_df <- sw_results %>%
  mutate(abs_diff = abs(avg_male - avg_herma)) %>%
  group_by(s, r, replicate, generation) %>%
  summarize(gini_val = Gini(abs_diff, na.rm = TRUE), .groups = "drop")


 # Create a directory for the Gini vs Time plots if it doesn't exist
 gini_plots_dir <- file.path(plots_dir, "gini_init_distribution_abs_diff_not_centered")
if (!dir.exists(gini_plots_dir)) {
  dir.create(gini_plots_dir, recursive = TRUE)
}

gen1_gini <- gini_time_df %>% filter(generation == 1)

p_gini_dist <- ggplot(gen1_gini, aes(x = gini_val)) +
  geom_histogram(bins = 20, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(x = "Gini Value", y = "Count",
       title = "Distribution of Gini Values at Generation 1")

# Save the plot
   ggsave(
     filename = paste0("Gini_distribution_Gen_1.png"),
     plot = p_gini_dist,
     path = gini_plots_dir,
     width = 8,
     height = 5,
     dpi = 300
   )

# # Summarize average Gini per s, r, generation across replicates
# gini_summary <- gini_time_df %>%
#   group_by(s, r, generation) %>%
#   summarize(
#     mean_gini = mean(gini_val, na.rm = TRUE),
#     se_gini   = sd(gini_val, na.rm = TRUE) / sqrt(n()),
#     .groups = "drop"
#   )
# 
# 
# # Create a directory for the Gini vs Time plots if it doesn't exist
# gini_plots_dir <- file.path(plots_dir, "gini_vs_time_abs_diff_not_centered")
# if (!dir.exists(gini_plots_dir)) {
#   dir.create(gini_plots_dir, recursive = TRUE)
# }
# 
# # For each (s, r) combination, plot per replicate curves with an average trend line
# param_pairs_gini <- gini_time_df %>% distinct(s, r)
# 
# for (i in seq_len(nrow(param_pairs_gini))) {
#   cs <- param_pairs_gini$s[i]
#   cr <- param_pairs_gini$r[i]
#   
#   # Filter data for current s and r
#   df_sr_gini <- gini_time_df %>%
#     filter(s == cs, r == cr)
#   
#   # Filter the average data
#   avg_line <- gini_summary %>%
#     filter(s == cs, r == cr)
#   
#   p_gini <- ggplot() +
#     # Plot the replicate-specific curves
#     geom_line(data = df_sr_gini,
#               aes(x = generation, y = gini_val, color = factor(replicate)),
#               size = 1) +
#     geom_point(data = df_sr_gini,
#                aes(x = generation, y = gini_val, color = factor(replicate)),
#                size = 2) +
#     # Overlay the average line as a bold black line
#     geom_line(data = avg_line,
#               aes(x = generation, y = mean_gini),
#               size = 1.5, color = "black") +
#     # (Optional) Overlay error bars for the average, if desired
#     geom_errorbar(data = avg_line,
#                   aes(x = generation, ymin = mean_gini - se_gini, ymax = mean_gini + se_gini),
#                   width = 2, color = "black") +
#     labs(
#       x = "Generation",
#       y = "Gini Coefficient",
#       color = "Replicate",
#       title = paste("Gini Coefficient vs Generation (s =", cs, ", r =", cr, ")")
#     ) +
#     theme_bw(base_size = 14) + 
#     coord_cartesian(ylim = c(0, 1))
#   
#   
#   # Save the plot
#   ggsave(
#     filename = paste0("Gini_vs_Time_s", cs, "_r", cr, ".png"),
#     plot = p_gini,
#     path = gini_plots_dir,
#     width = 8,
#     height = 5,
#     dpi = 300
#   )
# }
# 
# 
# # ###############################################################################
# # Produce Comparision (different s,r parameters) Plots of gini coefficients
# # ###############################################################################
# 
# # Set a directory to save these plots
# comparison_plots_dir <- file.path(plots_dir, "gini_comparison_abs_diff_not_centered")
# if (!dir.exists(comparison_plots_dir)) {
#   dir.create(comparison_plots_dir, recursive = TRUE)
# }
# 
# #############################################
# # 1. For each constant s, vary r
# #############################################
# unique_s <- sort(unique(gini_summary$s))
# for (s_val in unique_s) {
#   df_s <- gini_summary %>% filter(s == s_val)
#   
#   p_s <- ggplot(df_s, aes(x = generation, y = mean_gini, color = factor(r))) +
#     geom_line(size = 1) +
#     geom_point(size = 2) +
#     geom_errorbar(aes(ymin = mean_gini - se_gini, ymax = mean_gini + se_gini), 
#                   width = 2, color = "gray40", alpha = 0.5) +
#     labs(
#       title = paste("Mean Gini vs Generation (s =", s_val, ")"),
#       x = "Generation",
#       y = "Mean Gini Coefficient",
#       color = "r"
#     ) +
#     theme_bw(base_size = 14) + 
#     coord_cartesian(ylim = c(0, 1))
#   
#   
#   ggsave(
#     filename = paste0("Mean_Gini_vs_Generation_s", s_val, ".png"),
#     plot = p_s,
#     path = comparison_plots_dir,
#     width = 8,
#     height = 5,
#     dpi = 300
#   )
# }
# 
# #############################################
# # 2. For each constant r, vary s
# #############################################
# unique_r <- sort(unique(gini_summary$r))
# for (r_val in unique_r) {
#   df_r <- gini_summary %>% filter(r == r_val)
#   
#   p_r <- ggplot(df_r, aes(x = generation, y = mean_gini, color = factor(s))) +
#     geom_line(size = 1) +
#     geom_point(size = 2) +
#     geom_errorbar(aes(ymin = mean_gini - se_gini, ymax = mean_gini + se_gini), 
#                   width = 2, color = "gray40", alpha = 0.5) +
#     labs(
#       title = paste("Mean Gini vs Generation (r =", r_val, ")"),
#       x = "Generation",
#       y = "Mean Gini Coefficient",
#       color = "s"
#     ) +
#     theme_bw(base_size = 14) +
#     coord_cartesian(ylim = c(0, 1))
#   
#   
#   ggsave(
#     filename = paste0("Mean_Gini_vs_Generation_r", r_val, ".png"),
#     plot = p_r,
#     path = comparison_plots_dir,
#     width = 8,
#     height = 5,
#     dpi = 300
#   )
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #############################################################
# #Gini coefficients with time (Absolute difference  centered curves)
# #############################################################



gini_time_df <- sw_results %>%
  # 1) Group by the keys so we can compute means within each replicate & generation
  group_by(s, r, replicate, generation) %>%
  # 2) Subtract the mean from avg_male & avg_herma within each group
  mutate(
    male_centered  = avg_male  - mean(avg_male,  na.rm = TRUE),
    herma_centered = avg_herma - mean(avg_herma, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  # 3) Now define the abs_diff from these centered values
  mutate(abs_diff = abs(male_centered - herma_centered)) %>%
  # 4) Finally compute Gini per group, as usual
  group_by(s, r, replicate, generation) %>%
  summarize(
    gini_val = Gini(abs_diff, na.rm = TRUE),
    .groups  = "drop"
  )

   
# Create a directory for the Gini vs Time plots if it doesn't exist
gini_plots_dir <- file.path(plots_dir, "gini_init_distribution_abs_diff_centered")
if (!dir.exists(gini_plots_dir)) {
     dir.create(gini_plots_dir, recursive = TRUE)
 }
   
gen1_gini <- gini_time_df %>% filter(generation == 1)
   
p_gini_dist <- ggplot(gen1_gini, aes(x = gini_val)) +
     geom_histogram(bins = 20, fill = "skyblue", color = "black", alpha = 0.7) +
     labs(x = "Gini Value", y = "Count",
          title = "Distribution of Gini Values at Generation 1")
   
# Save the plot
ggsave(
     filename = paste0("Gini_distribution_Gen_1.png"),
     plot = p_gini_dist,
     path = gini_plots_dir,
     width = 8,
     height = 5,
     dpi = 300
   )
   
   
# 
# 
# # Summarize average Gini per s, r, generation across replicates
# gini_summary <- gini_time_df %>%
#   group_by(s, r, generation) %>%
#   summarize(
#     mean_gini = mean(gini_val, na.rm = TRUE),
#     se_gini   = sd(gini_val, na.rm = TRUE) / sqrt(n()),
#     .groups = "drop"
#   )
# 
# 
# # Create a directory for the Gini vs Time plots if it doesn't exist
# gini_plots_dir <- file.path(plots_dir, "gini_vs_time_abs_diff_centered")
# if (!dir.exists(gini_plots_dir)) {
#   dir.create(gini_plots_dir, recursive = TRUE)
# }
# 
# # For each (s, r) combination, plot per replicate curves with an average trend line
# param_pairs_gini <- gini_time_df %>% distinct(s, r)
# 
# for (i in seq_len(nrow(param_pairs_gini))) {
#   cs <- param_pairs_gini$s[i]
#   cr <- param_pairs_gini$r[i]
#   
#   # Filter data for current s and r
#   df_sr_gini <- gini_time_df %>%
#     filter(s == cs, r == cr)
#   
#   # Filter the average data
#   avg_line <- gini_summary %>%
#     filter(s == cs, r == cr)
#   
#   p_gini <- ggplot() +
#     # Plot the replicate-specific curves
#     geom_line(data = df_sr_gini,
#               aes(x = generation, y = gini_val, color = factor(replicate)),
#               size = 1) +
#     geom_point(data = df_sr_gini,
#                aes(x = generation, y = gini_val, color = factor(replicate)),
#                size = 2) +
#     # Overlay the average line as a bold black line
#     geom_line(data = avg_line,
#               aes(x = generation, y = mean_gini),
#               size = 1.5, color = "black") +
#     # (Optional) Overlay error bars for the average, if desired
#     geom_errorbar(data = avg_line,
#                   aes(x = generation, ymin = mean_gini - se_gini, ymax = mean_gini + se_gini),
#                   width = 2, color = "black") +
#     labs(
#       x = "Generation",
#       y = "Gini Coefficient",
#       color = "Replicate",
#       title = paste("Gini Coefficient vs Generation (s =", cs, ", r =", cr, ")")
#     ) +
#     theme_bw(base_size = 14) + 
#     coord_cartesian(ylim = c(0, 1))
#   
#   
#   # Save the plot
#   ggsave(
#     filename = paste0("Gini_vs_Time_s", cs, "_r", cr, ".png"),
#     plot = p_gini,
#     path = gini_plots_dir,
#     width = 8,
#     height = 5,
#     dpi = 300
#   )
# }
# 
# # ###############################################################################
# # Produce Comparision (different s,r parameters) Plots of gini coefficients
# # ###############################################################################
# 
# # Set a directory to save these plots
# comparison_plots_dir <- file.path(plots_dir, "gini_comparison_abs_diff_centered")
# if (!dir.exists(comparison_plots_dir)) {
#   dir.create(comparison_plots_dir, recursive = TRUE)
# }
# 
# #############################################
# # 1. For each constant s, vary r
# #############################################
# unique_s <- sort(unique(gini_summary$s))
# for (s_val in unique_s) {
#   df_s <- gini_summary %>% filter(s == s_val)
#   
#   p_s <- ggplot(df_s, aes(x = generation, y = mean_gini, color = factor(r))) +
#     geom_line(size = 1) +
#     geom_point(size = 2) +
#     geom_errorbar(aes(ymin = mean_gini - se_gini, ymax = mean_gini + se_gini), 
#                   width = 2, color = "gray40", alpha = 0.5) +
#     labs(
#       title = paste("Mean Gini vs Generation (s =", s_val, ")"),
#       x = "Generation",
#       y = "Mean Gini Coefficient",
#       color = "r"
#     ) +
#     theme_bw(base_size = 14) + 
#     coord_cartesian(ylim = c(0, 1))
#   
#   
#   ggsave(
#     filename = paste0("Mean_Gini_vs_Generation_s", s_val, ".png"),
#     plot = p_s,
#     path = comparison_plots_dir,
#     width = 8,
#     height = 5,
#     dpi = 300
#   )
# }
# 
# #############################################
# # 2. For each constant r, vary s
# #############################################
# unique_r <- sort(unique(gini_summary$r))
# for (r_val in unique_r) {
#   df_r <- gini_summary %>% filter(r == r_val)
#   
#   p_r <- ggplot(df_r, aes(x = generation, y = mean_gini, color = factor(s))) +
#     geom_line(size = 1) +
#     geom_point(size = 2) +
#     geom_errorbar(aes(ymin = mean_gini - se_gini, ymax = mean_gini + se_gini), 
#                   width = 2, color = "gray40", alpha = 0.5) +
#     labs(
#       title = paste("Mean Gini vs Generation (r =", r_val, ")"),
#       x = "Generation",
#       y = "Mean Gini Coefficient",
#       color = "s"
#     ) +
#     theme_bw(base_size = 14) +
#     coord_cartesian(ylim = c(0, 1))
#   
#   
#   ggsave(
#     filename = paste0("Mean_Gini_vs_Generation_r", r_val, ".png"),
#     plot = p_r,
#     path = comparison_plots_dir,
#     width = 8,
#     height = 5,
#     dpi = 300
#   )
# }










#############################################################
#Gini coefficients with time (Absolute difference  w.r.t mean for each curve)
#############################################################

gini_time_df <- sw_results %>%
  group_by(s, r, replicate, generation) %>%
  mutate(
    # 1) Center each: subtract the mean from avg_male or avg_herma
    #    then take absolute value if desired. Typically, "centered" means
    #    (value - mean(value)), but you can adapt as needed:
    male_centered  = abs(avg_male  - mean(avg_male,  na.rm = TRUE)),
    herma_centered = abs(avg_herma - mean(avg_herma, na.rm = TRUE))
  ) %>%
  ungroup() %>%
  # 2) Summarize Gini for male_centered and herma_centered separately
  group_by(s, r, replicate, generation) %>%
  summarize(
    gini_male  = Gini(male_centered,  na.rm = TRUE),
    gini_herma = Gini(herma_centered, na.rm = TRUE),
    .groups    = "drop"
  )


# Create a directory for the Gini vs Time plots if it doesn't exist
gini_plots_dir <- file.path(plots_dir, "gini_init_distribution_abs_separate_curves")
if (!dir.exists(gini_plots_dir)) {
  dir.create(gini_plots_dir, recursive = TRUE)
}

gen1_gini <- gini_time_df %>% filter(generation == 1)

gen1_gini_long <- gen1_gini %>%
  pivot_longer(
    cols = c(gini_male, gini_herma),
    names_to = "sex",
    values_to = "gini_value"
  )

p_gini_dist <- ggplot(gen1_gini_long, aes(x = gini_value)) +
  geom_histogram(bins = 20, fill = "skyblue", color = "black", alpha = 0.7) +
  facet_wrap(~ sex, scales = "free_y") +
  labs(
    x = "Gini Value",
    y = "Count",
    title = "Distribution of Gini Values at Generation 1"
  )


# Save the plot
ggsave(
  filename = paste0("Gini_distribution_Gen_1.png"),
  plot = p_gini_dist,
  path = gini_plots_dir,
  width = 8,
  height = 5,
  dpi = 300
)

gini_summary <- gini_time_df %>%
  group_by(s, r, generation) %>%
  summarize(
    mean_gini_male  = mean(gini_male, na.rm = TRUE),
    se_gini_male    = sd(gini_male,   na.rm = TRUE)/sqrt(n()),
    
    mean_gini_herma = mean(gini_herma, na.rm = TRUE),
    se_gini_herma   = sd(gini_herma,   na.rm = TRUE)/sqrt(n()),
    .groups         = "drop"
  )

# Directory for the new Gini plots
gini_plots_dir <- file.path(plots_dir, "gini_vs_time_abs_separate_curves_1")
if(!dir.exists(gini_plots_dir)) {
  dir.create(gini_plots_dir, recursive = TRUE)
}

# Distinct (s, r)
param_pairs_gini <- gini_time_df %>%
  distinct(s, r)

for(i in seq_len(nrow(param_pairs_gini))) {
  cs <- param_pairs_gini$s[i]
  cr <- param_pairs_gini$r[i]
  
  # 1) Filter replicate-level data for just this (s, r)
  df_sr_gini <- gini_time_df %>%
    filter(s == cs, r == cr)
  
  # New: Pivot replicate-level data into long format
  df_sr_long <- df_sr_gini %>%
    pivot_longer(
      cols = c(gini_male, gini_herma),
      names_to = "sex",
      values_to = "gini",
      names_prefix = "gini_"
    )
  
  
  # 2) Filter the average lines for that (s, r)
  avg_line <- gini_summary %>%
    filter(s == cs, r == cr)
  
  # New: Pivot average summary data into long format for each sex
  avg_line_long_male <- avg_line %>%
    mutate(sex = "male",
           mean_gini = mean_gini_male,
           se = se_gini_male)
  avg_line_long_herma <- avg_line %>%
    mutate(sex = "herma",
           mean_gini = mean_gini_herma,
           se = se_gini_herma)
  avg_line_long <- bind_rows(avg_line_long_male, avg_line_long_herma)
  
  
  # Build ggplot using the long-format data and facet_wrap to create separate subplots
  p_gini <- ggplot() +
    
    # Plot replicate-level lines and points for both sexes
    geom_line(
      data = df_sr_long,
      aes(x = generation, y = gini, color = factor(replicate)),
      size = 1
    ) +
    geom_point(
      data = df_sr_long,
      aes(x = generation, y = gini, color = factor(replicate)),
      size = 2
    ) +
    
    # Overlay average lines and error bars for each sex
    geom_line(
      data = avg_line_long,
      aes(x = generation, y = mean_gini),
      size = 1.5, color = "black"
    ) +
    geom_errorbar(
      data = avg_line_long,
      aes(x = generation,
          ymin = mean_gini - se,
          ymax = mean_gini + se),
      width = 2, color = "black"
    ) +
    
    # Facet by sex to create separate subplots
    facet_wrap(~ sex, ncol = 1) +
    
    labs(
      x = "Generation",
      y = "Gini Coefficient (Centered)",
      color = "Replicate",
      title = paste0("Gini Coefficient by Sex, s=", cs, ", r=", cr)
    ) +
    theme_bw(base_size = 14) +
    coord_cartesian(ylim = c(0, 1))
  
  # Save
  outname <- paste0("DoubleGini_s", cs, "_r", cr, ".png")
  ggsave(
    filename = outname,
    plot = p_gini,
    path = gini_plots_dir,
    width = 8,
    height = 5,
    dpi = 300
  )
}

# ###############################################################################
# Produce Comparision (different s,r parameters) Plots of gini coefficients
# ###############################################################################

# Set a directory to save these plots
comparison_plots_dir <- file.path(plots_dir, "gini_comparison_abs_separate_curves_1")
if (!dir.exists(comparison_plots_dir)) {
  dir.create(comparison_plots_dir, recursive = TRUE)
}

#############################################
# 1. For each constant s, vary r (with two curves: male & herma)
#############################################

unique_s <- sort(unique(gini_summary$s))
for (s_val in unique_s) {
  df_s <- gini_summary %>% filter(s == s_val)
  
  # Pivot to long format for separate male and herma Gini values and errors
  df_long <- df_s %>%
    pivot_longer(cols = c(mean_gini_male, mean_gini_herma),
                 names_to = "sex", values_to = "mean_gini") %>%
    mutate(sex = if_else(sex == "mean_gini_male", "male", "herma"),
           se = if_else(sex == "male", se_gini_male, se_gini_herma))
  
  p_s <- ggplot(df_long, aes(x = generation, y = mean_gini, color = factor(r))) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = mean_gini - se, ymax = mean_gini + se),
                  width = 2, color = "gray40", alpha = 0.5) +
    facet_wrap(~ sex, ncol = 1) +  # NEW: split by sex
    labs(
      title = paste("Mean Gini vs Generation (s =", s_val, ")"),
      x = "Generation",
      y = "Mean Gini Coefficient",
      color = "r"
    ) +
    theme_bw(base_size = 14) +
    coord_cartesian(ylim = c(0, 1))
  
  ggsave(
    filename = paste0("Double_Mean_Gini_vs_Generation_s", s_val, ".png"),
    plot = p_s,
    path = comparison_plots_dir,
    width = 8,
    height = 5,
    dpi = 300
  )
}

#############################################
# 2. For each constant r, vary s (with two curves: male & herma)
#############################################
unique_r <- sort(unique(gini_summary$r))
for (r_val in unique_r) {
  df_r <- gini_summary %>% filter(r == r_val)
  
  # Pivot to long format for separate male and herma Gini values and errors
  df_long <- df_r %>%
    pivot_longer(cols = c(mean_gini_male, mean_gini_herma),
                 names_to = "sex", values_to = "mean_gini") %>%
    mutate(sex = if_else(sex == "mean_gini_male", "male", "herma"),
           se = if_else(sex == "male", se_gini_male, se_gini_herma))
  
  p_r <- ggplot(df_long, aes(x = generation, y = mean_gini, color = factor(s))) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = mean_gini - se, ymax = mean_gini + se),
                  width = 2, color = "gray40", alpha = 0.5) +
    facet_wrap(~ sex, ncol = 1) +  # NEW: split by sex
    labs(
      title = paste("Mean Gini vs Generation (r =", r_val, ")"),
      x = "Generation",
      y = "Mean Gini Coefficient",
      color = "s"
    ) +
    theme_bw(base_size = 14) +
    coord_cartesian(ylim = c(0, 1))
  
  ggsave(
    filename = paste0("Double_Mean_Gini_vs_Generation_r", r_val, ".png"),
    plot = p_r,
    path = comparison_plots_dir,
    width = 8,
    height = 5,
    dpi = 300
  )
}







#######################################################################
######Summation absolute area of centred curve for each curve (not gini)
#######################################################################


abs_area_time_df <- sw_results %>%
  group_by(s, r, replicate, generation) %>%
  mutate(
    # 1) Center each: subtract the mean from avg_male or avg_herma
    #    then take absolute value if desired. Typically, "centered" means
    #    (value - mean(value)), but you can adapt as needed:
    male_centered  = abs(avg_male  - mean(avg_male,  na.rm = TRUE)),
    herma_centered = abs(avg_herma - mean(avg_herma, na.rm = TRUE))
  ) %>%
  ungroup() %>%
  # 2) Summarize Gini for male_centered and herma_centered separately
  group_by(s, r, replicate, generation) %>%
  summarize(
    abs_area_male  = sum(male_centered, na.rm = TRUE),
    abs_area_herma = sum(herma_centered, na.rm = TRUE),
    .groups    = "drop"
  )



abs_area_summary <- abs_area_time_df %>%
  group_by(s, r, generation) %>%
  summarize(
    mean_abs_area_male  = mean(abs_area_male, na.rm = TRUE),
    se_abs_area_male    = sd(abs_area_male, na.rm = TRUE)/sqrt(n()),
    
    mean_abs_area_herma  = mean(abs_area_herma, na.rm = TRUE),
    se_abs_area_herma    = sd(abs_area_herma, na.rm = TRUE)/sqrt(n()),
    .groups         = "drop"
  )

# Directory for the new Gini plots
abs_area_plots_dir <- file.path(plots_dir, "abs_area_vs_time_separate_curves")
if(!dir.exists(abs_area_plots_dir)) {
  dir.create(abs_area_plots_dir, recursive = TRUE)
}

# Distinct (s, r)
param_pairs_abs_area <- abs_area_time_df %>%
  distinct(s, r)

for(i in seq_len(nrow(param_pairs_abs_area))) {
  cs <- param_pairs_abs_area$s[i]
  cr <- param_pairs_abs_area$r[i]
  
  # 1) Filter replicate-level data for just this (s, r)
  df_sr_abs_area <- abs_area_time_df %>%
    filter(s == cs, r == cr)
  
  # New: Pivot replicate-level data into long format
  df_sr_long <- df_sr_abs_area %>%
    pivot_longer(
      cols = c(abs_area_male, abs_area_herma),
      names_to = "sex",
      values_to = "abs_area",
      names_prefix = "abs_area_"
    )
  
  
  # 2) Filter the average lines for that (s, r)
  avg_line <- abs_area_summary %>%
    filter(s == cs, r == cr)
  
  # New: Pivot average summary data into long format for each sex
  avg_line_long_male <- avg_line %>%
    mutate(sex = "male",
           mean_abs_area = mean_abs_area_male,
           se = se_abs_area_male)
  avg_line_long_herma <- avg_line %>%
    mutate(sex = "herma",
           mean_abs_area = mean_abs_area_herma,
           se = se_abs_area_herma)
  avg_line_long <- bind_rows(avg_line_long_male, avg_line_long_herma)
  
  # Build ggplot using the long-format data and facet_wrap to create separate subplots
  p_abs_area <- ggplot() +
    
    # Plot replicate-level lines and points for both sexes
    geom_line(
      data = df_sr_long,
      aes(x = generation, y = abs_area, color = factor(replicate)),
      size = 1
    ) +
    geom_point(
      data = df_sr_long,
      aes(x = generation, y = abs_area, color = factor(replicate)),
      size = 2
    ) +
    
    # Overlay average lines and error bars for each sex
    geom_line(
      data = avg_line_long,
      aes(x = generation, y = mean_abs_area),
      size = 1.5, color = "black"
    ) +
    geom_errorbar(
      data = avg_line_long,
      aes(x = generation,
          ymin = mean_abs_area - se,
          ymax = mean_abs_area + se),
      width = 2, color = "black"
    ) +
    
    # Facet by sex to create separate subplots
    facet_wrap(~ sex, ncol = 1) +
    
    labs(
      x = "Generation",
      y = "abs_area (Centered)",
      color = "Replicate",
      title = paste0("abs_area by Sex, s=", cs, ", r=", cr)
    ) +
    theme_bw(base_size = 14) +
    coord_cartesian(ylim = c(0, 2))
  
  # Save
  outname <- paste0("Doubleabs_area_s", cs, "_r", cr, ".png")
  ggsave(
    filename = outname,
    plot = p_abs_area,
    path = abs_area_plots_dir,
    width = 8,
    height = 5,
    dpi = 300
  )
}




###################################################################
################### Fracn herma curve > male curve and net difference 
###################################################################
########################################
# 1) Compute the two measures per replicate & generation
########################################
index_df <- sw_results %>%
  group_by(s, r, replicate, generation) %>%
  summarize(
    # Fraction of window_center rows where avg_herma > avg_male
    frac_herma_greater = mean(avg_herma > avg_male, na.rm = TRUE),
    
    # "Net difference" can be sum(...) or mean(...); here we use sum
    net_diff = sum(avg_herma - avg_male, na.rm = TRUE),
    .groups = "drop"
  )

########################################
# 2) Summarize across replicates to get the mean trend + SE
########################################
index_summary <- index_df %>%
  group_by(s, r, generation) %>%
  summarize(
    mean_frac     = mean(frac_herma_greater, na.rm = TRUE),
    se_frac       = sd(frac_herma_greater,   na.rm = TRUE) / sqrt(n()),
    
    mean_net_diff = mean(net_diff, na.rm = TRUE),
    se_net_diff   = sd(net_diff,   na.rm = TRUE) / sqrt(n()),
    .groups       = "drop"
  )

# Directory for the new Gini plots
frac_plots_dir <- file.path(plots_dir, "frac_vs_time")
if(!dir.exists(frac_plots_dir)) {
  dir.create(frac_plots_dir, recursive = TRUE)
}

param_pairs <- index_df %>%
  distinct(s, r) %>%
  arrange(s, r)




for (i in seq_len(nrow(param_pairs))) {
  cs <- param_pairs$s[i]
  cr <- param_pairs$r[i]
  
  # Subset replicate-level data for this (s, r)
  df_sub <- index_df %>%
    filter(s == cs, r == cr)
  
  # Subset summary data (mean + SE) for this (s, r)
  avg_sub <- index_summary %>%
    filter(s == cs, r == cr)
  
  ########### (A) Plot fraction_herma_greater ############
  p_frac <- ggplot() +
    geom_line(
      data = df_sub,
      aes(x = generation, y = frac_herma_greater, color = factor(replicate)),
      size = 1
    ) +
    geom_point(
      data = df_sub,
      aes(x = generation, y = frac_herma_greater, color = factor(replicate)),
      size = 2
    ) +
    geom_line(
      data = avg_sub,
      aes(x = generation, y = mean_frac),
      size = 1.5, color = "black"
    ) +
    geom_errorbar(
      data = avg_sub,
      aes(x = generation, ymin = mean_frac - se_frac, ymax = mean_frac + se_frac),
      width = 2, color = "black"
    ) +
    labs(
      x = "Generation",
      y = "Fraction(avg_herma > avg_male)",
      color = "Replicate",
      title = paste("Fraction(herma>male)\ns =", cs, ", r =", cr)
    ) +
    theme_bw(base_size = 14) +
    coord_cartesian(ylim = c(0, 1))
  
  ########### (B) Plot net_diff = sum(avg_herma - avg_male) ############
  p_diff <- ggplot() +
    geom_line(
      data = df_sub,
      aes(x = generation, y = net_diff, color = factor(replicate)),
      size = 1
    ) +
    geom_point(
      data = df_sub,
      aes(x = generation, y = net_diff, color = factor(replicate)),
      size = 2
    ) +
    geom_line(
      data = avg_sub,
      aes(x = generation, y = mean_net_diff),
      size = 1.5, color = "black"
    ) +
    geom_errorbar(
      data = avg_sub,
      aes(x = generation, ymin = mean_net_diff - se_net_diff, ymax = mean_net_diff + se_net_diff),
      width = 2, color = "black"
    ) +
    labs(
      x = "Generation",
      y = "Net Diff (herma - male)",
      color = "Replicate",
      title = paste("Net Diff(herma - male)\ns =", cs, ", r =", cr)
    ) +
    theme_bw(base_size = 14)
  
  ########### (C) Combine side-by-side with patchwork ###########
  # The '+' operator from patchwork assembles plots:
  combined_plot <- p_frac + p_diff + plot_layout(ncol = 2)
  
  ggsave(
    filename = paste0("FracAndDiff_s", cs, "_r", cr, ".png"),
    plot     = combined_plot,
    path = frac_plots_dir,
    width    = 12,  # wider, since we have 2 plots horizontally
    height   = 5,
    dpi      = 300
  )
}








# ###############################################################################
# # 3) Lorenz Curves with Gini Annotation for each (s, r, replicate), faceted by Generation
# ###############################################################################
# 
# # For each record in sw_results, compute the absolute difference and then, by group, 
# # compute the Lorenz curve coordinates and Gini value.
# # We'll use purrr::map to store the Lorenz curve coordinates as list-columns.
# library(tidyr)
# 
# lorenz_data <- sw_results %>%
#   mutate(abs_diff = abs(avg_male - avg_herma)) %>%
#   group_by(s, r, replicate, generation) %>%
#   summarize(
#     gini_val = Gini(abs_diff, na.rm = TRUE),
#     lorenz_obj = list(Lc(abs_diff)),
#     .groups = "drop"
#   ) %>%
#   mutate(
#     p = map(lorenz_obj, ~ .x$p),
#     L = map(lorenz_obj, ~ .x$L)
#   ) %>%
#   unnest(cols = c(p, L))
# 
# # Create a directory for Lorenz curve plots if it doesn't exist
# lorenz_plots_dir <- file.path(plots_dir, "lorenz_curves")
# if (!dir.exists(lorenz_plots_dir)) {
#   dir.create(lorenz_plots_dir, recursive = TRUE)
# }
# 
# # Get distinct combinations of s, r, and replicate for looping
# lorenz_paramrep <- lorenz_data %>% distinct(s, r, replicate)
# 
# for (i in seq_len(nrow(lorenz_paramrep))) {
#   cs   <- lorenz_paramrep$s[i]
#   cr   <- lorenz_paramrep$r[i]
#   crep <- lorenz_paramrep$replicate[i]
#   
#   df_lorenz <- lorenz_data %>%
#     filter(s == cs, r == cr, replicate == crep)
#   
#   # Create the Lorenz curves plot faceted by generation.
#   # The 45° red dashed line represents perfect equality.
#   p_lorenz <- ggplot(df_lorenz, aes(x = p, y = L)) +
#     geom_line(color = "blue", size = 1) +
#     geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
#     facet_wrap(~ generation, labeller = label_both) +
#     labs(
#       x = "Cumulative share of windows",
#       y = "Cumulative share of absolute difference",
#       title = paste0("Lorenz Curves with Gini Annotation (s=", cs, ", r=", cr, ", rep=", crep, ")")
#     ) +
#     theme_bw(base_size = 14)
#   
#   # For annotation, extract one Gini value per generation
#   gini_annotations <- df_lorenz %>%
#     group_by(generation) %>%
#     summarize(gini_val = unique(gini_val), .groups = "drop")
#   
#   # Annotate each facet with the corresponding Gini coefficient
#   p_lorenz <- p_lorenz +
#     geom_text(
#       data = gini_annotations,
#       aes(x = 0.8, y = 0.2, label = paste("Gini:", round(gini_val, 2))),
#       inherit.aes = FALSE,
#       size = 3,
#       color = "black"
#     )
#   
#   # Save the Lorenz curves plot for this (s, r, replicate) combination
#   ggsave(
#     filename = paste0("Lorenz_Curves_s", cs, "_r", cr, "_rep", crep, ".png"),
#     plot = p_lorenz,
#     path = lorenz_plots_dir,
#     width = 10,
#     height = 6,
#     dpi = 300
#   )
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ###############################################################################
# # 4) Heterozygosity vs Generation for All Replicates
# ###############################################################################
# # We assume combined_results has columns (s, r, replicate, generation, heterozygosity)
# 
# het_plots_dir <- file.path(plots_dir, "het_vs_time")
# if (!dir.exists(het_plots_dir)) {
#   dir.create(het_plots_dir, recursive = TRUE)
# }
# 
# # We'll produce one plot per (s, r), lines for each replicate
# het_param_pairs <- combined_results %>%
#   distinct(s, r) %>%
#   arrange(s, r)
# 
# for (j in seq_len(nrow(het_param_pairs))) {
#   cs <- het_param_pairs$s[j]
#   cr <- het_param_pairs$r[j]
#   
#   df_het_sr <- combined_results %>%
#     filter(s==cs, r==cr)
#   
#   p_het <- ggplot(df_het_sr, aes(x=generation, y=heterozygosity, color=factor(replicate))) +
#     geom_line(size=1) +
#     labs(
#       x="Generation",
#       y="Heterozygosity",
#       color="Replicate",
#       title=paste("Heterozygosity vs Generation (s=", cs, ", r=", cr, ")")
#     ) +
#     theme_bw(base_size=14)
#   
#   ggsave(
#     filename=paste0("Heterozygosity_vs_Time_s", cs, "_r", cr, ".png"),
#     plot=p_het,
#     path=het_plots_dir,
#     width=8,
#     height=5,
#     dpi=300
#   )
# }
# 
# ###############################################################################
# # 5) Produce Sliding-Window Curves (One Plot per (s, r, replicate), Facet by Gen)
# ###############################################################################
# # We already have `sw_results` containing columns:
# #   (s, r, replicate, generation, window_center, avg_male, avg_herma)
# # from your existing code.
# 
# # Create a subdirectory for these new sliding-window facet plots
# sw_facets_dir <- file.path(plots_dir, "sliding_window_facets")
# if (!dir.exists(sw_facets_dir)) {
#   dir.create(sw_facets_dir, recursive = TRUE)
# }
# 
# # Distinct (s, r, replicate)
# sw_paramrep <- sw_results %>%
#   distinct(s, r, replicate) %>%
#   arrange(s, r, replicate)
# 
# for (i in seq_len(nrow(sw_paramrep))) {
#   curr_s   <- sw_paramrep$s[i]
#   curr_r   <- sw_paramrep$r[i]
#   curr_rep <- sw_paramrep$replicate[i]
#   
#   df_sub <- sw_results %>%
#     filter(s == curr_s, r == curr_r, replicate == curr_rep)
#   
#   # We'll plot 'avg_male' (blue) and 'avg_herma' (red) vs window_center
#   # Facet by generation
#   p_sw <- ggplot(df_sub, aes(x = window_center)) +
#     geom_line(aes(y = avg_male, color = "Male"), size = 1) +
#     geom_line(aes(y = avg_herma, color = "Hermaphrodite"), size = 1) +
#     facet_wrap(~ generation, labeller = label_both) +
#     scale_color_manual(values = c("Male" = "blue", "Hermaphrodite" = "red")) +
#     labs(
#       x = "Window Center (Chromosome Position)",
#       y = "Average (Sel * Freq)",
#       color = "Sex",
#       title = paste0(
#         "Sliding-Window Curves, s=", curr_s,
#         ", r=", curr_r,
#         ", rep=", curr_rep
#       )
#     ) +
#     theme_bw(base_size = 14)
#   
#   ggsave(
#     filename = paste0("SW_Curves_s", curr_s, "_r", curr_r, "_rep", curr_rep, ".png"),
#     plot = p_sw,
#     path = sw_facets_dir,
#     width = 10,
#     height = 6,
#     dpi = 300
#   )
# }
# 
# 


########################################
###Fitness plots
########################################


# Filter for generation 200 and scale fitness values for each row
fitness_df <- combined_results %>%
  filter(generation == 200) %>%
  rowwise() %>%
  mutate(
    scale_factor = mean(c(fitness_males, fitness_hermas), na.rm = TRUE),
    # Divide each fitness vector by the scale factor
    fitness_males_scaled = list(fitness_males / scale_factor),
    fitness_hermas_scaled = list(fitness_hermas / scale_factor),
    mean_male  = mean(fitness_males_scaled[[1]], na.rm = TRUE),
    mean_herma = mean(fitness_hermas_scaled[[1]], na.rm = TRUE),
    ratio      = mean_herma / mean_male
  ) %>%
  ungroup()

fitness_plots_dir <- file.path(plots_dir, "fintess_ratio_plots")
if (!dir.exists(fitness_plots_dir)) {
  dir.create(fitness_plots_dir, recursive = TRUE)
}

######################### unique r
# For each unique s value, create a separate plot
unique_r <- sort(unique(fitness_df$r))
for(r_val in unique_r) {
  df_r <- filter(fitness_df, r == r_val)
  p <- ggplot(df_r, aes(x = factor(s), y = ratio)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.2, size = 2) +
    labs(
      x = "Selfing Rate (s)",
      y = "Mean Male / Mean Hermaphrodite Fitness",
      title = paste("Fitness Ratio at Generation 200 (r =", r_val, ")")
    ) +
    theme_bw(base_size = 14)
  ggsave(
    filename = paste0("Fitness_Ratio_Boxplot_r", r_val, ".png"),
    plot = p,
    path = fitness_plots_dir,
    width = 10,
    height = 6,
    dpi = 300
  )
}

######################### unique s
# For each unique s value, create a separate plot
unique_s <- sort(unique(fitness_df$s))
for(s_val in unique_s) {
  df_s <- filter(fitness_df, s == s_val)
  p <- ggplot(df_s, aes(x = factor(r), y = ratio)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.5) +
    geom_jitter(width = 0.2, size = 2) +
    labs(
      x = "Recombination (r)",
      y = "Mean Male / Mean Hermaphrodite Fitness",
      title = paste("Fitness Ratio at Generation 200 (s =", s_val, ")")
    ) +
    theme_bw(base_size = 14)
  ggsave(
    filename = paste0("Fitness_Ratio_Boxplot_s", s_val, ".png"),
    plot = p,
    path = fitness_plots_dir,
    width = 10,
    height = 6,
    dpi = 300
  )
}
