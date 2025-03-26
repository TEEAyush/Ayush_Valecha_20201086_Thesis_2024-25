
###############################################################################
# function to  parse the SLiM output
###############################################################################

parse_slim2_output <- function(slim_lines) {
  # slim_lines should be a character vector containing lines from SLiM's output
  
  n <- length(slim_lines)
  i <- 1
  
  parsed_data <- list()
  
  while (i <= n) {
    line_i <- slim_lines[i]
    
    # Look for lines starting with "generation <number>"
    if (str_detect(line_i, "^generation\\s+\\d+")) {
      # Extract the generation number
      gen_num <- as.integer(str_extract(line_i, "\\d+"))
      
      # Move to the next line, which should contain the polymorphic count (or some integer)
      i <- i + 1
      poly_count <- if (i <= n) as.integer(slim_lines[i]) else NA
      
      # Initialize placeholders for the data we want to capture
      heterozyg_val <- NA_real_
      fitness_males_vec  <- NA
      fitness_hermas_vec <- NA
      position_poly_loci <- NA
      freqs_poly_loci    <- NA
      poly_male_sels     <- NA
      poly_herma_sels    <- NA
      freq_1_fixed_loci  <- NA
      freq_1_male_sels   <- NA
      freq_1_herma_sels  <- NA
      freq_0_fixed_loci  <- NA
      freq_0_male_sels   <- NA
      freq_0_herma_sels  <- NA
      
      # # We'll store the per‐genome variant positions in a list
      # # e.g. genome_positions[[ genome_idx ]] = c(1,2,10, ...)
      # genome_positions <- list()
      
      
      i <- i + 1
      done_block <- FALSE
      
      # Read additional lines until we encounter the next "generation"
      while (!done_block && i <= n) {
        tmp_line <- slim_lines[i]
        
        if (str_detect(tmp_line, "^generation\\s+\\d+")) {
          # We have reached the next generation block, so we stop accumulating data.
          done_block <- TRUE
        } else {
          # Now parse the recognized patterns
          if (str_detect(tmp_line, "^fitness of males:")) {
            # e.g. "fitness of males: 0.664878,0.664878,0.697127,..."
            content_str <- str_replace(tmp_line, "^fitness of males:\\s*", "")
            fitness_males_vec <- as.numeric(str_split(content_str, ",")[[1]])
            
          } else if (str_detect(tmp_line, "^fitness of hermas:")) {
            # e.g. "fitness of hermas: 2.22726,2.31537, ..."
            content_str <- str_replace(tmp_line, "^fitness of hermas:\\s*", "")
            fitness_hermas_vec <- as.numeric(str_split(content_str, ",")[[1]])
            
          } else if (str_detect(tmp_line, "^positions of polymorphic loci:")) {
            # e.g. "positions of polymorphic loci: 3, 1"
            content_str <- str_replace(tmp_line, "^positions of polymorphic loci:\\s*", "")
            position_poly_loci <- as.numeric(str_split(content_str, ",")[[1]])
            
          } else if (str_detect(tmp_line, "^frequencies of polymorphic loci:")) {
            # e.g. "frequencies of polymorphic loci: 0.835,0.945"
            content_str <- str_replace(tmp_line, "^frequencies of polymorphic loci:\\s*", "")
            freqs_poly_loci <- as.numeric(str_split(content_str, ",")[[1]])
            
          } else if (str_detect(tmp_line, "^polymorphic positions male selection coefficients:")) {
            content_str <- str_replace(tmp_line, "^polymorphic positions male selection coefficients:\\s*", "")
            poly_male_sels <- as.numeric(str_split(content_str, ",")[[1]])
            
          } else if (str_detect(tmp_line, "^polymorphic positions herma selection coefficients:")) {
            content_str <- str_replace(tmp_line, "^polymorphic positions herma selection coefficients:\\s*", "")
            poly_herma_sels <- as.numeric(str_split(content_str, ",")[[1]])
            
          } else if (str_detect(tmp_line, "^positions of freq 1 loci:")) {
            content_str <- str_replace(tmp_line, "^positions of freq 1 loci:\\s*", "")
            freq_1_fixed_loci <- as.numeric(str_split(content_str, ",")[[1]])
          
          } else if (str_detect(tmp_line, "^Freq 1 positions male selection coefficients:")) {
            content_str <- str_replace(tmp_line, "^Freq 1 positions male selection coefficients:\\s*", "")
            freq_1_male_sels <- as.numeric(str_split(content_str, ",")[[1]])
            
          } else if (str_detect(tmp_line, "^Freq 1 positions herma selection coefficients:")) {
            content_str <- str_replace(tmp_line, "^Freq 1 positions herma selection coefficients:\\s*", "")
            freq_1_herma_sels <- as.numeric(str_split(content_str, ",")[[1]])
         
          } else if (str_detect(tmp_line, "^positions of freq 0 loci:")) {
            content_str <- str_replace(tmp_line, "^positions of freq 0 loci:\\s*", "")
            freq_0_fixed_loci <- as.numeric(str_split(content_str, ",")[[1]])
            
          } else if (str_detect(tmp_line, "^Freq 0 positions male selection coefficients:")) {
            content_str <- str_replace(tmp_line, "^Freq 0 positions male selection coefficients:\\s*", "")
            freq_0_male_sels <- as.numeric(str_split(content_str, ",")[[1]])
            
          } else if (str_detect(tmp_line, "^Freq 0 positions herma selection coefficients:")) {
            content_str <- str_replace(tmp_line, "^Freq 0 positions herma selection coefficients:\\s*", "")
            freq_0_herma_sels <- as.numeric(str_split(content_str, ",")[[1]])
            
          }  else if (str_detect(tmp_line, "^Heterozygosity\\s+")) {
            # e.g. "Heterozygosity 0.5243"
            content_str <- str_replace(tmp_line, "^Heterozygosity\\s+", "")
            heterozyg_val <- as.numeric(content_str)
          }
          
          
#          else if (str_detect(tmp_line, "^Genome\\s+\\d+:")) {
#            # Lines like: "Genome 0: 0,1,4,5"
#            # 1) Extract the genome index
#            genome_idx <- as.integer(str_extract(tmp_line, "(?<=Genome\\s)\\d+"))
#            # 2) Extract the positions
#            pos_str <- str_replace(tmp_line, "^Genome\\s+\\d+:\\s*", "")
#            if (str_trim(pos_str) == "") {
#              # in case there's no positions, means no variants in that genome
#              genome_positions[[genome_idx + 1]] <- integer(0)
#            } else {
#              genome_positions[[genome_idx + 1]] <- as.numeric(str_split(pos_str, ",")[[1]])
#            }
#            
#          }
          
          i <- i + 1
        }
      }
      
      # # ----------------------------------------------------------------------
      # # Build the haplotype matrix from the 'genome_positions' list
      # # 1) Figure out how many "Genome i" lines we have => # of rows
      # # 2) Figure out the max position => # of columns (positions go from 0..maxPos)
      # # ----------------------------------------------------------------------
      # if (length(genome_positions) > 0) {
      #   max_pos <- 0
      #   for (gp in genome_positions) {
      #     if (length(gp) > 0) {
      #       max_pos <- max(max_pos, max(gp))
      #     }
      #   }
      #   # matrix has (number_of_genomes) rows, columns from 0..max_pos => (max_pos+1) columns
      #   hap_mat <- matrix(0, nrow = length(genome_positions), ncol = max_pos + 1)
      #   
      #   # fill it in
      #   # genome_positions[[row_i]] is a vector of positions
      #   for (row_i in seq_along(genome_positions)) {
      #     positions_here <- genome_positions[[row_i]]
      #     if (length(positions_here) > 0) {
      #       hap_mat[row_i, positions_here + 1] <- 1  # +1 since R is 1-based
      #     }
      #   }
      # } else {
      #   # if no Genome lines for some reason
      #   hap_mat <- matrix(, nrow=0, ncol=0)
      # }
      # 
      
      # Record the parsed information for this generation
      # Some of these will be vectors, so we store them as list-columns
      parsed_data[[length(parsed_data) + 1]] <-
        tibble(
          generation           = gen_num,
          polymorphic_count    = poly_count,
          fitness_males        = list(fitness_males_vec/mean(c(fitness_males_vec,fitness_hermas_vec))),
          fitness_hermas       = list(fitness_hermas_vec/mean(c(fitness_males_vec,fitness_hermas_vec))),
          freq_poly_loci       = list(freqs_poly_loci),
          position_poly_loci   = list(position_poly_loci),
          poly_pos_male_sels   = list(poly_male_sels),
          poly_pos_herma_sels  = list(poly_herma_sels),
          freq_1_fixed_loci    = list(freq_1_fixed_loci),
          freq_1_male_sels     = list(freq_1_male_sels),
          freq_1_herma_sels    = list(freq_1_herma_sels),
          freq_0_fixed_loci    = list(freq_0_fixed_loci),
          freq_0_male_sels     = list(freq_0_male_sels),
          freq_0_herma_sels    = list(freq_0_herma_sels),
          heterozygosity       = heterozyg_val
          # fixed_pos_male_sels  = list(fixed_male_sels),
          # fixed_pos_herma_sels = list(fixed_herma_sels),
          # haplotype_matrix     = list(hap_mat)
        )
    } else {
      i <- i + 1
    }
  }
  
  # Combine all parsed data into a single data frame
  final_df <- bind_rows(parsed_data)
  return(final_df)
}


#################################################
#function to run slim 
#################################################


run_slim_simulation2 <- function(
    s,
    r,
    K,
    L,
    mut_init_freq 
){
  # 1) Build the SLiM command (adjust path to your own script)
  slimcmd <- paste(
    "slim",
    paste0("-d K=", K),
    paste0("-d s=", s),
    paste0("-d L=", L),
    paste0("-d mut_init_freq=", mut_init_freq),
    paste0("-d r=", r),
    # Example: new SLiM script path
    "/mnt/data2/avalecha/slim/SLiM/My_scripts/Andro_4/Androdioecy_4_lab_pc.slim"
  )
  
  # 2) Run SLiM and capture output (run in R, capturing console output as a character vector)
  slimoutput <- system(slimcmd, intern = TRUE)
  
  # 3) Parse the output using the new parser
  results_df <- parse_slim2_output(slimoutput)
  
  # 4) Optionally store the simulation parameters for reference
  results_df <- results_df %>%
    mutate(
      s = s,
      r = r,
      K = K,
      L = L,
      mut_init_freq = mut_init_freq
    )
  
  return(results_df)
}


# ############################################## functions for analysis of haplotype matrix to plot the LD vs chromosome position plots
# 
# ### 1.  A helper function to compute pairwise LD (D_{i,j}) between two columns/loci
# 
# compute_D_ij <- function(col_i, col_j) {
#   # col_i, col_j are 0/1 vectors, each length = number_of_haplotypes
#   p_i  <- mean(col_i)  # frequency of "1" at locus i
#   p_j  <- mean(col_j)  # frequency of "1" at locus j
#   p_ij <- mean(col_i & col_j)  # frequency of haplotypes that have "1" at both i and j
#   D_ij <- p_ij - (p_i * p_j)
#   return(D_ij)
# }
# 
# 
# ### 2. A function to compute average LD in a window
# 
# compute_LD_in_window <- function(window_submatrix) { # window_submatrix is a submatrix of the hap_matrix
#   # poly_indices = vector of polymorphic indices in the window
#   # Return average D_ij across all distinct pairs
#   col_is_poly <- apply(window_submatrix, 2, function(x) {
#     any(x == 1) && any(x == 0)
#   })
#   
#   if (sum(col_is_poly) < 2) {
#     return(NA_real_)  # can't compute pairwise LD with <2 loci
#   }
#   
#   # Extract just those columns
#   window_submatrix_poly <- window_submatrix[, col_is_poly, drop = FALSE]
#   n_loci <- length(col_is_poly)
#   
#   sum_D  <- 0
#   count  <- 0
#   for (a in 1:(n_loci-1)) {
#     for (b in (a+1):n_loci) {
#       D_val <- compute_D_ij(window_submatrix_poly[, a], window_submatrix_poly[, b])
#       sum_D <- sum_D + D_val
#       count <- count + 1
#     }
#   }
#   return(sum_D / count)
# }
# 
# 
# ### 3. The main sliding‐window function
# sliding_window_LD <- function(hap_mat, window_size = 20, step = 1) {
#   min_pos <- 1
#   max_pos <- dim(hap_mat)[[2]] - window_size + 1
#   
#   window_starts <- seq(min_pos, max_pos, by = step)
#   
#   for (ws in window_starts) {
#     window_submatrix <- hap_mat[,ws:ws + window_size]
#     mean_D <- compute_LD_in_window(window_submatrix)
#     window_center <- ws + (window_size / 2)
#     results_list[[length(results_list)+1]] <- tibble(
#       window_start  = ws,
#       window_end    = ws + window_size,
#       window_center = window_center,
#       mean_LD       = mean_D
#       )
#   }
#   return(bind_rows(results_list))
#   
# }

  
################################# Sliding window analysis (calculate male_sel*freq and male_sel*freq )


################ calculates these values within a window
calc_window_avgs <- function(freq_vec, male_sel_vec, herma_sel_vec) {
  # freq_vec, male_sel_vec, herma_sel_vec are numeric vectors of the same length (the sites in the window)
  # Returns a named list with:
  #   avg_male  = mean( male_sel * freq )
  #   avg_herma = mean( herma_sel * freq )
  
  # Filter out NA or mismatch lengths if needed
  stopifnot(
    length(freq_vec) == length(male_sel_vec),
    length(freq_vec) == length(herma_sel_vec)
  )
  
  avg_male  <- mean(male_sel_vec  * freq_vec, na.rm = TRUE)
  avg_herma <- mean(herma_sel_vec * freq_vec, na.rm = TRUE)
  
  list(avg_male = avg_male, avg_herma = avg_herma)
}



sliding_window_sel_freq <- function(positions, freq, male_sel, herma_sel,
                                    window_size = 20, step = 1)
  
  # positions: sorted numeric vector of locus positions (e.g. 0..L-1)
  # freq: a numeric vector, same length as positions (the frequencies for each locus)
  # male_sel, herma_sel: numeric vectors, same length as positions (the selection coefficients)
  #
  # window_size: size of the window in base pairs
  # step: how much to slide the window each time
  #
  # Output: a data frame with columns:
  #    window_center, avg_male, avg_herma
{
  stopifnot(
    length(positions) == length(freq),
    length(positions) == length(male_sel),
    length(positions) == length(herma_sel)
  )
  
  # Make sure data are sorted by position
  idx <- order(positions)
  positions <- positions[idx]
  freq      <- freq[idx]
  male_sel  <- male_sel[idx]
  herma_sel <- herma_sel[idx]
  
  results_list <- list()
  
  min_pos <- min(positions, na.rm = TRUE)
  max_pos <- max(positions, na.rm = TRUE)
  
  # Window "start" from min_pos up to (max_pos - window_size)
  window_starts <- seq(min_pos, max_pos - window_size, by = step)
  
  for (ws in window_starts) {
    we <- ws + window_size
    in_window_idx <- which(positions >= ws & positions <= we)
    
    if (length(in_window_idx) == 0) {
      # No loci in this window
      results_list[[length(results_list) + 1]] <- data.frame(
        window_center = ws + window_size/2,
        avg_male      = NA_real_,
        avg_herma     = NA_real_
      )
    } else {
      freq_sub  <- freq[in_window_idx]
      male_sub  <- male_sel[in_window_idx]
      herma_sub <- herma_sel[in_window_idx]
      
      avgs <- calc_window_avgs(freq_sub, male_sub, herma_sub)
      
      results_list[[length(results_list) + 1]] <- data.frame(
        window_center = ws + window_size/2,
        avg_male      = avgs$avg_male,
        avg_herma     = avgs$avg_herma
      )
    }
  }
  
  do.call(rbind, results_list)
}


# 2) A helper function for the deterministic formula
#    Using the closed-form:
#      H(t) = alpha + (S/2)^(t - 2) * (H2 - alpha)
#    where alpha = [2(1-S)/(2-S)], and H2 is the actual Het at generation=2.
det_het <- function(gen, S, H2, gen2=2) {
  alpha <- 2 * (1 - S) / (2 - S)
  alpha + (S/2)^(gen - gen2) * (H2 - alpha)
}

