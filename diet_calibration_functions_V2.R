#' Format Diet Matrix from PRM Files
#' 
#' Reads and formats a diet matrix from an Atlantis PRM file, handling both the original
#' and latest diet matrices with options for scaling and age-specific outputs.
#' 
#' @param prm_file Path to the PRM file
#' @param for_order_FG Vector of functional groups in desired order
#' @param scale_factor Numeric value to scale the diet values (default = 1)
#' @param keep_age_specific Logical indicating whether to return age-specific data before grouping (default = FALSE)
#' @param return_names_only Logical indicating whether to only return the names_pprey vector (default = FALSE)
#' 
#' @return If return_names_only = TRUE, returns the names_pprey vector
#'         If keep_age_specific = FALSE, returns a grouped data frame with columns:
#'         Predator, PredStage, Prey, Value
#'         If keep_age_specific = TRUE, returns the age-specific data frame before grouping
#' 
format_diet_matrix <- function(prm_file, for_order_FG, scale_factor = 1, 
                               keep_age_specific = FALSE, return_names_only = FALSE) {
  # Read PRM file
  bio_prm <- readLines(prm_file)
  
  # Identify and index the PPREY matrix
  diets_start <- grep("pPREY1KWT1", bio_prm)
  pprey_ind <- which(startsWith(x = bio_prm, "pPREY") == TRUE)
  diets_end <- max(pprey_ind) + 2
  
  # Extract relevant lines
  names_pprey <- bio_prm[pprey_ind]
  
  # Return early if only names are requested
  if (return_names_only) {
    return(names_pprey)
  }
  
  val_pprey <- bio_prm[pprey_ind + 1]
  
  # Extract consumer groups
  FG <- gsub(" ", "", unique(gsub("pPREY", "", 
                                  gsub('[[:digit:]]+', '', 
                                       gsub("\\   .*", "", names_pprey)))))
  
  # Format diet matrix
  DM_to_format <- t(
    sapply(seq(1, length(val_pprey)),
           function(x) {
             vec <- unlist(strsplit(val_pprey[x], " "))
             return(vec)
           })
  )
  
  # Set column names
  colnames(DM_to_format) <- c(for_order_FG, c("DCsed", "DLsed", "DRsed"))
  
  # Create formatted data frame
  formatted_DM <- DM_to_format %>%
    as_tibble() %>%
    cbind(label = gsub("\\ .*", "", names_pprey)) %>%
    cbind(PredatorCODE = gsub("pPREY", "",
                              gsub('[[:digit:]]+', '', 
                                   gsub("\\ .*", "", names_pprey)))) %>%
    cbind(PreyAgeClass = ifelse(substr(gsub("pPREY", "", 
                                            gsub("\\ .*", "", names_pprey)), 1, 1) %in% c(1, 2),
                                substr(gsub("pPREY", "", 
                                            gsub("\\ .*", "", names_pprey)), 1, 1),
                                "1")) %>%
    cbind(PredatorAgeClass = ifelse(substr(gsub("pPREY", "", 
                                                gsub("\\ .*", "", names_pprey)),
                                           nchar(gsub("pPREY", "", 
                                                      gsub("\\ .*", "", names_pprey))),
                                           nchar(gsub("pPREY", "", 
                                                      gsub("\\ .*", "", names_pprey)))) %in% c(1, 2),
                                    substr(gsub("pPREY", "", 
                                                gsub("\\ .*", "", names_pprey)),
                                           nchar(gsub("pPREY", "", 
                                                      gsub("\\ .*", "", names_pprey))),
                                           nchar(gsub("pPREY", "", 
                                                      gsub("\\ .*", "", names_pprey)))),
                                    "2")) %>%
    mutate(PredatorAgeClass = ifelse(PredatorAgeClass == 1, "Juvenile", "Adult"),
           PreyAgeClass = ifelse(PreyAgeClass == 2, "Adult", "Juvenile")) %>%
    dplyr::select(c("label", "PredatorCODE", "PreyAgeClass", 
                    "PredatorAgeClass", colnames(DM_to_format)))
  
  # Convert to long format with age-specific data
  formatted_DM_long_age <- formatted_DM %>%
    pivot_longer(5:ncol(formatted_DM), 
                 values_to = "diet_value", 
                 names_to = "Prey") %>%
    mutate(diet_value = as.numeric(diet_value) * scale_factor) %>%
    ungroup()
  
  # Return age-specific data if requested
  if (keep_age_specific) {
    return(formatted_DM_long_age)
  }
  
  # Otherwise, collapse prey stages and return grouped data
  formatted_DM_long <- formatted_DM_long_age %>%
    group_by(PredatorCODE, PredatorAgeClass, Prey) %>%
    summarize(diet_value = mean(diet_value, na.rm = TRUE)) %>%
    ungroup() %>%
    select(PredatorCODE, PredatorAgeClass, Prey, diet_value) %>%
    rename(
      Predator = PredatorCODE,
      PredStage = PredatorAgeClass,
      Value = diet_value
    )
  
  return(formatted_DM_long)
}


# target_diet = formatted_DM_original_long
# input_diet = formatted_DM_latest_long
# output_diet = formatted_DM_pred_long
# selected_predators = tier3
# alpha = my_alpha
# tolerance = 1e-6
# dampening = 0.5
# small_threshold = 0.01

#' Diet Matrix Calibration with Small Value Handling
#' 
#' Calibrates diet composition matrices in ecosystem models using a two-step approach:
#' first adjusting for proportions while handling small values appropriately, then
#' scaling for magnitude. Returns multipliers to be applied to the input matrix.
#'
#' @param target_diet Data frame of target diet proportions (matrix A). Must contain columns:
#'        Predator, PredStage, Prey, Value. Values must sum to 1 for each Predator-PredStage.
#'
#' @param input_diet Data frame of current model inputs (matrix B). Same structure as
#'        target_diet. Values do not sum to 1 as they represent interaction strengths.
#'
#' @param output_diet Data frame of current model outputs (matrix C). Same structure as
#'        target_diet. Values must sum to 1 for each Predator-PredStage.
#'
#' @param selected_predators Vector of predator codes to calibrate (e.g., c("ATF", "COD")).
#'        Other predators retain original values. If NULL, calibrates all predators.
#'
#' @param alpha Numeric [0,1] controlling magnitude preservation. alpha = 0 applies no
#'        magnitude scaling, alpha = 1 fully restores original magnitudes. Default 0.5.
#'
#' @param dampening Numeric [0,1] controlling adjustment aggressiveness. Lower values
#'        make more conservative adjustments. Default 0.5.
#'
#' @param small_threshold Numeric value. Prey items with proportions below this in
#'        target_diet are not calibrated. Default 0.001.
#'
#' @param tolerance Numeric value for floating point comparisons. Default 1e-6.
#'
#' @return List with three elements:
#'         multipliers: Data frame of multipliers to apply to input_diet
#'         diagnostics: Summary statistics by Predator-PredStage
#'         verification: Detailed tracking of adjustment process
#'
#' @details
#' Step 1 - Proportion Adjustment:
#' - Calculates how model output differs from target
#' - Scales adjustments by prey importance in target diet
#' - Skips adjustment for very small prey items
#' - Applies dampening to prevent overshooting
#'
#' Step 2 - Magnitude Scaling:
#' - Scales values to help preserve original magnitudes
#' - Controlled by alpha parameter
#' - Maintains proportional relationships from step 1
#'
#' Special cases:
#' - PTE prey items retain original values
#' - Non-selected predators retain original values
#' - Very small prey items (<small_threshold) retain original values
#' - Values capped at 0.99
#' - Zero values handled specially
#'
#' @examples
#' result <- calibrate_diet_matrix_grouped(
#'   target_diet = read.csv("matrixA.csv"),
#'   input_diet = read.csv("matrixB.csv"),
#'   output_diet = read.csv("matrixC.csv"),
#'   selected_predators = c("ATF"),
#'   alpha = 0.5,
#'   dampening = 0.5,
#'   small_threshold = 0.001
#' )

calibrate_diet_matrix_grouped <- function(target_diet, input_diet, output_diet, 
                                          selected_predators = NULL,
                                          alpha = 0.5, 
                                          dampening = 0.5,
                                          small_threshold = 0.01,
                                          tolerance = 1e-6) {
  # Validate input data
  validate_data <- function(df, name) {
    required_cols <- c("Predator", "PredStage", "Prey", "Value")
    if (!all(required_cols %in% names(df))) {
      stop(paste("Missing required columns in", name))
    }
  }
  
  validate_data(target_diet, "target_diet")
  validate_data(input_diet, "input_diet")
  validate_data(output_diet, "output_diet")
  
  # Verify row sums for target and output diets
  check_proportions <- function(df, name) {
    sums <- df %>%
      group_by(Predator, PredStage) %>%
      summarise(total = sum(Value), .groups = 'drop')
    
    if (!all(abs(sums$total - 1) < tolerance)) {
      warning(paste("Not all groups in", name, "sum to 1"))
      print(filter(sums, abs(total - 1) >= tolerance))
    }
  }
  
  check_proportions(target_diet, "target_diet")
  check_proportions(output_diet, "output_diet")
  
  # Step 1: Calculate proportion-focused adjustments with small value handling
  B_proportional <- target_diet %>%
    left_join(input_diet, 
              by = c("Predator", "PredStage", "Prey"),
              suffix = c("_A", "_B")) %>%
    left_join(output_diet %>% select(Predator, PredStage, Prey, Value_C = Value),
              by = c("Predator", "PredStage", "Prey")) %>%
    group_by(Predator, PredStage) %>%
    mutate(
      # Calculate relative importance of each prey in target diet
      relative_importance = Value_A / max(Value_A),
      
      # Calculate basic proportion adjustment
      raw_adjustment = case_when(
        Value_C > 0 & Value_A > 0 ~ Value_A / Value_C,
        TRUE ~ 1
      ),
      
      # Apply adjustments based on importance and threshold
      proportion_adjustment = case_when(
        # Ignore very small prey items in target diet
        Value_A < small_threshold ~ 1,
        # For significant prey items, scale adjustment by relative importance
        TRUE ~ 1 + (raw_adjustment - 1) * relative_importance
      ),
      
      # Apply dampening to the final adjustment
      B_adjusted = Value_B * (dampening * proportion_adjustment),
      
      # Store original B magnitudes for later
      B_magnitude = sum(Value_B, na.rm = TRUE)
    ) %>%
    ungroup()
  
  # Step 2: Scale magnitudes while preserving proportions
  B_final <- B_proportional %>%
    group_by(Predator, PredStage) %>%
    mutate(
      # Get the current magnitude after proportion adjustment
      new_magnitude = sum(B_adjusted, na.rm = TRUE),
      
      # Scale according to original magnitude, controlled by alpha
      magnitude_scalar = alpha * (B_magnitude / new_magnitude) + (1 - alpha),
      
      # Apply magnitude scaling
      B_final = B_adjusted * magnitude_scalar,
      
      # Cap values at 0.99
      B_final = pmin(B_final, 0.99),
      
      # Handle special cases (PTE and non-selected predators)
      B_final = case_when(
        Prey == "PTE" ~ Value_B,
        is.null(selected_predators) ~ B_final,
        Predator %in% selected_predators ~ B_final,
        TRUE ~ Value_B
      ),
      
      # Calculate multipliers
      multiplier = case_when(
        abs(Value_B) < tolerance & abs(B_final) < tolerance ~ 1,
        abs(Value_B) < tolerance & B_final > 0 ~ 999,
        abs(Value_B) < tolerance ~ 0,
        TRUE ~ B_final / Value_B
      )
    ) %>%
    ungroup()
  
  # Calculate diagnostics
  diagnostics <- B_final %>%
    group_by(Predator, PredStage) %>%
    summarise(
      proportion_adjustment_mean = mean(proportion_adjustment, na.rm = TRUE),
      proportion_adjustment_max = max(proportion_adjustment, na.rm = TRUE),
      original_magnitude = first(B_magnitude),
      final_magnitude = sum(B_final, na.rm = TRUE),
      magnitude_ratio = final_magnitude / original_magnitude,
      n_small_values = sum(Value_A < small_threshold),
      mean_relative_importance = mean(relative_importance, na.rm = TRUE),
      max_raw_adjustment = max(raw_adjustment, na.rm = TRUE),
      max_final_adjustment = max(proportion_adjustment, na.rm = TRUE),
      max_multiplier = max(multiplier, na.rm = TRUE),
      min_multiplier = min(multiplier, na.rm = TRUE),
      calibrated = Predator[1] %in% selected_predators,
      pte_values_present = any(Prey == "PTE"),
      .groups = 'drop'
    )
  
  # Return results with verification data
  return(list(
    multipliers = B_final %>% 
      select(Predator, PredStage, Prey, Value = multiplier),
    diagnostics = diagnostics,
    verification = B_final %>% 
      select(Predator, PredStage, Prey,
             original = Value_B,
             target_prop = Value_A,
             relative_importance,
             raw_adjustment,
             proportion_adjustment,
             final = B_final,
             multiplier)
  ))
}