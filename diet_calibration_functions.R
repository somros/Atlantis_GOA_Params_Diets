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

#' Calibrate Diet Matrix with Multipliers
#' 
#' This function calculates multipliers to adjust diet composition matrices in ecosystem models.
#' It aims to reconcile three matrices: target diet proportions (A), current model inputs (B),
#' and model outputs (C). The function produces multipliers that, when applied to matrix B,
#' will generate new input values that should help the model produce outputs closer to the
#' target proportions while preserving the original magnitude of interactions.
#'
#' @param target_diet Data frame containing target diet proportions (matrix A). Must contain
#'        columns: Predator, PredStage, Prey, Value. Values should sum to 1 for each
#'        Predator-PredStage combination.
#'
#' @param input_diet Data frame containing current model input values (matrix B). Must have
#'        the same structure as target_diet. Values do not need to sum to 1, as they
#'        represent interaction strengths rather than proportions.
#'
#' @param output_diet Data frame containing current model outputs (matrix C). Must have
#'        the same structure as target_diet. Values should sum to 1 for each
#'        Predator-PredStage combination.
#'
#' @param selected_predators Vector of predator codes (e.g., c("ATF", "COD")) indicating
#'        which predators should be calibrated. Predators not in this list will retain
#'        their original values (multiplier = 1). If NULL, all predators are calibrated.
#'
#' @param alpha Numeric value between 0 and 1 controlling the balance between preserving
#'        matrix B's magnitudes (alpha = 0) and achieving matrix A's proportions (alpha = 1).
#'        Default is 0.5.
#'
#' @param tolerance Small numeric value used for comparing floating point numbers and
#'        checking row sums. Default is 1e-6.
#'
#' @return List containing three elements:
#'         - multipliers: Data frame with the same structure as input matrices, containing
#'           multipliers to be applied to input_diet
#'         - diagnostics: Summary statistics and checks for each Predator-PredStage combination
#'         - verification: Data frame showing original values, multipliers, and reconstructed values
#'
#' @details
#' The function operates through several steps:
#' 
#' 1. Input Validation:
#'    - Checks for required columns in all input matrices
#'    - Verifies that target_diet and output_diet sum to 1 by predator-stage
#'    - Warns if row sums are not approximately 1 (within tolerance)
#'
#' 2. Calibration Process:
#'    - Calculates initial adjustments based on how output_diet differs from target_diet
#'    - Applies dampening (0.5) to prevent overshooting
#'    - Uses alpha parameter to balance between:
#'      a) Scaled target proportions (B_scaled = target proportions * original magnitude)
#'      b) Dampened adjustments to current inputs
#'
#' 3. Special Cases:
#'    - PTE prey items: Always retain original values (multiplier = 1)
#'    - Non-selected predators: Retain original values (multiplier = 1)
#'    - Zero values in input_diet:
#'      * If target is also zero: multiplier = 1
#'      * If target is non-zero: multiplier = 999
#'      * If wanting zero: multiplier = 0
#'
#' 4. Value Constraints:
#'    - Caps final input values at 0.99 before calculating multipliers
#'    - Handles zero-division cases appropriately
#'
#' @diagnostics The diagnostics data frame includes:
#'    - max_multiplier: Largest multiplier for each predator-stage
#'    - min_multiplier: Smallest multiplier for each predator-stage
#'    - mean_multiplier: Average multiplier for each predator-stage
#'    - calibrated: Boolean indicating if predator was in selected_predators
#'    - pte_values_present: Boolean indicating if group contains PTE prey
#'    - n_large_multipliers: Count of multipliers > 10
#'    - n_zero_multipliers: Count of zero multipliers
#'
#' @examples
#' # Basic usage with single predator
#' result <- calibrate_diet_matrix_grouped(
#'   target_diet = read.csv("matrixA.csv"),
#'   input_diet = read.csv("matrixB.csv"),
#'   output_diet = read.csv("matrixC.csv"),
#'   selected_predators = c("ATF"),
#'   alpha = 0.5
#' )
#'
#' # Examine results
#' head(result$multipliers)
#' print(result$diagnostics)
#'
#' @notes
#' - The function assumes all three matrices contain the same predator-prey combinations
#' - Multipliers should be carefully reviewed before application, especially large values
#' - The alpha parameter may need tuning based on model behavior
#' - The dampening factor (0.5) is currently hardcoded but could be made into a parameter
#' - The function preserves original values for PTE prey and non-selected predators
#'
#' @warnings
#' - Large multipliers (>10) should be reviewed manually
#' - Zero values in input_diet require special attention
#' - Row sums that don't equal 1 in target_diet or output_diet may indicate issues
calibrate_diet_matrix_grouped <- function(target_diet, input_diet, output_diet, 
                                          selected_predators = NULL,
                                          alpha = 0.5, tolerance = 1e-6) {
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
  
  # Calculate multipliers
  multipliers <- target_diet %>%
    # Join with input_diet to get magnitudes
    left_join(input_diet, 
              by = c("Predator", "PredStage", "Prey"),
              suffix = c("_A", "_B")) %>%
    # Join with output_diet to get current model output
    left_join(output_diet %>% select(Predator, PredStage, Prey, Value_C = Value),
              by = c("Predator", "PredStage", "Prey")) %>%
    group_by(Predator, PredStage) %>%
    mutate(
      # Calculate group magnitudes
      B_magnitude = sum(Value_B, na.rm = TRUE),
      
      # Combine target proportions with input magnitudes
      B_scaled = Value_A * B_magnitude,
      
      # Calculate adjustment based on output vs target difference
      adjustment = case_when(
        Value_C > 0 & Value_A > 0 ~ Value_A / Value_C,
        TRUE ~ 1
      ),
      
      # Apply dampened adjustment
      dampening = 0.5,
      B_prime = Value_B * (1 + dampening * (adjustment - 1)),
      
      # Weighted combination using alpha
      B_final = alpha * B_scaled + (1 - alpha) * B_prime,
      
      # Cap values at 0.99
      B_final = pmin(B_final, 0.99),
      
      # Keep original B values for:
      # 1. Non-selected predators
      # 2. PTE prey items
      B_final = case_when(
        Prey == "PTE" ~ Value_B,
        is.null(selected_predators) ~ B_final,
        Predator %in% selected_predators ~ B_final,
        TRUE ~ Value_B
      ),
      
      # Calculate multiplier
      # For zero B values, set multiplier to 1 if B_final is also 0, otherwise use a large number
      multiplier = case_when(
        abs(Value_B) < tolerance & abs(B_final) < tolerance ~ 1,  # Both zero
        abs(Value_B) < tolerance & B_final > 0 ~ 999,  # B zero but want non-zero
        abs(Value_B) < tolerance ~ 0,  # B zero and want zero
        TRUE ~ B_final / Value_B  # Normal case
      ),
      
      # Ensure multiplier is 1 for preserved cases
      multiplier = case_when(
        Prey == "PTE" ~ 1,
        is.null(selected_predators) ~ multiplier,
        Predator %in% selected_predators ~ multiplier,
        TRUE ~ 1
      )
    ) %>%
    ungroup() %>%
    select(Predator, PredStage, Prey, Value = multiplier)
  
  # Calculate diagnostics
  diagnostics <- multipliers %>%
    left_join(input_diet %>% select(Predator, PredStage, Prey, Value_B = Value),
              by = c("Predator", "PredStage", "Prey")) %>%
    group_by(Predator, PredStage) %>%
    summarise(
      max_multiplier = max(Value, na.rm = TRUE),
      min_multiplier = min(Value, na.rm = TRUE),
      mean_multiplier = mean(Value, na.rm = TRUE),
      calibrated = Predator[1] %in% selected_predators,
      pte_values_present = any(Prey == "PTE"),
      n_large_multipliers = sum(Value > 10, na.rm = TRUE),
      n_zero_multipliers = sum(Value == 0, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Verify multipliers
  verification <- multipliers %>%
    left_join(input_diet, by = c("Predator", "PredStage", "Prey")) %>%
    mutate(
      reconstructed_value = Value.y * Value.x  # original * multiplier
    )
  
  # Return results
  return(list(
    multipliers = multipliers,
    diagnostics = diagnostics,
    verification = verification %>% 
      select(Predator, PredStage, Prey, 
             original = Value.y, 
             multiplier = Value.x, 
             reconstructed = reconstructed_value)
  ))
}