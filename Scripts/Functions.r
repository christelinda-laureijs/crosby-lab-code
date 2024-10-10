
# Import and modify data ----

add_new_cells <- function(new_raw_data_csv,
                          cell_characteristics_csv,
                          old_raw_data_csv,
                          current_type) {
  # Obtain argument values as strings
  # Required to check if the filenames and current_type match (e.g. User enters raw-sEPSC-data.csv for current_type = "sEPSC")
  
  new_raw_data_name <- deparse(substitute(new_raw_data_csv))
  cell_characteristics_csv_name <- deparse(substitute(cell_characteristics_csv))
  old_raw_data_csv_name <- deparse(substitute(old_raw_data_csv))
  current_type_name <- deparse(substitute(current_type))
  
  list_of_argument_names <- c(
    new_raw_data_name,
    cell_characteristics_csv_name,
    old_raw_data_csv_name,
    current_type_name
  )
  
  if (is.null(new_raw_data_csv) ||
      !is.character(new_raw_data_csv)) {
    stop(
      "'new_raw_data_csv' must be a character (e.g. \"Data/Raw-CSVs/eEPSC-data/20240912-raw-data.csv\")"
    )
  }
  
  if (is.null(cell_characteristics_csv) ||
      !is.character(cell_characteristics_csv)) {
    stop(
      "'cell_characteristics_csv' must be a character (e.g. \"Data/Plaintext-Cell-Characteristics.csv\")"
    )
  }
  
  if (is.null(current_type) ||
      length(current_type) != 1L ||
      !current_type %in% c("eEPSC", "sEPSC")) {
    stop("'current_type' argument must be one of: 'eEPSC' or 'sEPSC'")
  }
  
  # Required to see if new cells have associated data for synapses, treatment, sex, age, etc.
  cell_characteristics <-
    read.csv(here(cell_characteristics_csv)) %>%
    rename_with(tolower) %>%
    rename(X = x, Y = y)
  
  new_raw_data <- read.csv(here(new_raw_data_csv))
  
  new_raw_data <- new_raw_data %>% 
    rename_with(tolower) %>% 
    mutate(id = factor(id))
  
  if (current_type == "eEPSC") {
    if (any(grepl("sEPSC", list_of_argument_names))) {
      stop(
        "current_type = \"",
        current_type,
        "\" but some filenames have ",
        "\"sEPSC\".",
        "\n Are you sure that you selected the correct current type?"
      )
    }
    new_raw_data <- new_raw_data %>%
      rename_with(tolower) %>%
      rename(ID = id, P1 = p1, P2 = p2) %>%
      group_by(letter) %>%
      mutate(time = (row_number() - 1) / 12)
  }
  
  if (current_type == "sEPSC") {
    if (any(grepl("eEPSC", list_of_argument_names))) {
      stop(
        "current_type = \"",
        current_type,
        "\" but some filenames have ",
        "\"eEPSC\".",
        "\n Are you sure that you selected the correct current type?"
      )
    }
    
    new_raw_data <- new_raw_data %>%
      rename_with(tolower) %>%
      rename(ID = id) %>%
      group_by(letter) %>%
      mutate(amplitude = (-1) * amplitude,
             time = ((recording_num - 1) * 300 + (trace - 1) * 5 + (time_of_peak /
                                                                      1000)) / 60)
  }
  
  warning("Renamed dataframe columns to lowercase")
  
  new_letters <- unique(new_raw_data$letter)
  new_letters_spaces <- paste0(new_letters, " ")
  
  if (all(new_letters %in% cell_characteristics$letter)) {
    message(
      "Letter check passed: All letters are present in both \"",
      cell_characteristics_csv,
      "\" \nand \"",
      new_raw_data_csv,
      "\".",
      "\n \nThe matching cells are:",
      "\n",
      new_letters_spaces
    )
    
  } else {
    stop(
      "Missing letters detected in \"",
      cell_characteristics_csv,
      "\". \nDid you add cell characteristics for ALL the new cells in \"",
      new_raw_data_csv,
      "\"?"
    )
  }
  
  new_raw_data_complete <-
    merge(new_raw_data, cell_characteristics, by = "letter")
  
  # Import old Raw-Data sheet
  old_raw_data <- read.csv(here(old_raw_data_csv), header = T)
  
  old_raw_data <- old_raw_data %>% 
    rename_with(tolower) %>% 
    mutate(id = factor(id))
  
  if (current_type == "eEPSC") {
    old_raw_data <- old_raw_data %>%
      rename_with(tolower) %>%
      rename(
        ID = id,
        P1 = p1,
        P2 = p2,
        X = x,
        Y = y
      )
  }
  
  if (current_type == "sEPSC") {
    old_raw_data <- old_raw_data %>%
      rename_with(tolower) %>%
      rename(ID = id, X = x, Y = y)
  }
  
  if (any(grepl("cells", colnames(old_raw_data)))) {
    old_raw_data <- old_raw_data %>%
      rename(cell = cells)
    
    warning("Renamed column 'cells' to 'cell'")
  }
  
  letters_in_old_raw_data <- unique(old_raw_data$letter)
  letters_in_new_raw_data <- unique(new_raw_data_complete$letter)
  letters_in_new_raw_data_spaces <- paste0(letters_in_new_raw_data, " ")
  
  # Check for accidental letter duplication (e.g. running same code twice)
  if (any(letters_in_new_raw_data %in% letters_in_old_raw_data)) {
    overlapping_letters <- intersect(letters_in_old_raw_data, letters_in_new_raw_data)
    overlapping_letters_spaces <- paste0(overlapping_letters, " ")
    
    stop(
      "\n!! Overlapping letter(s) detected!!    --->  ",
      overlapping_letters_spaces,
      "\nThese cells in \"",
      new_raw_data_csv,
      "\" already exist in ",
      "\"",
      old_raw_data_csv,
      "\"",
      "\n \nDuplicate letters will produce plotting errors and distorted statistics.",
      "\nDid you accidentally run this code twice?"
    )
    
  } else {
    message(
      "\n \nLetter check passed: All letters in \"",
      new_raw_data_csv,
      "\" are new relative to ",
      "\"",
      old_raw_data_csv,
      "\"",
      "\n \nAdding the new following new cells:",
      "\n",
      letters_in_new_raw_data_spaces
    )
    
    file.rename(from = here(old_raw_data_csv),
                to = here(paste0(
                  str_sub(old_raw_data_csv, 1, -5), "-old.csv"
                )))
    
    full_raw_data <- bind_rows(old_raw_data, new_raw_data_complete)
    
    write.csv(full_raw_data, here(old_raw_data_csv), row.names = F)
    
    file.remove(here(paste0(
      str_sub(old_raw_data_csv, 1, -5), "-old.csv"
    )))
  }
}

#' Normalize raw current data

#' `make_normalized_EPSC_df()` creates a dataframe of evoked or spontaneous current data from a raw .csv file.
#' 
#' The function will create a new column containing the evoked/spontaneous current amplitudes
#' normalized relative to the mean baseline amplitude.
#' For evoked current data, the function also adds a column for the paired pulse ratio (P2/P1, where P2 is the amplitude of the second evoked current)
#' 
#' @param filename A filepath to a .csv file. See [add_new_cells()] for the function that will merge raw data (a .csv with 4 columns: `letter`, `ID`, `P1`, and `P2`)
#' and a `cell-characteristics.csv` file (with columns for factors like `animal`, `age`, `sex`, `synapses`).
#' 
#' If the data are evoked currents (current_type == "eEPSC"), the .csv file must contain the following columns and data types:
#' 
#' \itemize{
#'  \item `letter` A character value that is a unique identifier for a single recording. Used to link data sets for evoked or spontaneous currents and cell-characteristics.
#'  \item `synapses` A character value (e.g. "Glutamate" or "GABA").
#'  \item `sex` A character value (e.g. "Male" or "Female").
#'  \item `treatment` A character value (e.g. "Control", "HNMPA").
#'  \item `time` A numeric value that represents time in minutes. This column is autogenerated in [add_new_cells()].
#'  \item `ID` A character value for the recording filename.
#'  \item `P1` A numeric value representing the amplitude of the first evoked current in pA.
#'  \item `P2` A numeric value representing the amplitude of the second evoked current in pA.
#'  \item `X` A numeric value representing the x-value of the cell's location in µm.
#'  \item `Y` A numeric value representing the y-value of the cell's location in µm.
#'  \item `age` A numeric value representing the animal's age. Can be any value as long as the time units are consistent throughout (e.g. don't mix up days and months when reporting animal ages).
#'  \item `animal` A numeric value representing the animal's ID or number.
#'  \item `category` A numeric value representing the experiment type. Used to assign top-level groups for further analyses, with `treatment` as subgroups.
#'  \item `cell` A character or numeric value representing the cell. For example, use `3.1.1` for animal #3, slice #1, cell #1.
#' }
#' If the data are spontaneous currents (current_type == "sEPSC"), the .csv file must contain the following columns and data types:
#' 
#' \itemize{
#'  \item `letter` A character value that is a unique identifier for a single recording. Used to link data sets for evoked or spontaneous currents and cell-characteristics.
#'  \item `synapses` A character value (e.g. "Glutamate" or "GABA").
#'  \item `sex` A character value (e.g. "Male" or "Female").
#'  \item `treatment` A character value (e.g. "Control", "HNMPA").
#'  \item `time` A numeric value that represents time in minutes. This column is autogenerated in [add_new_cells()].
#'  \item `ID` A character value for the recording filename.
#'  \item `recording_num` A numeric value representing the recording number. This was incorporated before we switched to concatenating all recordings into one, but it needs to remain here to prevent breaking previous projects. It should be set to 1.
#'  \item `time_of_peak` A numeric value representing the time of the peak in milliseconds relative to trace number, which is in turn relative to the start of the recording. See [add_new_cells()] for a description of how the true time value (`time`) is calculated from the `recording_num` and `trace.`
#'  \item `amplitude` A numeric value representing the amplitude of the evoked current in pA.
#'  \item `X` A numeric value representing the x-value of the cell's location in µm.
#'  \item `Y` A numeric value representing the y-value of the cell's location in µm.
#'  \item `age` A numeric value representing the animal's age. Can be any value as long as the time units are consistent throughout (e.g. don't mix up days and months when reporting animal ages).
#'  \item `animal` A numeric value representing the animal's ID or number.
#'  \item `category` A numeric value representing the experiment type. Used to assign top-level groups for further analyses, with `treatment` as subgroups.
#'  \item `cell` A character or numeric value representing the cell. For example, use `3.1.1` for animal #3, slice #1, cell #1.
#' }
#' 
#' 
#' @param current_type A character describing the current type. Allowed values are "eEPSC" or "sEPSC".
#' @param min_time_value Minimum time value (in minutes), which defaults to 0.
#' @param max_time_value Maximum recording length (in minutes). All data points will be filtered to time values less than or equal to this value. Defaults to 25.
#' @param baseline_length Length of the baseline (in minutes). Refers to data collected before applying a hormone, antagonist, or a protocol like high frequency stimulation. Defaults to 5.
#' @param interval_length Length of each interval (in minutes). Used to divide the dataset into broad ranges for statistical analysis. Important! `max_recording_length` must be divisible by `interval_length`.
#' Defaults to 5. 
#' @param negative_transform_currents A character ("yes" or "no") describing if `P1` and `P2` should be negative transformed. If "yes", the values will be multiplied by (-1).
#' @returns A dataframe that can be viewed and used for further analyses in R. New or modified columns include:
#' 
#' \itemize{
#'  \item `P1` (for evoked currents only) May be negative-transformed if `negative_transform` == "yes"
#'  \item `P2` (for evoked currents only) May be negative-transformed if `negative_transform` == "yes"
#'  \item `PPR` (for evoked currents only) A numeric value that represents the paired pulse ratio (PPR) of the evoked currents, generated using `mutate(PPR = P2/P1)`.
#'  \item `baseline_range` A logical value required for the baseline transformation. It is set to TRUE when time is within the baseline period (e.g. Time <= 5) and FALSE at all other times.
#'  \item `baseline_mean` A numeric value representing the mean evoked current amplitude during the baseline period. There is a different baseline_mean for each letter.
#'  \item `P1_transformed` A numeric value representing the first evoked current amplitude (pA) normalized relative to the mean amplitude during the recording's baseline.
#'  \item `P2-transformed` A numeric value representing the second evoked current amplitude (pA) normalized relative to the mean amplitude during the recording's baseline.
#'  \item `amplitude_transformed` (for spontaneous currents only) A numeric value representing the spontaneous current amplitude (pA) normalized relative to the mean amplitude during the recording's baseline.
#' }
#' 
#' 
#' If save_RDS_files is set to "yes" somewhere in the document, the function will also save an .RDS file.
#' @examples
#' make_normalized_EPSC_data(
#'   filename = "Data/Sample-eEPSC-data.csv",
#'   current_type = "eEPSC",
#'   min_time_value = 0,
#'   max_time_value = 25,
#'   interval_length = 5,
#'   baseline_length = 5,
#'   negative_transform_currents = "yes"
#' )
#' @export


make_normalized_EPSC_data <- function(filename = "Data/Sample-eEPSC-data.csv",
                                      current_type = "eEPSC",
                                      min_time_value = 0,
                                      max_time_value = 25,
                                      baseline_length = 5,
                                      interval_length = 5,
                                      negative_transform_currents = "yes") {
  
  list_of_argument_names <- c(filename, current_type)
  
  if (is.null(current_type) ||
      length(current_type) != 1L ||
      !current_type %in% c("eEPSC", "sEPSC")) {
    stop("'current_type' argument must be one of: 'eEPSC' or 'sEPSC'")
  }
  
  if (!negative_transform_currents %in% c("yes", "no")) {
    stop("'negative_transform_currents' argument must be one of: 'yes' or 'no'")
  }
  
  
  if (current_type == "eEPSC") {
    if (any(grepl("sEPSC", list_of_argument_names))) {
      stop(
        "current_type = \"",
        current_type,
        "\" but some arguments have the text ",
        "\"sEPSC\".",
        "\n Are you sure that you selected the correct current type?"
      )
    }
  }
  
  if (current_type == "sEPSC") {
    if (any(grepl("eEPSC", list_of_argument_names))) {
      stop(
        "current_type = \"",
        current_type,
        "\" but some arguments have the text ",
        "\"eEPSC\".",
        "\n Are you sure that you selected the correct current type?"
      )
    }
  }
  
  if (max_time_value %% baseline_length != 0) {
    stop("max_time_value is ", max_time_value, ", which is not divisible by interval_length, ", interval_length)
  }
  
  raw_df <- read.csv(here(filename), header = TRUE) %>%
    rename_with(tolower)
  
  if (current_type == "eEPSC") {
    raw_df <- raw_df %>%
      rename(
        ID = id,
        P1 = p1,
        P2 = p2,
        X = x,
        Y = y
      )
  }
  
  if (current_type == "sEPSC") {
    raw_df <- raw_df %>%
      rename(ID = id, X = x, Y = y)
  }
  
  raw_df <- raw_df %>%
    mutate(across(c(
      ID, letter, category, treatment, sex, synapses
    ), as.factor)) %>%
    filter(time <= max_time_value)
  
  
  if (current_type == "eEPSC") {
    if (negative_transform_currents == "yes") {
    raw_df <- raw_df %>%
      mutate(P1 = P1 * -1,
             # Need positive current amplitude values to make plots more intuitive
             P2 = P2 * -1,
             PPR = P2 / P1)
    } else {
      raw_df <- raw_df %>%
        mutate(PPR = P2 / P1)
    }
  }
  
  # Divide data into intervals (e.g. 5-min intervals)
  
  
  time_sequence <- seq(from = min_time_value, to = max_time_value, by = interval_length)
  time_labels <- head(paste0("t", time_sequence, "to", time_sequence + interval_length),
                      -1)
  
  raw_df <- raw_df %>%
    mutate(interval = cut(
      time,
      breaks = time_sequence,
      include.lowest = TRUE,
      labels = time_labels
    )) %>%
    group_by(letter)
  
  # Within each cell, normalize all of the eEPSC amplitudes
  # relative to the mean baseline amplitude
  
  if (current_type == "eEPSC") {
    raw_df <- raw_df %>%
      mutate(
        baseline_range = (time <= baseline_length),
        baseline_mean = sum(P1 * baseline_range) / sum(baseline_range),
        P1_transformed = (P1 / baseline_mean) * 100,
        P2_transformed = (P2 / baseline_mean) * 100
      )
  }
  
  if (current_type == "sEPSC") {
    raw_df <- raw_df %>%
      mutate(
        baseline_range = (time <= baseline_length),
        baseline_mean = sum(amplitude * baseline_range) / sum(baseline_range),
        amplitude_transformed = (amplitude / baseline_mean) * 100
      )
  }
  
  assign(paste0("raw_", current_type, "_df"), raw_df, envir = .GlobalEnv)
  
  if (save_RDS_files == "yes") {
    saveRDS(raw_df, file = here(
      paste0("Data/Output-Data-from-R/raw_", current_type, "_df.RDS")
    ))
  }
}

#' Prune and summarize current data per minute
#' 
#' make_pruned_eEPSC_df() creates a dataframe of evoked or spontaneous current data summarized to a user-defined interval length.
#' This is equivalent to GraphPad Prism's "prune rows" function to reduce data to summary values for every n rows.
#' @param data A `data.frame` object generated using [make_normalized_EPSC_data()]. It must contain the following columns:
#' 
#' \itemize{
#'  \item `letter` A character value that is a unique identifier for a single recording. Used to link data sets for evoked or spontaneous currents and cell-characteristics.
#'  \item `synapses` A character value (e.g. "Glutamate" or "GABA").
#'  \item `sex` A character value (e.g. "Male" or "Female").
#'  \item `treatment` A character value (e.g. "Control", "HNMPA").
#'  \item `time` A numeric value that represents time in minutes. This column is autogenerated in [add_new_cells()].
#'  \item `ID` A character value for the recording filename.
#'  \item `P1` A numeric value representing the amplitude of the first evoked current in pA.
#'  \item `P2` A numeric value representing the amplitude of the second evoked current in pA.
#'  \item `PPR` A numeric value that represents the paired pulse ratio (PPR) of the evoked currents, generated using `mutate(PPR = P2/P1)` in [make_normalized_EPSC_data()]
#'  \item `X` A numeric value representing the x-value of the cell's location in µm.
#'  \item `Y` A numeric value representing the y-value of the cell's location in µm.
#'  \item `age` A numeric value representing the animal's age. Can be any value as long as the time units are consistent throughout (e.g. don't mix up days and months when reporting animal ages).
#'  \item `animal` A numeric value representing the animal's ID or number.
#'  \item `category` A numeric value representing the experiment type. Used to assign top-level groups for further analyses, with `treatment` as subgroups.
#'  \item `cell` A character or numeric value representing the cell. For example, use `3.1.1` for animal #3, slice #1, cell #1.
#'  \item `interval` A character value indicating the interval that the data belong to (e.g. "t0to5" for the first 5 minutes, "t5to10"). Generated automatically in [make_normalized_EPSC_data()].
#'  \item `baseline_mean` A numeric value representing the mean evoked current amplitude during the baseline period. There is a different baseline_mean for each letter.
#'  \item `P1_transformed` A numeric value representing the first evoked current amplitude (pA) normalized relative to the mean amplitude during the recording's baseline.
#'  \item `P2-transformed` A numeric value representing the second evoked current amplitude (pA) normalized relative to the mean amplitude during the recording's baseline.
#'  \item `amplitude_transformed` (for spontaneous currents only) A numeric value representing the spontaneous current amplitude (pA) normalized relative to the mean amplitude during the recording's baseline.
#' }
#' 
#' 
#' @param current_type A character describing the current type. Allowed values are "eEPSC" or "sEPSC".
#' @param min_time_value Minimum time value (in minutes), which defaults to 0.
#' @param max_time_value Maximum recording length (in minutes). All data points will be filtered to time values less than or equal to this value. Defaults to 25.
#' @param baseline_length Length of the baseline (in minutes). Refers to data collected before applying a hormone, antagonist, or a protocol like high frequency stimulation. Defaults to 5.
#' @param interval_length Length of each interval (in minutes). Used to divide the dataset into broad ranges for statistical analysis. Defaults to 1 for one summary point per minute.
#' @returns Three dataframes that can be viewed and used for further analyses in R. These are:
#' 
#' \itemize{
#'  \item `pruned_eEPSC_df_individual_cells`/`pruned_sEPSC_df_individual_cells` A dataframe containing current data for each individual cell, but the data are reduced to a summary  point per per minute (or another value if a different `interval_length` is set).
#'  New columns include mean amplitude (`mean_P1` in pA), standard deviation (`sd_P1`), standard error (`se`), coefficient of variation (`cv`) and, inverse coefficient of variation squared (`cv_inverse_square`).
#'  \item `pruned_eEPSC_df_all_cells`/`pruned_sEPSC_df_all_cells` A dataframe consisting of data grouped per __
#'  \item `pruned_eEPSC_df_for_table`/`pruned_sEPSC_df_for_table` A dataframe containing two columns: letter and `P1_transformed` (for `eEPSC`) or `spont_amplitude_transformed` (for `sEPSC`). The current data (`P1_transformed``spont_amplitude_transformed`) is collapsed into a single row for each letter,
#'  with the current data for each letter stored as a list. This is required to create sparklines of current amplitude over time within the cell summary table. See [make_cell_summary_df()] and [make_interactive_summary_table()].
#'  }
#' 
#' @export
#' @examples
#' make_pruned_EPSC_data(data = raw_eEPSC_df,
#'   current_type = "eEPSC",
#'   min_time_value = 0,
#'   max_time_value = 25,
#'   baseline_length = 5,
#'   interval_length = 1)

make_pruned_EPSC_data <- function(data = raw_eEPSC_df,
                                  current_type = "eEPSC",
                                  min_time_value = 0,
                                  max_time_value = 25,
                                  baseline_length = 5,
                                  interval_length = 1) {
  time_sequence <- seq(from = min_time_value, to = max_time_value, by = interval_length)
  time_labels <- head(paste0("t", time_sequence, "to", time_sequence + interval_length),-1)
  
  if (is.null(current_type) ||
      length(current_type) != 1L ||
      !current_type %in% c("eEPSC", "sEPSC")) {
    stop("'current_type' argument must be one of: 'eEPSC' or 'sEPSC'")
  }
  
  # Prune within individual cells
  pruned_df_individual_cells <- data %>%
    mutate(
      interval_pruned = cut(
        time,
        breaks = time_sequence,
        include.lowest = TRUE,
        labels = time_labels
      )
    ) %>%
    group_by(category, letter, sex, treatment, interval_pruned)
  
  if (current_type == "eEPSC") {
    pruned_df_individual_cells <- pruned_df_individual_cells %>%
      reframe(
        mean_P1 = mean(P1, na.rm = TRUE),
        # Mean amplitude per minute across all cells
        sd_P1 = sd(P1, na.rm = TRUE),
        n = n(),
        se = sd_P1 / sqrt(n),
        cv = sd_P1 / mean_P1,
        cv_inverse_square = 1 / (cv ^ 2),
        letter = unique(letter),
        category = unique(category),
        time = last(time),
        baseline_mean = unique(baseline_mean),
        synapses = unique(synapses)
      )
    
    pruned_df_for_table <- pruned_df_individual_cells %>%
      group_by(letter) %>% 
      summarize(P1_transformed = list(mean_P1))
  }
  
  if (current_type == "sEPSC") {
    pruned_df_individual_cells <- pruned_df_individual_cells %>%
      reframe(
        mean_amplitude = mean(amplitude_transformed, na.rm = TRUE),
        mean_raw_amplitude = mean(amplitude, na.rm = TRUE),
        sd_amplitude = sd(amplitude_transformed, na.rm = TRUE),
        n = n(),
        # Gets number of currents within each minute
        frequency = n / 60,
        # Frequency in Hz
        se = sd_amplitude / sqrt(n),
        letter = unique(letter),
        category = unique(category),
        interval = unique(interval),
        synapses = unique(synapses),
        time = last(time) # Time value at the end of the interval; used for plots
      ) %>%
      group_by(letter) %>%
      # Obtain normalized frequency
      mutate(
        baseline_range = (time <= baseline_length),
        baseline_mean_frequency = sum(frequency * baseline_range) / sum(baseline_range),
        frequency_transformed = (frequency / baseline_mean_frequency) * 100
      )
    
    pruned_df_for_table <- pruned_df_individual_cells %>% 
      group_by(letter) %>% 
      summarize(spont_amplitude_transformed = list(mean_amplitude))
  }
  
  assign(
    paste0("pruned_", current_type, "_df_individual_cells"),
    pruned_df_individual_cells,
    envir = .GlobalEnv
  )
  
  assign(
    paste0("pruned_", current_type, "_df_for_table"),
    pruned_df_for_table,
    envir = .GlobalEnv
  )
  
  if (save_RDS_files == "yes") {
    saveRDS(pruned_df_individual_cells, file = here(
      paste0(
        "Data/Output-Data-from-R/pruned_",
        current_type,
        "_df_individual_cells.RDS"
      )
    ))
  }
  
  # Prune all cells
  if (current_type == "eEPSC") {
    pruned_df_all_cells <- data %>%
      mutate(
        interval_pruned = cut(
          time,
          breaks = time_sequence,
          include.lowest = TRUE,
          labels = time_labels
        )
      ) %>%
      group_by(category, letter, sex, treatment, interval_pruned) %>%
      reframe(
        mean_P1 = mean(P1_transformed, na.rm = TRUE),
        sd_P1 = sd(P1_transformed, na.rm = TRUE),
        n = n(),
        se = sd_P1 / sqrt(n),
        cv = sd_P1 / mean_P1 * 100,
        letter = unique(letter),
        category = unique(category),
        time = last(time)
      ) %>%
      group_by(category, interval_pruned, sex, treatment) %>%
      reframe(
        mean_P1_all_cells = mean(mean_P1, na.rm = TRUE),
        sd_P1_all_cells = sd(mean_P1, na.rm = TRUE),
        n = n(),
        se_P1_all_cells = sd_P1_all_cells / sqrt(n),
        cv_P1_all_cells = sd_P1_all_cells / mean_P1_all_cells * 100,
        time = last(time),
        category = unique(category),
        
      )
  }
  
  if (current_type == "sEPSC") {
    pruned_df_all_cells <- pruned_df_individual_cells %>%
      ungroup() %>%
      group_by(category, interval_pruned, sex, treatment) %>%
      reframe(
        mean_all_amplitude = mean(mean_amplitude, na.rm = TRUE),
        mean_all_raw_amplitude = mean(mean_raw_amplitude, na.rm = TRUE),
        sd_all_amplitude = sd(mean_amplitude, na.rm = TRUE),
        n = n(),
        se_amplitude = sd_all_amplitude / sqrt(n),
        sd_all_raw_amplitude = sd(mean_raw_amplitude, na.rm = TRUE),
        se_raw_amplitude = sd_all_raw_amplitude / sqrt(n),
        mean_all_frequency = mean(frequency_transformed, na.rm = TRUE),
        sd_all_frequency = sd(frequency_transformed, na.rm = TRUE),
        se_frequency = sd_all_frequency / sqrt(n),
        mean_all_raw_frequency = mean(frequency, na.rm = TRUE),
        sd_all_raw_frequency = sd(frequency, na.rm = TRUE),
        se_raw_frequency = sd_all_raw_frequency / sqrt(n),
        time = last(time),
        interval = unique(interval),
        category = unique(category)
      )
  }
  assign(paste0("pruned_", current_type, "_df_all_cells"),
         pruned_df_all_cells,
         envir = .GlobalEnv)
  
  if (save_RDS_files == "yes") {
    saveRDS(pruned_df_all_cells, file = here(
      paste0(
        "Data/Output-Data-from-R/pruned_",
        current_type,
        "_df_all_cells.RDS"
      )
    ))
  }
}


# TODO : expand description for make_summary data -----
#' Divide and summarize current data into large intervals
#' 
#' `make_summary_EPSC_data()` allows you to divide data from a long recording (e.g. 30 minutes) into evenly-spaced intervals (e.g. 5 minutes).
#' This can be useful for future statistical testing to compare effect sizes across time.
#' 
#' @param data A dataset
#' @param current_type A character describing the current type. Allowed values are "eEPSC" or "sEPSC".
#' 
#' @returns A dataframe with the following columns
#' 
#' @export
#' @examples
#' 
#' make_summary_EPSC_data(data = raw_eEPSC_df, current_type = "eEPSC")
#' 
# Create summary data for each five minute interval
# Useful for future statistical testing
make_summary_EPSC_data <- function(data, current_type) {
  if (is.null(current_type) ||
      length(current_type) != 1L ||
      !current_type %in% c("eEPSC", "sEPSC")) {
    stop("'current_type' argument must be one of: 'eEPSC' or 'sEPSC'")
  }
  
  if (current_type == "eEPSC") {
    summary_df <- data %>%
      group_by(category, letter, sex, treatment, interval) %>%
      summarize(
        mean_P1_transformed = mean(P1_transformed, na.rm = TRUE),
        mean_P1_raw = mean(P1, na.rm = TRUE),
        n = n(),
        sd = sd(P1_transformed, na.rm = TRUE),
        cv = sd / mean_P1_transformed,
        se = sd(P1_transformed, na.rm = TRUE) / sqrt(n),
        cv_inverse_square = 1 / (cv ^ 2),
        variance = var(P1_transformed, na.rm = TRUE),
        VMR = variance / mean_P1_transformed,
        age = unique(age),
        # unique() retains unmodified columns that still should be included
        animal = unique(animal),
        X = unique(X),
        Y = unique(Y),
        time = last(time),
        synapses = unique(synapses)
      ) %>%
      ungroup()
  }
  
  
  if (current_type == "sEPSC") {
    summary_df <- data %>%
      group_by(category, letter, sex, treatment, interval) %>%
      summarize(
        mean_transformed_amplitude = mean(mean_amplitude, na.rm = TRUE),
        mean_raw_amplitude = mean(mean_raw_amplitude, na.rm = TRUE),
        sd_all_amplitude = sd(mean_amplitude, na.rm = TRUE),
        n = n(),
        se_amplitude = sd_all_amplitude / sqrt(n),
        mean_transformed_frequency = mean(frequency_transformed, na.rm = TRUE),
        sd_transformed_frequency = sd(frequency_transformed, na.rm = TRUE),
        se_frequency = sd_transformed_frequency / sqrt(n),
        mean_raw_frequency = mean(frequency, na.rm = TRUE),
        time = last(time),
        interval = unique(interval),
        category = unique(category),
        synapses = last(synapses)
      )
  }
  assign(paste0("summary_", current_type, "_df"), summary_df, envir = .GlobalEnv)
}


# Summary table ----

make_cell_summary_df <- function(cell_characteristics_df,
                                 include_all_treatments = "yes",
                                 list_of_treatments = NULL,
                                 include_all_categories = "yes",
                                 list_of_categories = NULL) {
  cell_characteristics_for_table <- cell_characteristics_df %>%
    mutate(
      R_a = lapply(str_split(R_a, pattern = ", "), FUN = as.numeric),
      R_a = lapply(R_a, FUN = replace_na, replace = 0),
      letter = factor(letter)
    )
  
  table_data <-
    merge(pruned_eEPSC_df_for_table, pruned_sEPSC_df_for_table, by = "letter") %>%
    merge(., cell_characteristics_for_table, by = "letter") %>%
    merge(., treatment_names_and_colours, by = "treatment") %>%
    select(
      c(
        letter,
        display_names,
        treatment,
        synapses,
        sex,
        P1_transformed,
        spont_amplitude_transformed,
        R_a,
        X,
        Y,
        age,
        animal,
        category,
        cell,
        colours
      )
    ) %>%
    rename_with(str_to_title) %>%
    mutate(X = round(X, -1),
           Y = round(Y, -1))
  
  if (include_all_treatments == "yes") {
    if (!is.null(list_of_treatments)) {
      warning(
        "include_all_treatments = \"yes\", but you included a list of treatments to filter. All treatments will be used."
      )
    }
    
  } else {
    if (is.null(list_of_treatments)) {
      stop(
        "include_all_treatments = \"",
        include_all_treatments,
        "\", but list_of_treatments is NULL.",
        "\nDid you forget to add a list of treatments?"
      )
    }
    
    if (!is.character(list_of_treatments)) {
      stop(
        "include_all_treatments = \"",
        include_all_treatments,
        "\", but list_of_treatments is not a character object.",
        "\nDid you forget to add a list of treatments?"
      )
    }
    
    table_data <- table_data %>%
      filter(Treatment %in% list_of_treatments)
  }
  # Category filter
  
  if (include_all_categories == "yes") {
    if (!is.null(list_of_categories)) {
      warning(
        "include_all_categories = \"yes\", but you included a list of categories to filter. All categories will be used."
      )
    }
    
  } else {
    if (is.null(list_of_categories)) {
      stop(
        "include_all_categories = \"",
        include_all_categories,
        "\", but list_of_categories is NULL.",
        "\nDid you forget to add a list of categories?"
      )
    }
    
    if (!is.character(list_of_categories)) {
      stop(
        "include_all_categories = \"",
        include_all_categories,
        "\", but list_of_categories is not a character object.",
        "\nDid you forget to add a list of categories?"
      )
    }
    
    table_data <- table_data %>%
      filter(Category %in% list_of_categories)
  }
  
  assign("summary_table_data", table_data, envir = .GlobalEnv)
}


make_interactive_summary_table <- function(dataframe) {
  cell_table <- reactable(
    data = dataframe,
    defaultSorted = c("Category", "Treatment", "Animal"),
    filterable = TRUE,
    showPageSizeOptions = TRUE,
    elementId = "cell-table",
    defaultPageSize = 15,
    defaultColDef = colDef(vAlign = "center", headerVAlign = "center"),
    columns = list(
      Colours = colDef(show = FALSE),
      Letter = colDef(
        name = "Letter",
        sticky = "left",
        style = list(borderRight = "1px solid #eee"),
        headerStyle = list(borderRight = "1px solid #eee")
      ),
      Display_names = colDef(name = "Treatment"),
      Treatment = colDef(show = FALSE),
      Sex = colDef(
        name = "Sex",
        filterMethod = JS(
          "function(rows, columnId, filterValue) {
        return rows.filter(function(row) {
          return String(row.values[columnId]).toUpperCase() === filterValue.toUpperCase()
        })
      }"
      )
      ),
      P1_transformed = colDef(
        name = "eEPSC amplitude (pA)",
        filterable = FALSE,
        cell = react_sparkline(
          dataframe,
          line_color_ref = "Colours",
          show_area = TRUE,
          area_opacity = 1
        )
      ),
      Spont_amplitude_transformed = colDef(
        name = "sEPSC amplitude (pA)",
        filterable = FALSE,
        cell = react_sparkline(
          dataframe,
          line_color_ref = "Colours",
          show_area = TRUE,
          area_opacity = 1
        )
      ),
      R_a = colDef(
        name = "Ra (MΩ)",
        filterable = FALSE,
        cell = react_sparkline(
          dataframe,
          line_color_ref = "Colours",
          labels = c("first", "last"),
          decimals = 1
        )
      )
    )
  )
  cell_table
}



# Plot data ----
make_baseline_comparison_plot <- function(data,
                                          include_all_treatments = "yes",
                                          list_of_treatments = NULL,
                                          baseline_interval,
                                          filename_suffix = "",
                                          current_type,
                                          parameter,
                                          large_axis_text = "no",
                                          plot_width) {
  if (is.null(current_type) ||
      length(current_type) != 1L ||
      !current_type %in% c("eEPSC", "sEPSC")) {
    stop("'current_type' argument must be one of: 'eEPSC' or 'sEPSC'")
  }
  
  if (is.null(baseline_interval) ||
      !is.character(baseline_interval)) {
    stop("'baseline_interval' must be a character (e.g. \"t0to5\", \"t0to3\")")
  }
  
  if (include_all_treatments == "yes") {
    treatment_info <- treatment_names_and_colours
    plot_data <- data %>%
      filter(treatment %in% treatment_names_and_colours$treatment) %>%
      droplevels()
    
    if (!is.null(list_of_treatments)) {
      warning(
        "include_all_treatments = \"yes\", but you included a list of treatments to filter. All treatments will be used."
      )
    }
    
  } else {
    if (is.null(list_of_treatments)) {
      stop(
        "include_all_treatments = \"",
        include_all_treatments,
        "\", but list_of_treatments is NULL.",
        "\nDid you forget to add a list of treatments?"
      )
    }
    
    if (!is.character(list_of_treatments)) {
      stop(
        "include_all_treatments = \"",
        include_all_treatments,
        "\", but list_of_treatments is not a character object.",
        "\nDid you forget to add a list of treatments?"
      )
    }
    
    treatment_info <- treatment_names_and_colours %>%
      filter(treatment %in% list_of_treatments)
    plot_data <- data %>%
      filter(treatment %in% list_of_treatments) %>%
      droplevels()
  }
  
  if (current_type == "sEPSC") {
    filepath <- "Figures/Spontaneous-currents/Output-summary-plots"
    
    allowed_parameters_list <- "\"raw_amplitude\", or \"raw_frequency\""
    
    if (!parameter %in% c("raw_amplitude", "raw_frequency")) {
      stop(
        "parameter must be ",
        allowed_parameters_list,
        " for current_type \"",
        current_type,
        "\" \nbecause transformed data are all 100 during the baseline."
      )
    }
    
    if (parameter == "raw_amplitude") {
      y_var <- "mean_raw_amplitude"
      y_title <- "Baseline sEPSC Amplitude (pA)"
    }
    
    if (parameter == "raw_frequency") {
      y_var <- "mean_raw_frequency"
      y_title <- "Baseline sEPSC Frequency (Hz)"
    }
  }
  
  if (current_type == "eEPSC") {
    allowed_parameters_list <- c("\"raw_amplitude\"")
    
    if (!parameter %in% c("raw_amplitude")) {
      stop(
        "parameter must be ",
        allowed_parameters_list,
        " for current_type \"",
        current_type,
        "\"",
        " \nbecause transformed data are all 100 during the baseline."
      )
    }
    
    filepath <- "Figures/Evoked-currents/Output-summary-plots"
    
    if (parameter == "raw_amplitude") {
      y_var <- "mean_P1_raw"
      y_title = "Baseline eEPSC Amplitude (pA)"
    }
  }
  
  baseline_comparison_plot <- plot_data %>%
    filter(interval == baseline_interval) %>%
    mutate(
      treatment = str_replace_all(
        treatment,
        setNames(treatment_info$display_names, treatment_info$treatment)
      ),
      treatment = factor(treatment, levels = treatment_info$display_names)
    ) %>%
    ggplot(aes(
      x = treatment,
      y = .data[[y_var]],
      color = treatment,
      shape = sex
    )) +
    labs(x = NULL, y = y_title) +
    geom_sina(
      bw = 12,
      alpha = 0.8,
      maxwidth = 0.5,
      size = 2
    ) +
    stat_summary(
      fun.data = mean_se,
      geom = "pointrange",
      color = mean_point_colour,
      size = mean_point_size + 0.2,
      alpha = 0.8
    ) +
    scale_color_manual(breaks = treatment_info$display_names,
                       values = treatment_info$colours) +
    scale_shape_manual(values = c(17, 16)) +
    theme(legend.position = "right",
          legend.background = element_rect(fill = NA),
    ) +
    guides(color = "none", shape = guide_legend(reverse = TRUE))
  
  if (large_axis_text == "yes") {
    baseline_comparison_plot <- baseline_comparison_plot +
      theme(
        axis.text.x = element_text(size = 24, margin = margin(t = 10)),
        axis.title.y = element_text(size = 24, face = "plain")
      )
  }
  
  if (save_plot_pngs == "yes") {
    ggsave(
      baseline_comparison_plot,
      path = here(filepath),
      file = paste0(
        "Baseline-",
        parameter,
        "-comparison",
        filename_suffix,
        ".png"
      ),
      width = plot_width,
      height = 5,
      units = "in",
      dpi = 300
    )
  }
  baseline_comparison_plot
}

make_raw_plots <-
  function(data,
           plot_treatment,
           plot_category,
           current_type,
           parameter,
           pruned,
           hormone_added,
           hormone_or_HFS_start_time) {
    list_of_plots <- list()
    
    if (is.null(current_type) ||
        length(current_type) != 1L ||
        !current_type %in% c("eEPSC", "sEPSC")) {
      stop("'current_type' argument must be one of: 'eEPSC' or 'sEPSC'")
    }
    
    if (is.null(hormone_added) ||
        length(hormone_added) != 1L ||
        !is.character(hormone_added)) {
      stop("\"hormone_added\" must be a character. Use \"HFS\" for high-frequency stimulation or \"Insulin\", \"CCK\", etc. for any hormones")
    }
    
    
    if (is.null(hormone_or_HFS_start_time) ||
        !is.numeric(hormone_or_HFS_start_time)) {
      stop(
        "\"hormone_or_HFS_start_time\" must be numeric (e.g. 5 for HFS or a hormone applied at five minutes into the recording)."
      )
    }
    
    df <- data %>%
      filter(category == plot_category) %>%
      filter(treatment == plot_treatment)
    
    
    letters <- as.character(unique(unlist(df$letter)))
    
    plot_colour <- treatment_names_and_colours %>%
      filter(treatment == plot_treatment) %>%
      pull(colours)
    
    treatment_label <- treatment_names_and_colours %>%
      filter(treatment == plot_treatment) %>%
      pull(treatment)
    
    if (current_type == "eEPSC") {
      # The plots should go to specific folders depending on current type
      filepath <- "Figures/Evoked-currents/Output-individual-plots"
      
      allowed_parameters_list <- "\"P1\", \"mean_P1\", or \"PPR\""
      
      if (!parameter %in% c("P1", "mean_P1", "PPR")) {
        stop(
          "parameter must be ",
          allowed_parameters_list,
          " for current_type \"",
          current_type,
          "\". \nCheck parameter, current_type or data."
        )
      }
      
      if (parameter == "P1") {
        y_title <- "eEPSC Amplitude (pA)"
        
        if (pruned == "yes") {
          stop("pruned = \"yes\"), but parameter = \"P1\". ",
               "\nDid you mean \"mean_P1\"?. ")
        }
      }
      
      if (parameter == "mean_P1") {
        y_title <- "eEPSC Amplitude (pA)"
        
        
        if (pruned == "no") {
          stop(
            "parameter = \"mean_P1\", but pruned = \"no\". ",
            "Are you trying to create a pruned plot? \nIf so, change pruned to \"yes\" ",
            "and ensure that you have specified the correct dataframe in the data argument"
          )
        }
      }
      
      if (parameter == "PPR") {
        y_title <- "Paired-pulse Ratio"
      }
    }
    
    
    if (current_type == "sEPSC") {
      filepath <- "Figures/Spontaneous-currents/Output-individual-plots"
      
      allowed_parameters_list <- "\"amplitude\" or \"frequency\""
      
      if (!parameter %in% c("amplitude", "frequency")) {
        stop(
          "parameter must be ",
          allowed_parameters_list,
          " for current_type \"",
          current_type,
          "\". \nCheck parameter, current_type or data."
        )
      }
      
      
      if (parameter == "amplitude") {
        y_title <- "sEPSC Amplitude (pA)"
      }
      
      if (parameter == "frequency") {
        y_title <- "sEPSC Frequency (Hz)"
      }
    }
    
    if (pruned == "yes") {
      annotation <- "_pruned"
    } else {
      annotation <- ""
    }
    
    # Pruned sEPSC amplitude plots use mean +/- SE, unlike the other plots
    
    
    for (i in letters) {
      plot_df <- df %>% filter(letter == i)
      if (current_type == "sEPSC" &
          pruned == "yes" &
          parameter == "amplitude") {
        y_title <- "sEPSC Amplitude (pA)"
        
        list_of_plots[[i]] <- ggplot(
          plot_df,
          aes(
            x = time,
            y = mean_amplitude,
            ymin = mean_amplitude - se,
            ymax = mean_amplitude + se
          )
        )
        
      } else {
        list_of_plots[[i]] <- ggplot(plot_df, aes(x = time, y = .data[[parameter]]))
        
      }
      
      list_of_plots[[i]] <- list_of_plots[[i]] +
        ggtitle(
          paste("Recording", i),
          subtitle = paste(
            "Treatment:",
            treatment_label,
            " Sex:",
            unique(plot_df$sex)
          )
        ) +
        labs(x = "Time (min)", y = y_title)
      
      if (current_type == "sEPSC" &
          pruned == "yes" &
          parameter == "amplitude") {
        list_of_plots[[i]] <- list_of_plots[[i]] +
          geom_pointrange(
            shape = if (unique(plot_df$sex) == "Male") {
              male_shape
            } else {
              female_shape
            },
            colour = plot_colour,
            size = 1,
            alpha = 1
          )
      } else {
        list_of_plots[[i]] <- list_of_plots[[i]] +
          geom_point(
            shape = if (unique(plot_df$sex) == "Male") {
              male_shape
            } else {
              female_shape
            },
            colour = plot_colour,
            size = if (current_type == "sEPSC" & pruned == "no") {
              1
            } else {
              3.5
            },
            alpha = if (pruned == "yes") {
              1
            } else {
              0.7
            }
          )
      }
      
      # Get limits of x- and y-axes
      ymax <- layer_scales(list_of_plots[[i]])$y$get_limits()[2]
      xmax <- layer_scales(list_of_plots[[i]])$x$get_limits()[2]
      ymax2 <- layer_scales(list_of_plots[[i]])$y$get_limits()[2]
      
      # If hormone_added = Insulin, CCK, i.e. anything other than "HFS" (high frequency stimulation),
      # add an annotated line over the application period:
      
      if (hormone_added != "HFS") {
        list_of_plots[[i]] <-
          list_of_plots[[i]] +
          annotate(
            geom = "segment",
            x = hormone_or_HFS_start_time,
            xend = xmax,
            y = ymax + 0.1 * ymax,
            yend = ymax + 0.1 * ymax,
            colour = line_col,
            linewidth = 0.6
          )
        
        list_of_plots[[i]] <-
          list_of_plots[[i]] +
          annotate(
            geom = "text",
            x = hormone_or_HFS_start_time,
            y = ymax2 + 0.16 * ymax2,
            label = hormone_added,
            size = 4,
            hjust = 0,
            family = plot_font_family
          )
      }
      
      # If hormone_added == HFS (experiments involving HFS)
      # add a labelled arrow at 5 minutes:
      
      if (hormone_added == "HFS") {
        list_of_plots[[i]] <-
          list_of_plots[[i]] +
          annotate(
            geom = "segment",
            x = hormone_or_HFS_start_time,
            y = ymax + 0.22 * ymax,
            xend = hormone_or_HFS_start_time,
            yend = ymax + 0.10 * ymax,
            arrow = arrow(type = "closed", length = unit(0.02, "npc"))
          )
        
        
        list_of_plots[[i]] <-
          list_of_plots[[i]] +
          annotate(
            geom = "text",
            x = hormone_or_HFS_start_time,
            y = ymax + 0.27 * ymax,
            label = "HFS",
            size = 3.5,
            hjust = 0.5,
            family = plot_font_family
          )
      }
      
      if (save_plot_pngs == "yes") {
        ggsave(
          list_of_plots[[i]],
          path = here(filepath),
          file = paste0(i, annotation, ".png"),
          width = 7,
          height = 5,
          units = "in",
          dpi = 300
        )
      }
    }
    
    # Assign plot_list to a global variable to access plots outside of this function
    
    assign(
      paste0(
        current_type,
        "_",
        parameter,
        annotation,
        "_plots_",
        plot_treatment,
        "_category_",
        plot_category
      ),
      list_of_plots,
      envir = .GlobalEnv
    )
  }


perform_t_tests_for_summary_plot <- function(data,
                                             test_category,
                                             include_all_treatments = "yes",
                                             list_of_treatments = NULL,
                                             baseline_interval,
                                             interval_length,
                                             parameter,
                                             current_type) {
  if (is.null(current_type) ||
      length(current_type) != 1L ||
      !current_type %in% c("eEPSC", "sEPSC")) {
    stop("'current_type' argument must be one of: 'eEPSC' or 'sEPSC'")
  }
  
  if (include_all_treatments == "yes") {
    treatment_info <- treatment_names_and_colours
    t_test_data <- data %>%
      semi_join(treatment_info, by = c("treatment")) %>%
      filter(treatment %in% treatment_info$treatment)
    
    if (!is.null(list_of_treatments)) {
      warning(
        "include_all_treatments = \"yes\", but you included a list of treatments to filter. All treatments will be used."
      )
    }
    
  } else {
    if (is.null(list_of_treatments)) {
      stop(
        "include_all_treatments = \"",
        include_all_treatments,
        "\", but list_of_treatments is NULL.",
        "\nDid you forget to add a list of treatments?"
      )
    }
    
    if (!is.character(list_of_treatments)) {
      stop(
        "include_all_treatments = \"",
        include_all_treatments,
        "\", but list_of_treatments is not a character object.",
        "\nDid you forget to add a list of treatments?"
      )
    }
    
    if (is.null(baseline_interval) ||
        !is.character(baseline_interval)) {
      stop("'baseline_interval' must be a character (e.g. \"t0to5\" or \"t0to3\")")
    }
    
    treatment_info <- treatment_names_and_colours %>%
      filter(treatment %in% list_of_treatments)
    t_test_data <- data %>%
      filter(treatment %in% list_of_treatments) %>%
      droplevels()
  }
  
  if (current_type == "eEPSC") {
    allowed_parameters_list <- "\"amplitude\""
    
    if (!parameter %in% c("amplitude")) {
      stop(
        "parameter must be ",
        allowed_parameters_list,
        " for current_type \"",
        current_type,
        "\". \nCheck parameter, current_type or data."
      )
    }
    
    if (parameter == "amplitude") {
      t_test_results <- t_test_data %>%
        filter(category == test_category) %>%
        group_by(treatment) %>%
        pairwise_t_test(
          mean_P1_transformed ~ interval,
          ref.group = baseline_interval,
          paired = TRUE,
          p.adjust.method = "holm"
        )
    }
  }
  
  if (current_type == "sEPSC") {
    allowed_parameters_list <- "\"amplitude\", \"raw_amplitude\", \"raw_frequency\", or \"frequency\""
    
    if (!parameter %in% c("amplitude",
                          "raw_amplitude",
                          "frequency",
                          "raw_frequency")) {
      stop(
        "parameter must be ",
        allowed_parameters_list,
        " for current_type \"",
        current_type,
        "\". \nCheck parameter, current_type or data."
      )
    }
    
    if (parameter == "amplitude") {
      t_test_results <- t_test_data %>%
        filter(category == test_category) %>%
        group_by(treatment) %>%
        pairwise_t_test(
          mean_transformed_amplitude ~ interval,
          ref.group = baseline_interval,
          paired = TRUE,
          p.adjust.method = "holm"
        )
    }
    
    if (parameter == "raw_amplitude") {
      t_test_results <- t_test_data %>%
        filter(category == test_category) %>%
        group_by(treatment) %>%
        pairwise_t_test(
          mean_raw_amplitude ~ interval,
          ref.group = baseline_interval,
          paired = TRUE,
          p.adjust.method = "holm"
        )
    }
    
    if (parameter == "frequency") {
      t_test_results <- t_test_data %>%
        filter(category == test_category) %>%
        group_by(treatment) %>%
        pairwise_t_test(
          mean_transformed_frequency ~ interval,
          ref.group = baseline_interval,
          paired = TRUE,
          p.adjust.method = "holm"
        )
    }
    
    if (parameter == "raw_frequency") {
      t_test_results <- t_test_data %>%
        filter(category == test_category) %>%
        group_by(treatment) %>%
        pairwise_t_test(
          mean_raw_frequency ~ interval,
          ref.group = baseline_interval,
          paired = TRUE,
          p.adjust.method = "holm"
        )
    }
  }
  
  t_test_table <- t_test_results %>%
    mutate(
      statistic = round(statistic, 2),
      p_string = pvalString(p.adj),
      significance_stars = case_when(p.adj.signif == "ns" ~ "", T ~ p.adj.signif)
    )
  
  # Need sequence of integers from 1 to the maximum number of intervals
  # to generate asterisk_time as a function of the number of intervals
  integer_sequence <- seq(1, length(unique(t_test_table$group2)), by = 1)
  
  positions_df <- data.frame(
    group2 = unique(t_test_table$group2),
    asterisk_time = (interval_length / 2) + interval_length * integer_sequence
  )
  
  t_test_table <- merge(positions_df, t_test_table, by = "group2") %>%
    select(
      treatment,
      .y.,
      group1,
      group2,
      n1,
      n2,
      statistic,
      df,
      p_string,
      significance_stars,
      asterisk_time
    ) %>%
    arrange(match(treatment, treatment_info$treatment))
  
  assign(paste0("t_test_", current_type, "_", parameter),
         t_test_table,
         envir = .GlobalEnv)
  
  if (save_RDS_files == "yes") {
    saveRDS(t_test_table, file = here(
      paste0("Data/Output-Data-from-R/t_test_", current_type, ".RDS")
    ))
  }
  
  t_test_table
}


make_summary_plots <- function(plot_category,
                               plot_treatment,
                               data,
                               current_type,
                               parameter,
                               hormone_added,
                               hormone_or_HFS_start_time,
                               include_representative_trace = "yes",
                               signif_stars = "no",
                               large_axis_text = "no",
                               shade_intervals = "no") {
  if (is.null(current_type) ||
      length(current_type) != 1L ||
      !current_type %in% c("eEPSC", "sEPSC")) {
    stop("'current_type' argument must be one of: 'eEPSC' or 'sEPSC'")
  }
  
  if (is.null(hormone_added) ||
      length(hormone_added) != 1L ||
      !is.character(hormone_added)) {
    stop("\"hormone_added\" must be a character. Use \"HFS\" for high-frequency stimulation or \"Insulin\", \"CCK\", etc. for any hormones")
  }
  
  
  if (is.null(hormone_or_HFS_start_time) ||
      !is.numeric(hormone_or_HFS_start_time)) {
    stop(
      "\"hormone_or_HFS_start_time\" must be numeric (e.g. 5 for HFS or a hormone applied at five minutes into the recording)."
    )
  }
  
  df <-
    data %>% filter(category == plot_category) %>% filter(treatment == plot_treatment)
  
  plot_colour <- treatment_names_and_colours %>%
    filter(treatment == plot_treatment) %>%
    pull(colours)
  
  plot_colour_pale <- treatment_names_and_colours %>%
    filter(treatment == plot_treatment) %>%
    pull(very_pale_colours)
  
  if (current_type == "eEPSC") {
    allowed_parameters_list <- "\"amplitude\""
    
    if (!parameter %in% c("amplitude")) {
      stop(
        "parameter must be ",
        allowed_parameters_list,
        " for current_type \"",
        current_type,
        "\". \nCheck parameter, current_type or data."
      )
    }
    
    if (parameter == "amplitude") {
      y_var <- "mean_P1_all_cells"
      se_var <- "se_P1_all_cells"
      filepath <- "Figures/Evoked-currents/Output-summary-plots"
      file_name_ending <- ""
      
      if (large_axis_text == "yes") {
        y_title <- "eEPSC Amplitude\n(% Baseline)"
      } else {
        y_title <- "eEPSC Amplitude (% Baseline)"
      }
    }
  }
  
  
  if (current_type == "sEPSC") {
    filepath <- "Figures/Spontaneous-currents/Output-summary-plots"
    
    allowed_parameters_list <- "\"amplitude\", \"raw_amplitude\", \"raw_frequency\", or \"frequency\""
    
    if (!parameter %in% c("amplitude",
                          "raw_amplitude",
                          "frequency",
                          "raw_frequency")) {
      stop(
        "parameter must be ",
        allowed_parameters_list,
        " for current_type \"",
        current_type,
        "\". \nCheck parameter, current_type or data."
      )
    }
    
    
    if (parameter == "amplitude") {
      y_var <- "mean_all_amplitude"
      se_var <- "se_amplitude"
      file_name_ending <- paste0("_", parameter)
      
      if (large_axis_text == "yes") {
        y_title <- "sEPSC Amplitude\n(% Baseline)"
      } else {
        y_title <- "sEPSC Amplitude (% Baseline)"
      }
      
    }
    
    if (parameter == "raw_amplitude") {
      y_var <- "mean_all_raw_amplitude"
      se_var <- "se_raw_amplitude"
      file_name_ending <- paste0("_", parameter)
      
      y_title <- "sEPSC Amplitude (pA)"
      
    }
    
    if (parameter == "frequency") {
      y_var <- "mean_all_frequency"
      se_var <- "se_frequency"
      file_name_ending <- paste0("_", parameter)
      
      if (large_axis_text == "yes") {
        y_title <- "sEPSC Frequency\n(% Baseline)"
      } else {
        y_title <- "sEPSC Frequency (% Baseline)"
      }
    }
    
    if (parameter == "raw_frequency") {
      y_var <- "mean_all_raw_frequency"
      se_var <- "se_frequency"
      file_name_ending <- paste0("_", parameter)
      
      y_title <- "sEPSC Frequency (Hz)"
      
    }
  }
  
  treatment_plot <- df %>%
    ggplot(aes(
      x = time,
      y = .data[[y_var]],
      ymin = .data[[y_var]] - .data[[se_var]],
      ymax = .data[[y_var]] + .data[[se_var]]
    ))
  
  if (shade_intervals == "yes") {
    treatment_plot <- treatment_plot +
      geom_rect(aes(
        xmin = 5,
        xmax = 10,
        ymin = -5,
        ymax = y_axis_limit
      ), fill = rectangle_shading_colour) +
      geom_rect(aes(
        xmin = 15,
        xmax = 20,
        ymin = -5,
        ymax = y_axis_limit
      ), fill = rectangle_shading_colour)
  }
  
  treatment_plot <- treatment_plot +
    geom_pointrange(
      aes(color = sex, shape = sex),
      size = if (large_axis_text == "yes") {
        1.3
      } else {
        0.9
      },
      alpha = 1,
      position = position_dodge(width = if (current_type == "eEPSC") {
        0.3
      } else {
        0
      })
    ) +
    geom_hline(yintercept = 100, linetype = "dashed") +
    coord_cartesian(ylim = c(0, y_axis_limit)) +
    labs(
      x = "Time (min)",
      y = y_title,
      shape = "Sex",
      color = "Sex"
    )
  
  if (large_axis_text == "yes") {
    treatment_plot <- treatment_plot +
      theme(
        axis.title = element_text(size = 24, face = "plain"),
        legend.title = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.17, 0.13),
        legend.text = element_text(size = 14),
        legend.key.spacing.y = unit(0.5, "cm"),
        legend.background = element_rect(fill = NA)
      )
  }
  
  # Nested if statements enable correct legend labels even if only one sex is present
  
  if (is.na(df$n[df$sex == "Female"][1])) {
    treatment_plot <- treatment_plot +
      scale_shape_manual(values = c(male_shape), labels = c((paste0(
        "Males, n = ", df$n[df$sex == "Male"][1]
      )))) +
      scale_color_manual(values = c(plot_colour), labels = c((paste0(
        "Males, n = ", df$n[df$sex == "Male"][1]
      ))))
    
  } else if (is.na(df$n[df$sex == "Male"][1])) {
    treatment_plot <- treatment_plot +
      scale_shape_manual(values = c(female_shape), labels = c((paste0(
        "Females, n = ", df$n[df$sex == "Female"][1]
      )))) +
      scale_color_manual(values = c(plot_colour_pale),
                         labels = c((paste0(
                           "Females, n = ", df$n[df$sex == "Female"][1]
                         ))))
  } else {
    treatment_plot <- treatment_plot +
      scale_shape_manual(values = c(female_shape, male_shape),
                         labels = c((paste0(
                           "Females, n = ", df$n[df$sex == "Female"][1]
                         )), (paste0(
                           "Males, n = ", df$n[df$sex == "Male"][1]
                         )))) +
      scale_color_manual(values = c(plot_colour_pale, plot_colour),
                         labels = c((paste0(
                           "Females, n = ", df$n[df$sex == "Female"][1]
                         )), (paste0(
                           "Males, n = ", df$n[df$sex == "Male"][1]
                         ))))
  }
  
  # Get limits of x- and y-axes
  ymax <- y_axis_limit - 25
  xmax <-
    layer_scales(treatment_plot)$x$get_limits()[2]
  
  # If hormone_added = Insulin, CCK, i.e. anything other than "HFS" (high frequency stimulation),
  # add an annotated line over the application period:
  
  if (hormone_added != "HFS") {
    treatment_plot <-
      treatment_plot +
      annotate(
        geom = "segment",
        x = hormone_or_HFS_start_time,
        xend = xmax,
        y = ymax,
        yend = ymax,
        colour = line_col,
        linewidth = 0.5
      ) +
      annotate(
        geom = "text",
        x = hormone_or_HFS_start_time,
        y = ymax + 0.06 * ymax,
        label = hormone_added,
        size = if (large_axis_text == "yes") {
          6
        } else {
          4
        },
        hjust = 0,
        family = plot_font_family
      )
  }
  
  # If plot_category = 1 or 3 (experiments involving HFS) add an annotation arrow at 5 minutes
  if (hormone_added == "HFS") {
    # Get limits of x- and y-axes
    ymax <- layer_scales(treatment_plot)$y$get_limits()[2]
    xmax <- layer_scales(treatment_plot)$x$get_limits()[2]
    ymax2 <- layer_scales(treatment_plot)$y$get_limits()[2]
    
    # Add an arrow showing when HFS was applied (x = 5 min)
    treatment_plot <-
      treatment_plot +
      annotate(
        geom = "segment",
        x = hormone_or_HFS_start_time,
        y = ymax + 0.22 * ymax,
        xend = hormone_or_HFS_start_time,
        yend = ymax + 0.10 * ymax,
        arrow = arrow(type = "closed", length = unit(0.02, "npc"))
      )
    
    # Add "HFS" text label
    treatment_plot <-
      treatment_plot +
      annotate(
        geom = "text",
        x = hormone_or_HFS_start_time,
        y = ymax + 0.27 * ymax,
        label = "HFS",
        size = 3.5,
        hjust = 0.5,
        family = plot_font_family
      )
  }
  
  if (signif_stars == "yes") {
    t_test_df <- get(paste0("t_test_", current_type, "_", parameter))
    
    treatment_plot <- treatment_plot +
      geom_text(
        data = t_test_df %>% filter(treatment == plot_treatment),
        aes(
          x = asterisk_time,
          y = y_axis_limit - 50,
          label = significance_stars
        ),
        inherit.aes = FALSE,
        size = 8,
        family = plot_font_family
      )
  }
  
  if (include_representative_trace == "yes") {
    representative_trace_file <- paste0(
      "Figures/Representative-Traces/Category-",
      plot_category,
      "-",
      plot_treatment,
      "-Trace.png"
    )
    
    # Representative traces must be saved as PNGS with the following file name convention:
    # Category-[number]-[treatment]-Trace.png or a warning will display about a missing file
    # e.g. "Category-2-Control-Trace.png"
    
    if (file.exists(here(representative_trace_file))) {
      representative_trace <- png::readPNG(here(representative_trace_file)) %>% rasterGrob()
      
      treatment_plot <- treatment_plot +
        annotation_custom(
          representative_trace,
          xmin = 1,
          xmax = 8,
          ymin = 0,
          ymax = 40
        )
    } else {
      warning(
        "The file here(Figures/Representative-Traces/Category-",
        plot_category,
        "-",
        plot_treatment,
        "-Trace.png) does not exist. Plotting without a representative trace."
      )
    }
  }
  
  
  if (large_axis_text == "yes") {
    text_size <- "_LARGE_text"
  } else {
    text_size <- ""
  }
  
  if (save_plot_pngs == "yes") {
    ggsave(
      treatment_plot,
      path = here(filepath),
      file = paste0(
        "Summary-plot-",
        plot_treatment,
        "-category-",
        plot_category,
        file_name_ending,
        text_size,
        ".png"
      ),
      width = 10,
      height = 7,
      units = "in",
      dpi = 300,
      scaling = 1.25
    )
  }
  treatment_plot
}

make_sEPSC_comparison_plot <-
  function(data,
           plot_category,
           plot_treatment,
           parameter,
           baseline_interval,
           post_hormone_interval,
           large_axis_text = "no",
           hormone_added) {
    if (is.null(baseline_interval) ||
        !is.character(baseline_interval)) {
      stop("'baseline_interval' must be a character (e.g. \"t0to5\" or \"t0to3\")")
    }
    
    if (is.null(post_hormone_interval) ||
        !is.character(post_hormone_interval)) {
      stop("'post_hormone_interval' must be a character (e.g. \"t20to25\")")
    }
    
    allowed_parameters_list <- "\"raw_amplitude\", or \"raw_frequency\""
    
    if (!parameter %in% c("raw_amplitude", "raw_frequency")) {
      stop(
        "parameter must be ",
        allowed_parameters_list,
        " because transformed data are all 100 during the baseline."
      )
    }
    
    plot_colour <- treatment_names_and_colours %>%
      filter(treatment == plot_treatment) %>%
      pull(colours)
    
    sEPSC_comparison_plot_data <- data %>%
      filter(category == plot_category &
               treatment == plot_treatment) %>%
      filter(interval == baseline_interval |
               interval == post_hormone_interval)
    
    if (parameter == "raw_amplitude") {
      y_var <- "mean_raw_amplitude"
      y_title <- "sEPSC Amplitude (pA)"
      
    }
    
    if (parameter == "raw_frequency") {
      y_var <- "mean_raw_frequency"
      y_title <- "sEPSC Frequency (Hz)"
      
    }
    
    sEPSC_comparison_plot <- sEPSC_comparison_plot_data %>%
      ggplot(aes(x = interval, y = .data[[y_var]])) +
      labs(x = NULL, y = y_title) +
      scale_x_discrete(labels = c("Baseline", hormone_added)) +
      geom_violin(
        fill = gray_shading_colour,
        color = NA,
        scale = "width",
        width = 0.2
      ) +
      geom_sina(
        bw = 12,
        alpha = 0.8,
        maxwidth = 0.3,
        size = 2,
        color = plot_colour
      ) +
      geom_signif(
        comparisons = list(c(
          baseline_interval, post_hormone_interval
        )),
        test = "wilcox.test",
        test.args = list(paired = TRUE),
        map_signif_level = list_of_significance_stars,
        vjust = -0.3,
        family = significance_stars_font,
        textsize = geom_signif_text_size,
        size = 0.4
      ) +
      stat_summary(
        fun.data = mean_se,
        geom = "pointrange",
        color = mean_point_colour,
        size = mean_point_size + 0.2,
        alpha = 0.8
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.2, .2)))
    
    if (large_axis_text == "yes") {
      sEPSC_comparison_plot <- sEPSC_comparison_plot +
        theme(
          axis.text.x = element_text(size = 24, margin = margin(t = 10)),
          axis.title.y = element_text(size = 28, face = "plain")
        )
    }
    
    if (save_plot_pngs == "yes") {
      ggsave(
        plot = sEPSC_comparison_plot,
        path = here(
          "Figures/Spontaneous-currents/Baseline-vs-post-hormone-comparisons"
        ),
        file = paste0(
          "Baseline-vs-",
          hormone_added,
          "-",
          parameter,
          "-",
          plot_treatment,
          ".png"
        ),
        width = 7,
        height = 5,
        units = "in",
        dpi = 300
      )
    }
    
    sEPSC_comparison_plot
  }

make_facet_plot <-
  function(data,
           current_type = "eEPSC",
           plot_treatment,
           plot_sex,
           plot_category,
           y_var = P1,
           pruned) {
    if (is.null(current_type) ||
        length(current_type) != 1L ||
        !current_type %in% c("eEPSC")) {
      stop(
        "'current_type' argument must be 'eEPSC' because 'sEPSC' currently takes too long to run and can cause R to crash."
      )
    }
    
    plot_data <- data %>%
      filter(treatment == plot_treatment &
               sex == plot_sex & category == plot_category)
    
    plot_colour <- treatment_names_and_colours %>%
      filter(treatment == plot_treatment) %>%
      pull(colours)
    
    
    plot <- plot_data %>%
      ggplot(aes(x = time, {{ y_var }})) +
      geom_point(color = plot_colour, size = if (pruned == "yes") {
        2.5
      } else {
        1.5
      })
    
    if (pruned == "yes") {
      plot <- plot +
        geom_hline(aes(yintercept = baseline_mean), linetype = "dashed")
    }
    
    plot <- plot +
      modified_facet_theme +
      facet_wrap(. ~ letter, ncol = 3, scales = "free_y") +
      annotate(
        "segment",
        x = -Inf,
        xend = Inf,
        y = -Inf,
        yend = -Inf,
        color = "gray"
      ) +
      annotate(
        "segment",
        x = -Inf,
        xend = -Inf,
        y = -Inf,
        yend = Inf,
        color = "gray"
      ) +
      labs(
        title = paste0(
          unique(plot_data$sex),
          "s, Treatment: ",
          unique(plot_data$treatment)
        ),
        subtitle = paste0("Category: ", unique(plot_data$category)),
        x = "Time (min)"
      )
    
    plot
  }

## Variance analysis ----

make_variance_df <- function(data,
                             df_category,
                             include_all_treatments = "yes",
                             list_of_treatments = NULL,
                             baseline_interval,
                             post_hormone_interval) {
  if (include_all_treatments == "yes") {
    dataframe <- data %>%
      filter(treatment %in% treatment_names_and_colours$treatment) %>%
      droplevels()
    
    treatment_info <- treatment_names_and_colours
    
    if (!is.null(list_of_treatments)) {
      warning(
        "include_all_treatments = \"yes\", but you included a list of treatments to filter. All treatments will be used."
      )
    }
    
  } else {
    if (is.null(list_of_treatments)) {
      stop(
        "include_all_treatments = \"",
        include_all_treatments,
        "\", but list_of_treatments is NULL.",
        "\nDid you forget to add a list of treatments?"
      )
    }
    
    if (!is.character(list_of_treatments)) {
      stop(
        "include_all_treatments = \"",
        include_all_treatments,
        "\", but list_of_treatments is not a character object.",
        "\nDid you forget to add a list of treatments?"
      )
    }
    
    if (is.null(baseline_interval) ||
        !is.character(baseline_interval)) {
      stop("'baseline_interval' must be a character (e.g. \"t0to5\", \"t0to3\")")
    }
    
    if (is.null(post_hormone_interval) ||
        !is.character(post_hormone_interval)) {
      stop("'post_hormone_interval' must be a character (e.g. \"t20to25\")")
    }
    
    dataframe <- data %>%
      filter(treatment %in% list_of_treatments) %>%
      droplevels()
    
    treatment_info <- treatment_names_and_colours %>%
      filter(treatment %in% list_of_treatments)
  }
  
  
  
  variance_df <- dataframe %>%
    filter(category == df_category) %>%
    filter(interval == baseline_interval |
             interval == post_hormone_interval) %>%
    mutate(
      state = case_when(
        interval == baseline_interval ~ "Baseline",
        interval == post_hormone_interval ~ "Post-modification",
        T ~ interval
      )
    ) %>%
    group_by(treatment, state) %>%
    mutate(
      mean_cv_inverse_square = mean(cv_inverse_square),
      mean_VMR = mean(VMR)
    ) %>%
    arrange(match(treatment, treatment_info$display_names))
  
  assign("variance_df", variance_df, envir = .GlobalEnv)
}

# TODO add citation to variance analysis paper ------------
# TODO Add note about how this also uses a wilcox test to add a significance bracket; also add this to PPR plot single treatment -----
#' Make a plot of variance measures for a single treatment. The variance measures are inverse coefficient of variation squared (variance_measure == "cv") or variance-to-mean ratio (variance_measure == "VMR").
#' These can help determine if a mechanism is pre- or post-synaptic 
#' 
#' `make_variance_comparison_plot()` creates a connected scatter plot with 
#' creates a categorical scatter plot with experimental state (i.e. baseline/before and after) on the x-axis
#' and the paired-pulse ratio (PPR) on the y-axis. There are also lines connecting the "before" data point to the "after" data point for each letter.
#' You may customize the baseline and post-modification label to any value if "Baseline" and "Post-hormone" do not work. For example, you may want to use the hormone name instead of "Post-hormone"
#' 


make_variance_comparison_plot <- function(data,
                                          plot_category,
                                          plot_treatment,
                                          large_axis_text = "no",
                                          variance_measure,
                                          baseline_interval,
                                          post_hormone_interval) {
  if (is.null(post_hormone_interval) ||
      !is.character(post_hormone_interval)) {
    stop("'post_hormone_interval' must be a character (e.g. \"t20to25\")")
  }
  
  
  plot_colour <- treatment_names_and_colours %>%
    filter(treatment == plot_treatment) %>%
    pull(colours)
  
  variance_comparison_data <- data %>%
    filter(category == plot_category) %>%
    filter(treatment == plot_treatment)
  
  allowed_parameters_list <- "\"cv\" or \"VMR\""
  
  
  if (!variance_measure %in% c("cv", "VMR")) {
    stop("parameter must be ", allowed_parameters_list)
  }
  
  
  if (variance_measure == "cv") {
    variance_comparison_plot <- variance_comparison_data %>%
      ggplot(aes(x = interval, y = cv_inverse_square, group = letter)) +
      labs(y = "1/CV<sup>2</sup>") +
      geom_signif(
        comparisons = list(c(
          baseline_interval, post_hormone_interval
        )),
        test = "wilcox.test",
        test.args = list(paired = TRUE),
        map_signif_level = list_of_significance_stars,
        vjust = -0.3,
        family = significance_stars_font,
        textsize = geom_signif_text_size,
        size = 0.4
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.2, .2)))
  }
  
  if (variance_measure == "VMR") {
    variance_comparison_plot <- variance_comparison_data %>%
      ggplot(aes(x = interval, y = VMR, group = letter)) +
      labs(y = "VMR") +
      geom_signif(
        comparisons = list(c(
          baseline_interval, post_hormone_interval
        )),
        test = "wilcox.test",
        test.args = list(paired = TRUE),
        map_signif_level = list_of_significance_stars,
        vjust = -0.3,
        family = significance_stars_font,
        textsize = geom_signif_text_size,
        size = 0.4
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.2, .2)))
  }
  
  variance_comparison_plot <- variance_comparison_plot +
    geom_point(color = connecting_line_colour_aps, size = 1.8) +
    geom_line(color = connecting_line_colour_aps, linewidth = 0.4) +
    labs(x = NULL) +
    scale_x_discrete(labels = c("Baseline", "Post-modification")) +
    theme(axis.title.y = element_markdown(family = plot_font_family))
  
  if (variance_measure == "cv") {
    variance_comparison_plot <- variance_comparison_plot +
      annotate(
        geom = "segment",
        x = baseline_interval,
        xend = post_hormone_interval,
        y = variance_comparison_data$mean_cv_inverse_square[variance_comparison_data$interval == baseline_interval][1],
        yend = variance_comparison_data$mean_cv_inverse_square[variance_comparison_data$interval == post_hormone_interval][1],
        color = plot_colour,
        linewidth = 1.2
      ) +
      annotate(
        geom = "point",
        x = baseline_interval,
        y = variance_comparison_data$mean_cv_inverse_square[variance_comparison_data$interval == baseline_interval][1],
        color = plot_colour,
        size = 2.5
      ) +
      annotate(
        geom = "point",
        x = post_hormone_interval,
        y = variance_comparison_data$mean_cv_inverse_square[variance_comparison_data$interval == post_hormone_interval][1],
        color = plot_colour,
        size = 2.5
      )
  }
  
  if (variance_measure == "VMR") {
    variance_comparison_plot <- variance_comparison_plot +
      annotate(
        geom = "segment",
        x = baseline_interval,
        xend = post_hormone_interval,
        y = variance_comparison_data$mean_VMR[variance_comparison_data$interval == baseline_interval][1],
        yend = variance_comparison_data$mean_VMR[variance_comparison_data$interval == post_hormone_interval][1],
        color = plot_colour,
        linewidth = 1.2
      ) +
      annotate(
        geom = "point",
        x = baseline_interval,
        y = variance_comparison_data$mean_VMR[variance_comparison_data$interval == baseline_interval][1],
        color = plot_colour,
        size = 2.5
      ) +
      annotate(
        geom = "point",
        x = post_hormone_interval,
        y = variance_comparison_data$mean_VMR[variance_comparison_data$interval == post_hormone_interval][1],
        color = plot_colour,
        size = 2.5
      )
  }
  
  
  if (large_axis_text == "yes") {
    variance_comparison_plot <- variance_comparison_plot +
      theme(
        axis.text.x = element_text(size = 24, margin = margin(t = 10)),
        axis.title.y = element_markdown(
          size = 28,
          face = "plain",
          family = plot_font_family
        )
      )
  }
  
  if (save_plot_pngs == "yes") {
    ggsave(
      plot = variance_comparison_plot,
      path = here("Figures/Evoked-currents/CV"),
      file = paste0(
        "Variance-comparison-category-",
        plot_category,
        "-",
        plot_treatment,
        "-",
        variance_measure,
        ".png"
      ),
      width = 7,
      height = 5,
      units = "in",
      dpi = 300
    )
  }
  
  variance_comparison_plot
}


# TODO: Fill in required dataframe for cv function --------
# TODO: Find out how to cite variance analysis paper -------

#' Make a plot of coefficient of variation over time
#' 
#' `make_cv_plot()` enables you to save a plot of the coefficient of variation in evoked current amplitudes over time.
#' 
#' @param data A dataframe containing the following columns.
#' @param plot_treatment A character value describing the treatment group. Defaults to "Control".
#' 
#' @export
#' 
#' @returns A ggplot object. If save_plot_PNGs is defined as "yes" in the Global Environment, it will also generate a .png file in the folder `Figures/Evoked-currents/CV` relative to the project directory.
#' The treatment will be included with the filename.
#' 
#' @examples
#' make_cv_plot(
#'   data = pruned_eEPSC_df_all_cells,
#'   plot_treatment = "Control"
#')
#' 
#' @seealso [make_variance_comparison_plot()] to make plots of inverse coefficient of variation squared and VMR, which are useful to determine if a mechanism is pre- or post-synaptic.


make_cv_plot <- function(data, plot_treatment = "Control") {
  plot_colour <- treatment_names_and_colours %>%
    filter(treatment == plot_treatment) %>%
    pull(colours)
  
  cv_plot <- data %>%
    filter(treatment == plot_treatment) %>%
    ggplot(aes(x = time, y = cv_P1_all_cells)) +
    geom_point(color = plot_colour) +
    labs(x = "Time (min)", y = "CV")
  
  if (save_plot_pngs == "yes") {
    ggsave(
      plot = cv_plot,
      path = here("Figures/Evoked-currents/CV"),
      file = paste0("CV_plot", plot_treatment, ".png"),
      width = 7,
      height = 5,
      units = "in",
      dpi = 600
    )
    
  }
  cv_plot
}

## PPR Analysis ----
make_PPR_before_vs_post_hormone_data <- function(data,
                                                 include_all_treatments = "yes",
                                                 list_of_treatments = NULL,
                                                 PPR_min = 0,
                                                 PPR_max = 5,
                                                 baseline_interval,
                                                 post_hormone_interval) {
  if (include_all_treatments == "yes") {
    dataframe <- data %>%
      filter(treatment %in% treatment_names_and_colours$treatment) %>%
      droplevels()
    
    treatment_info <- treatment_names_and_colours
    
    if (!is.null(list_of_treatments)) {
      warning(
        "include_all_treatments = \"yes\", but you included a list of treatments to filter. All treatments will be used."
      )
    }
    
  } else {
    if (is.null(list_of_treatments)) {
      stop(
        "include_all_treatments = \"",
        include_all_treatments,
        "\", but list_of_treatments is NULL.",
        "\nDid you forget to add a list of treatments?"
      )
    }
    
    if (!is.character(list_of_treatments)) {
      stop(
        "include_all_treatments = \"",
        include_all_treatments,
        "\", but list_of_treatments is not a character object.",
        "\nDid you forget to add a list of treatments?"
      )
    }
    
    dataframe <- data %>%
      filter(treatment %in% list_of_treatments) %>%
      droplevels()
    
    treatment_info <- treatment_names_and_colours %>%
      filter(treatment %in% list_of_treatments)
  }
  
  PPR_df <- dataframe %>%
    filter(PPR < PPR_max & PPR > PPR_min) %>%
    filter(interval == baseline_interval |
             interval == post_hormone_interval) %>%
    mutate(
      state = case_when(
        interval == baseline_interval ~ "Baseline",
        interval == post_hormone_interval ~ "Post-modification",
        T ~ interval
      )
    ) %>%
    arrange(match(treatment, treatment_info$display_names))
  assign("PPR_before_vs_post_hormone_data", PPR_df, envir = .GlobalEnv)
}

# TODO: Fill in columns for PPR -------
# TODO: In Vignette describe treatment vs category, see example from make ppr plot for single treatment.

#' Make a PPR plot for a single treatment
#' 
#' `make_PPR_plot_single_treatment()` creates a categorical scatter plot with experimental state (i.e. baseline/before and after) on the x-axis
#' and the paired-pulse ratio (PPR) on the y-axis. There are also lines connecting the "before" data point to the "after" data point for each letter.
#' You may customize the baseline and post-modification label to any value if "Baseline" and "Post-hormone" do not work. For example, you may want to use the hormone name instead of "Post-hormone"
#' 
#' 
#' @param data Paired pulse ratio data, ideally generated from [make_PPR_before_vs_post_hormone_data()]. Must contain the following columns:
#' 
#' \itemize
#'  \item 
#'
#' @param plot_treatment A character value describing the treatment group. Defaults to "Control".
#' @param plot_category A numeric value describing the experimental category. In the sample dataset for this package, 2 represents
#' experiments where insulin was applied continuously after a 5-minute baseline period. Here, `plot_treatment` represents antagonists that were present on the brain slice, or the animals were fasted, etc.
#' @param large_axis_text A character value (must be "yes" or "no") describing if the axis text and labels should be large.
#' If "yes", this activates a ggplot theme layer that increases the axis text size and legend spacing to accommodate the larger text.
#' 
#' @param baseline_label A character value for the x-axis label applied to the pre-hormone state. Defaults to "Baseline".
#' @param post_modification_label A character value for x-axis label applied to the post-hormone or post-protocol state. Defaults to "Post-hormone" but you will likely change this to the hormone or protocol name.
#' 
#' @export
#' 
#' @returns A ggplot object. If save_plot_PNGs is defined as "yes" in the Global Environment, it will also generate a .png file in the folder `Figures/Evoked-currents/PPR` relative to the project directory.
#' The treatment will be included in the filename.
#' 
#' @examples
#' make_PPR_plot_single_treatment(data = PPR_before_vs_post_hormone_data,
#'   plot_treatment = "Control",
#'   plot_category = 2,
#'   large_axis_text = "no",
#'   baseline_label = "Baseline",
#'   post_modification_label = "Insulin")
#'   
#' @seealso [make_PPR_plot_multiple_treatments()] to plot changes in PPR for multiple treatments.

make_PPR_plot_single_treatment <- function(data,
                                           plot_treatment = "Control",
                                           plot_category,
                                           large_axis_text = "no",
                                           baseline_label = "Baseline",
                                           post_modification_label = "Post-hormone") {
  if (!large_axis_text %in% c("yes", "no")) {
    stop("'large_axis_text' argument must be one of: 'yes' or 'no'")
  }
  
  plot_colour <- treatment_names_and_colours %>%
    filter(treatment == plot_treatment) %>%
    pull(colours)
  
  
  PPR_single_plot <- data %>%
    filter(treatment == plot_treatment) %>%
    filter(category == plot_category) %>%
    mutate(state = case_match(state,
                              "Baseline" ~ baseline_label,
                              "Post-modification" ~ post_modification_label)) %>% 
    group_by(treatment, state, letter, sex) %>%
    summarise(mean_PPR_cell = mean(PPR), .groups = "drop") %>%
    ggplot(aes(x = state, y = mean_PPR_cell, shape = sex)) +
    geom_point(
      size = 4,
      color = plot_colour,
      position = position_jitter(width = 0.04, height = 0),
      alpha = 0.8
    ) +
    geom_line(
      aes(group = letter),
      color = plot_colour,
      linewidth = connecting_line_width_PPR,
      alpha = 0.3
    ) +
    stat_summary(
      fun.data = "mean_se",
      geom = "pointrange",
      color = mean_point_colour,
      size = mean_point_size + 0.2,
      alpha = 1,
      position = position_nudge(x = -0.04),
      show.legend = FALSE
    ) +
    theme(legend.position = "right") +
    coord_cartesian(ylim = c(0, 3)) +
    theme(axis.text.x = element_text(margin = margin(b = 5, t = 5))) +
    labs(x = NULL, y = "Paired pulse ratio", shape = NULL) +
    scale_shape_manual(values = c(female_shape, male_shape)) +
    geom_signif(
      comparisons = list(c(baseline_label, post_modification_label)),
      test = "t.test",
      test.args = list(paired = TRUE),
      map_signif_level = list_of_significance_stars,
      vjust = -0.3,
      family = significance_stars_font,
      textsize = geom_signif_text_size,
      size = 0.4,
      y_position = 2.5
    )

  if (large_axis_text == "yes") {
    PPR_single_plot <- PPR_single_plot +
      theme(
        axis.text.x = element_text(size = 24, margin = margin(t = 10)),
        axis.title.y = element_text(size = 28, face = "plain"),
        legend.text = element_text(size = 18),
        legend.key.spacing.y = unit(0.5, "cm")
      )
  }


  if (save_plot_pngs == "yes") {
    ggsave(
      plot = PPR_single_plot,
      path = here("Figures/Evoked-currents/PPR"),
      file = paste0("PPR_comparison-", plot_treatment, ".png"),
      width = 7,
      height = 5,
      units = "in",
      dpi = 300
    )
  }
  PPR_single_plot
}


# TODO: Fill in columns for PPR dataset -------
#' Make a PPR plot for multiple treatments
#' 
#' `make_PPR_plot_multiple_treatments()` creates a categorical scatter plot with experimental state (i.e. grouped as baseline/before and after) and treatment on the x-axis,
#' and the paired-pulse ratio (PPR) on the y-axis. There are also lines connecting the "before" data point to the "after" data point for each letter. It is the same as [make_PPR_plot_multiple_treatments()] but for more than one treatment.
#' You may customize the baseline and post-modification label to any value if "Baseline" and "Post-hormone" do not work. For example, you may want to use the hormone name instead of "Post-hormone"
#' 
#' 
#' @param data Paired pulse ratio data, ideally generated from [make_PPR_before_vs_post_hormone_data()]. Must contain the following columns:
#' 
#' \itemize
#'  \item 
#'
#' @param include_all_treatments A character value (must be "yes" or "no") describing if all treatments in the dataset should be included in the plot. Defaults to "yes", but if "no", you must define a `list_of_treatments`.
#' @param list_of_treatments A list of character values describing the treatments you would like to have in the plot. Examples include `c("Control", "Fasting")`. Defaults to `NULL` since `include_all_treatments` is "yes" by default.
#' @param plot_category A numeric value describing the experimental category. In the sample dataset for this package, 2 represents
#' experiments where insulin was applied continuously after a 5-minute baseline period. Here, `plot_treatment` represents antagonists that were present on the brain slice, or the animals were fasted, etc.
#' @param baseline_label A character value for the x-axis label applied to the pre-hormone state. Defaults to "Baseline".
#' @param post_modification_label A character value for x-axis label applied to the post-hormone or post-protocol state. Defaults to "Post-hormone" but you will likely change this to the hormone or protocol name.
#' 
#' @export
#' 
#' @returns A ggplot object. If save_plot_PNGs is defined as "yes" in the Global Environment, it will also generate a .png file in the folder `Figures/Evoked-currents/PPR` relative to the project directory.
#' 
#' @examples
#' make_PPR_plot_multiple_treatments(data = PPR_before_vs_post_hormone_data,
#'   include_all_treatments = "yes",
#'   plot_category = 2,
#'   baseline_label = "B",
#'   post_modification_label = "I")
#'   
#'   
#' @seealso [make_PPR_plot_single_treatment()] to plot changes in PPR for a single treatment.

make_PPR_plot_multiple_treatments <- function(data,
                                              include_all_treatments = "yes",
                                              list_of_treatments = NULL,
                                              plot_category,
                                              baseline_label,
                                              post_modification_label) {
  if (!include_all_treatments %in% c("yes", "no")) {
    stop("'include_all_treatments' argument must be one of: 'yes' or 'no'")
  }
  
  
  if (include_all_treatments == "yes") {
    treatment_info <- treatment_names_and_colours
    plot_data <- data %>%
      filter(treatment %in% treatment_names_and_colours$treatment) %>%
      droplevels()
    
    if (!is.null(list_of_treatments)) {
      warning(
        "include_all_treatments = \"yes\", but you included a list of treatments to filter. All treatments will be used."
      )
    }
    
  } else {
    if (is.null(list_of_treatments)) {
      stop(
        "include_all_treatments = \"",
        include_all_treatments,
        "\", but list_of_treatments is NULL.",
        "\nDid you forget to add a list of treatments?"
      )
    }
    
    if (!is.character(list_of_treatments)) {
      stop(
        "include_all_treatments = \"",
        include_all_treatments,
        "\", but list_of_treatments is not a character object.",
        "\nDid you forget to add a list of treatments?"
      )
    }
    
    treatment_info <- treatment_names_and_colours %>%
      filter(treatment %in% list_of_treatments)
    plot_data <- data %>%
      filter(treatment %in% list_of_treatments) %>%
      droplevels()
  }
  
  PPR_summary_plot <- plot_data %>%
    filter(category == plot_category) %>%
    mutate(
      state = case_match(state, "Baseline" ~ baseline_label, "Post-modification" ~ post_modification_label),
      treatment = str_replace_all(
        treatment,
        setNames(treatment_info$display_names, treatment_info$treatment)
      ),
      treatment = factor(treatment, levels = treatment_info$display_names)
    ) %>%
    group_by(treatment, state, letter, sex) %>%
    summarize(mean_PPR = mean(PPR)) %>%
    ungroup() %>%
    ggplot(aes(x = state, y = mean_PPR)) +
    geom_point(
      aes(color = treatment),
      size = 2,
      position = position_jitter(0.01),
      alpha = 0.9
    ) +
    geom_line(aes(color = treatment, group = letter),
              linewidth = connecting_line_width_PPR,
              alpha = 0.3) +
    geom_signif(
      comparisons = list(c(baseline_label, post_modification_label)),
      test = "t.test",
      test.args = list(paired = TRUE),
      map_signif_level = list_of_significance_stars,
      vjust = -0.3,
      family = significance_stars_font,
      textsize = 4,
      size = 0.3,
      margin_top = 0.1,
      extend_line = 0.03
    ) +
    facet_wrap(
      ~ treatment,
      ncol = length(treatment_info$treatment),
      strip.position = "bottom"
    ) +
    scale_color_manual(breaks = treatment_info$display_names,
                       values = treatment_info$colours) +
    theme(
      strip.text = element_text(size = 10),
      strip.placement = "outside",
      panel.spacing.x = unit(2, "mm"),
      legend.position = "none"
    ) +
    stat_summary(
      fun.data = "mean_se",
      geom = "pointrange",
      color = mean_point_colour,
      size = mean_point_size + 0.25,
      alpha = 0.9
    ) +
    labs(x = NULL, y = "Paired pulse ratio") +
    scale_y_continuous(expand = expansion(mult = c(0.2, .2)))
  
  if (save_plot_pngs == "yes") {
    ggsave(
      plot = PPR_summary_plot,
      path = here("Figures/Evoked-currents/PPR"),
      file = paste0("PPR_Summary_plot.png"),
      width = 7,
      height = 5,
      units = "in",
      dpi = 300
    )
  }
  PPR_summary_plot
}

# TODO: Find out how to cite another package -----
#' Read and plot raw .abf files ----
#' `import_ABF_file()` is a wrapper around `abf2_load()` from `read_abf`. It converts the array from `abf2_load()`
#' into a dataframe, and it also converts time to minutes
#' 
#' @param file_name Filepath to an .abf file
#' 
#' @returns A dataframe
#' 
#' @export
#' 
#' @examples
#' import_ABF_file("Data/23711004.abf")

import_ABF_file <-
  function(file_name) {
    abf2_load(here(file_name)) %>%
      MeltAbf() %>%
      rename("current" = chan1, "voltage" = chan2) %>%
      rename_with(tolower) %>%
      mutate(time_min = time / 10000)
  }



# TODO: add more time flexibility for the representative trace - and add this to @params------
# TODO: Add more flexibility for the scale bar, and add this to @params
#' Plot a representative spontaneous current trace
#' 
#' `make_representative_sp_trace()` generates a plot of current amplitude over time for a specified sweep from an ABF file.
#' It requires a dataframe generated from raw .abf data with [import_ABF_file()]. The function returns a ggplot object
#' with an optional scale bar.
#' 
#' @param recording_name A dataframe containing at least these columns: `time`, `episode`, `current`, `voltage`, `time_min`.
#' An easy way to obtain this is by importing a raw .abf file through the [import_ABF_file()] function.
#' @param plot_treatment The treatment applied during this recording. This is required to pull the corresponding treatment colour from the `treatment_names_and_colours` dataframe.
#' @param is_baseline A character value that will use the baseline colour (`baseline_group_colour`) if set to "yes". If not, the plot colour is the treatment colour.
#' @param include_scale_bar A character value that determines if a scale bar will be added to the plot. Allowed values are "yes" and "no".
#' @param plot_episode A character value describing the sweep (e.g. `epi1`) that will be used for the plot.
#' 
#' 
#' @returns A ggplot object. If save_plot_PNGs is defined as "yes" in the Global Environment, it will also generate a .png file in the folder `Figures/Spontaneous-currents/Representative-Traces` relative to the project directory.
#' 
#' `Figures/Spontaneous-currents/Representative-Traces`
#' 
#' @export
#' 
#' @examples
#' make_representative_sp_trace(
#'   recording_name = BN_baseline,
#'   plot_treatment = "Control",
#'   is_baseline = "yes",
#'   include_scale_bar = "yes",
#'   plot_episode = "epi1"
#' )
#' 

make_representative_sp_trace <-
  function(recording_name,
           plot_treatment,
           is_baseline = "yes",
           include_scale_bar = "yes",
           plot_episode = "epi1") {
    
    if (!is_baseline %in% c("yes", "no")) {
      stop("'is_baseline' argument must be one of: 'yes' or 'no'")
    }
    
    if (!include_scale_bar %in% c("yes", "no")) {
      stop("'include_scale_bar' argument must be one of: 'yes' or 'no'")
    }
    
    if (is_baseline == "no") {
      plot_colour <- treatment_names_and_colours %>%
        filter(treatment == plot_treatment) %>%
        pull(colours)
    } else {
      plot_colour <- baseline_group_colour
    }
    
    representative_traces_plot <- recording_name %>%
      filter(episode == plot_episode) %>%
      filter(between(time_min, 1, 5)) %>%
      ggplot(aes(x = time_min, y = current)) +
      coord_cartesian(ylim = c(-100, 35)) +
      geom_line(color = plot_colour) +
      theme_void()
    
    
    if (include_scale_bar == "yes") {
      representative_traces_plot <- representative_traces_plot +
        annotate(
          "segment",
          x = 1.25,
          xend = 1.75,
          y = 15,
          yend = 15,
          lwd = 0.4
        ) +
        annotate(
          "segment",
          x = 1.25,
          xend = 1.25,
          y = 15,
          yend = 35,
          lwd = 0.4
        ) +
        annotate(
          "text",
          x = 1.20,
          y = 25,
          label = "20 pA",
          hjust = 1,
          vjust = 0.5,
          family = plot_font_family
        ) +
        annotate(
          "text",
          x = 1.5,
          y = 10,
          label = "30 s",
          hjust = 0.5,
          family = plot_font_family
        )
    }
    
    if (save_plot_pngs == "yes") {
      ggsave(
        plot = representative_traces_plot,
        path = here("Figures/Spontaneous-currents/Representative-Traces"),
        file = paste0(substitute(recording_name), ".png"),
        width = 7,
        height = 5,
        units = "in",
        dpi = 600
      )
    }
    
    representative_traces_plot
  }


# TODO Action Potential Analysis scripts ----

perform_AP_frequency_wilcox_test <- function(category, treatment) {
  ap_frequency_wilcox_test_full <- full_ap_df %>%
    filter(category == category) %>%
    filter(treatment == treatment) %>%
    select(letter,
           State,
           AP_frequency,
           Current_injection,
           category,
           treatment) %>%
    group_by(Current_injection) %>%
    wilcox_test(
      AP_frequency ~ State,
      ref.group = "Baseline",
      paired = T,
      p.adjust.method = "holm"
    ) %>%
    mutate(
      statistic = round(statistic, 2),
      p_string = pvalString(p),
      significance_stars = case_when(
        p == NA ~ "",
        p > 0.05 ~ "",
        0.01 < p & p < 0.05 ~ "*",
        0.001 < p & p < 0.01 ~ "**",
        p < 0.001 ~ "***",
        T ~ ""
      )
    ) %>% merge(
      .,
      max_mean_AP_frequencies %>%
        filter(category == category) %>%
        filter(treatment == treatment),
      by = "Current_injection"
    )
  
  assign(
    paste0(
      "ap_frequency_wilox_test_category_",
      category,
      "_treatment_",
      treatment
    ),
    ap_frequency_wilcox_test_full,
    envir = .GlobalEnv
  )
  
}

make_AP_frequency_plot_single_treatment <- function(data,
                                                    treatment,
                                                    colour_choice,
                                                    signif_stars = "no",
                                                    large_axis_text) {
  single_treatment_AP_plot <- data %>%
    filter(treatment == treatment) %>%
    ggplot(
      aes(
        x = Current_injection,
        y = mean_AP_frequency,
        ymin = mean_AP_frequency - SE,
        ymax = mean_AP_frequency + SE,
        color = State
      )
    ) +
    geom_pointrange(size = 1, linewidth = 0.6) +
    labs(x = "Current Injection (pA)", y = "AP Frequency (Hz)", color = NULL) +
    scale_color_manual(
      values = c(baseline_group_colour, my_colours[colour_choice]),
      labels = c(
        paste0(
          "Baseline, n = ",
          ap_plot_counts %>%
            filter(treatment == treatment & State == "Baseline") %>%
            pull(n)
        ),
        paste0(
          hormone_added,
          ", n = ",
          ap_plot_counts %>%
            filter(treatment == treatment & State == "Insulin") %>%
            pull(n)
        )
      )
    ) +
    theme(
      legend.position = "inside",
      legend.position.inside = c(0.14, 0.8),
      legend.key.spacing.y = unit(1.5, "lines"),
      axis.title = element_text(face = "plain")
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.1, .1)))
  
  
  if (large_axis_text == "yes") {
    single_treatment_AP_plot <- single_treatment_AP_plot +
      theme(
        axis.title.x = element_text(
          size = 28,
          face = "plain",
          margin = margin(t = 10)
        ),
        axis.title.y = element_text(size = 28, face = "plain"),
        legend.text = element_text(size = 18),
        legend.position.inside = c(0.24, 0.8),
        legend.key.spacing.y = unit(0.5, "cm")
      )
  }
  
  
  if (signif_stars == "yes") {
    single_treatment_AP_plot <- single_treatment_AP_plot +
      geom_text(
        data = ap_frequency_wilcox_test_all_treatments %>% filter(treatment == treatment),
        aes(
          x = Current_injection,
          y = max_AP_frequency + max_se + 1,
          label = significance_stars
        ),
        inherit.aes = FALSE,
        size = 5,
        family = significance_stars_font
      )
  }
  
  
  if (save_plot_pngs == "yes") {
    ggsave(
      plot = single_treatment_AP_plot,
      path = here("Figures/Action-potentials"),
      file = paste0("Action-potential-comparison-", treatment, ".png"),
      width = 7,
      height = 5,
      units = "in",
      dpi = 600
    )
    
  }
  single_treatment_AP_plot
}


make_AP_plot <-
  function(data, y, y_axis_title) {
    data %>%
      ggplot(aes(
        x = State,
        y = {
          {
            y
          }
        },
        color = State,
        shape = State
      )) +
      geom_line(aes(group = letter), linewidth = connecting_line_width, color = connecting_line_colour_aps) +
      geom_point(alpha = 0.8,
                 size = geom_sina_size,
                 position = position_jitter(0.02)) +
      #geom_sina(bw = 7, alpha = 0.8, maxwidth = 0.25, size = geom_sina_size) +
      scale_color_manual(values = c(baseline_group_colour, insulin_group_colour)) +
      stat_summary(
        fun.data = mean_se,
        geom = "pointrange",
        color = mean_point_colour,
        size = mean_point_size
      ) +
      theme(
        legend.position = "none",
        axis.line = element_line(linewidth = 0.4),
        axis.title = element_text(face = "plain")
      ) +
      labs(x = NULL, y = y_axis_title) +
      geom_signif(
        comparisons = list(c("Baseline", "Insulin")),
        test = "wilcox.test",
        test.args = list(paired = TRUE),
        map_signif_level = list_of_significance_stars,
        vjust = -0.3,
        family = significance_stars_font,
        textsize = geom_signif_text_size,
        size = 0.4
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.2, .2)))
  }


save_AP_plot <- function(ap_plot, plot_name) {
  ggsave(
    plot = ap_plot,
    path = here("Figures/Action-potentials"),
    file = paste0(plot_name, "-summary-plot.png"),
    width = 7,
    height = 5,
    units = "in",
    dpi = 300
  )
}


make_qq_plot <- function(data, parameter_title, parameter) {
  ggplot(data, aes(sample = {
    {
      parameter
    }
  })) +
    stat_qq() +
    stat_qq_line() +
    labs(title = glue("QQ-Plot for {parameter_title}"))
}

# coord_cartesian() is required to lock' the coordinates of all AP traces to the same scales
make_AP_trace_plot <-
  function(file_ID, sweep1, sweep2, trace_color) {
    abf2_load(here(paste0("Data/ABF-Files/", file_ID, ".abf"))) %>%
      MeltAbf() %>%
      rename("Voltage" = chan1, "Current" = chan2) %>%
      filter(Episode %in% c(sweep1, sweep2)) %>%
      ggplot(aes(x = time, y = Voltage, group = Episode)) +
      geom_line(color = trace_color, linewidth = AP_trace_size) +
      labs(x = NULL, y = NULL) +
      coord_cartesian(
        xlim  = c(ap_traces_x_min, ap_traces_x_max),
        ylim = c(ap_traces_y_min, ap_traces_y_max)
      ) +
      theme_void()
  }