# HERO: Hierarchy-aware Expression Region Organization
# User-facing run script for whole-brain light-sheet SmartAnalytics outputs.
#
# INSTRUCTIONS:
# 1) Put this file in the same folder as functions_HERO.R, or edit functions_file below.
# 2) Edit the paths in Section 1.
# 3) Run Sections 2-6 as needed.

# -----------------------------------------------------------------------------
# 0. Load HERO functions
# -----------------------------------------------------------------------------
functions_file <- "GitHubScripts/functions_HERO.R"
if (!file.exists(functions_file)) {
  # If running from a different working directory, use an absolute path instead.
  stop("Could not find functions_HERO.R. Edit functions_file to the correct path.")
}
source(functions_file)

# Install packages once if needed:
# install.packages(c("readxl", "writexl", "dplyr", "tidyr", "stringr", "circlize"))
HERO_check_packages(c("readxl", "writexl", "dplyr", "tidyr", "stringr"))

# -----------------------------------------------------------------------------
# 1. User inputs: edit these paths and settings
# -----------------------------------------------------------------------------

# Folder containing raw SmartAnalytics CSV files for ONE experimental group.
original_csv_folder <- "GitHubScripts/PracticeScriptRun/Original CSV folder"

# Output folders/files for this group.
analysis_folder <- "GitHubScripts/PracticeScriptRun/HERO analysis"
filtered_folder <- file.path(analysis_folder, "Filtered")
averaged_file <- file.path(analysis_folder, "AveragedRegions.xlsx")
cleaned_top30_file <- file.path(analysis_folder, "Top30RegionsDepth6_Cleaned.xlsx")

# Analysis settings.
depth_start <- 4              # Starting parent depth; Depth 4 is recommended as a broad first pass.
target_depth <- 6             # Target child depth for detailed region selection.
laterality <- "combined"      # Choose "combined", "left", or "right".
top_n_parents <- 20           # Number of parent regions to expand.
top_n_final <- 30             # Number of final regions to report.
collapse_threshold <- 3       # Collapse parent if >3 of its children appear in the final list.
exclude_terms <- c("tract")   # Remove regions containing these terms.

# -----------------------------------------------------------------------------
# 2. Prepare: filter raw SmartAnalytics CSVs by depth and laterality
# -----------------------------------------------------------------------------
FilterbyDepth(
  input_folder = original_csv_folder,
  output_folder = filtered_folder,
  depth_level = depth_start,
  laterality = laterality,
  exclude_terms = exclude_terms
)

# -----------------------------------------------------------------------------
# 3. Average: compute group-level region densities
# -----------------------------------------------------------------------------
RegionsbyDensity(
  filtered_folder = filtered_folder,
  output_file = averaged_file
)

# -----------------------------------------------------------------------------
# 4. Select + clean: hierarchy-aware top region selection
# -----------------------------------------------------------------------------
Top30Cleaned(
  averaged_file = averaged_file,
  original_csv_folder = original_csv_folder,
  filtered_folder_for_id_map = filtered_folder,
  output_file = cleaned_top30_file,
  depth_start = depth_start,
  target_depth = target_depth,
  top_n_parents = top_n_parents,
  top_n_final = top_n_final,
  collapse_threshold = collapse_threshold,
  laterality = laterality
)

# -----------------------------------------------------------------------------
# 5. Optional focused analysis: start from one or more user-defined parent regions
# -----------------------------------------------------------------------------
# Example: analyze descendants from selected parent regions only.
# focused_file <- file.path(analysis_folder, "Top30_From_SelectedRegions_Cleaned.xlsx")
# Top30FromRegionCleaned(
#   original_csv_folder = original_csv_folder,
#   start_region = c("Hypothalamus", "Pallidum"),
#   start_depth = c(Hypothalamus = 3, Pallidum = 4),
#   output_file = focused_file,
#   target_depth = target_depth,
#   top_n_final = top_n_final,
#   collapse_threshold = collapse_threshold,
#   laterality = laterality
# )

# -----------------------------------------------------------------------------
# 6. Optional visualization: chord diagram comparing multiple groups
# -----------------------------------------------------------------------------
# To compare groups, run Sections 1-4 separately for each group, then point to the
# resulting Top30RegionsDepth6_Cleaned.xlsx files below.
#
# group_files <- c(
#   "path/to/GroupA/Top30RegionsDepth6_Cleaned.xlsx",
#   "path/to/GroupB/Top30RegionsDepth6_Cleaned.xlsx"
# )
# group_names <- c("Group A", "Group B")
#
# make_density_chord_parentgroups(
#   group_files = group_files,
#   group_names = group_names,
#   top_regions_per_group = 15,  # plots the top N selected within each group
#   show_parent_labels = FALSE,
#   show_legends = TRUE
# )

# -----------------------------------------------------------------------------
# End of HERO run script
# -----------------------------------------------------------------------------
