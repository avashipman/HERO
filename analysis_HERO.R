# HERO: Hierarchy-aware Expression Region Organization
# User-facing run script for LifeCanvas SmartAnalytics cell detection outputs.
#
# INSTRUCTIONS:
# 1) Put this file in the same folder as functions_HERO.R, or edit functions_file below.
# 2) Edit the paths in Section 1.
# 3) Run Sections 2-6 as needed.

# -----------------------------------------------------------------------------
# 0. Load HERO functions
# -----------------------------------------------------------------------------
functions_file <- "path/to/functions_HERO.R"
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
original_csv_folder <- "path/to/folder/with/cell-detection_datasets"

# Output folders/files for this group. Names can be changed to user's desire
analysis_folder <- "path/to/folder/for/HERO/output/files"
filtered_folder <- file.path(analysis_folder, "nameoffolderforHEROfilteredstep")
averaged_file <- file.path(analysis_folder, "AveragedRegions_Depth6.xlsx")
uncleaned_top30_file <- file.path(analysis_folder, "Top30RegionsDepth6.xlsx")
cleaned_top30_file <- file.path(analysis_folder, "Cleaned_Top30RegionsDepth6.xlsx")

# Analysis settings.
depth_start <- 4              # Starting parent depth; Depth 4 is recommended as a broad first pass.
target_depth <- 6             # Target child depth for detailed region selection.
laterality <- "combined"      # Choose "combined", "left", or "right".
top_n_parents <- 30           # Number of parent regions to expand.
top_n_final <- 30             # Number of final regions to report.
collapse_threshold <- 2       # Collapse parent if >3 of its children appear in the final list.
exclude_terms <- c("tract", "fascicle", "nerve", "tracts", "nerves")   # Remove regions containing these terms.

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
# 4. Select: uncleaned hierarchy-aware top region selection
# -----------------------------------------------------------------------------
# This creates Top30RegionsDepth6.xlsx, which contains the raw top target-depth
# regions before overrepresented parent regions are collapsed.
Top30Regions(
  averaged_file = averaged_file,
  original_csv_folder = original_csv_folder,
  filtered_folder_for_id_map = filtered_folder,
  output_file = uncleaned_top30_file,
  depth_start = depth_start,
  target_depth = target_depth,
  top_n_parents = top_n_parents,
  top_n_final = top_n_final,
  laterality = laterality,
  exclude_terms = exclude_terms
)

# -----------------------------------------------------------------------------
# 5. Clean: collapse overrepresented parents and refill final list
# -----------------------------------------------------------------------------
# This creates Top30RegionsDepth6_Cleaned.xlsx, which applies the HERO
# collapse/refill logic to reduce parent-child overrepresentation.
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
  laterality = laterality,
  exclude_terms = exclude_terms
)

# -----------------------------------------------------------------------------
# 6. Optional focused analysis: start from one or more user-defined parent regions
# -----------------------------------------------------------------------------
# Example: analyze descendants from selected parent regions only.
focused_file <- file.path(analysis_folder, "Top30_From_SelectedRegions_Cleaned.xlsx")
Top30FocusedCleaned(
   original_csv_folder = original_csv_folder,
   start_region = c("Striatum", "Midbrain, motor related", "Pallidum", "Hypothalamus", "Cortical Plate", "Cortical Subplate", "Hippocampal Formation"),
   output_file = focused_file,
   start_depth = 4,
   target_depth = 6,
   top_n_final = top_n_final,
   collapse_threshold = collapse_threshold,
   laterality = "combined"
 )

# -----------------------------------------------------------------------------
# 7. Optional visualization: chord diagram comparing multiple groups
# -----------------------------------------------------------------------------
# To compare groups, run Sections 1-4 separately for each group, then point to the
# resulting Top30RegionsDepth6_Cleaned.xlsx files below.
#
make_density_chord_parentgroups(
  group_files = c("path/to/group1/file/to/visualize",
                  "path/to/group2/file/to/visualize"),
  group_names = c("Insula_rg", "BST_rg"),
  sheet_avg = "Averages",
  top_regions_per_group = 30,
  region_label_cex = 0.5,
  show_parent_labels = FALSE,
  show_legends = FALSE,
  unique_group_fills = c(
  Insula_rg = "#00b0be",
  BST_rg = "#bdd373"
  )
)
# -----------------------------------------------------------------------------
# End of HERO run script
# -----------------------------------------------------------------------------
