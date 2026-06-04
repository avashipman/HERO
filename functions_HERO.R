# HERO: Hierarchy-aware Expression Region Organization
# Functions for hierarchy-aware analysis of SmartAnalytics/LifeCanvas whole-brain light-sheet expression CSVs.

# This file contains reusable functions only. Edit and run analysis_HERO.R to apply
# these functions to your own data.
#
# Author: Ava Shipman / Centanni Lab

# -----------------------------------------------------------------------------
# Required packages
# -----------------------------------------------------------------------------
# HERO uses: readxl, writexl, dplyr, tidyr, stringr, circlize
# Install missing packages with:
# install.packages(c("readxl", "writexl", "dplyr", "tidyr", "stringr", "circlize"))

HERO_check_packages <- function(pkgs = c("readxl", "writexl", "dplyr", "tidyr", "stringr")) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Install required packages first: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

HERO_strip_lr <- function(x) {
  x <- stringr::str_replace(x, stringr::regex("^left\\s+",  ignore_case = TRUE), "")
  x <- stringr::str_replace(x, stringr::regex("^right\\s+", ignore_case = TRUE), "")
  x
}

HERO_strip_acr_lr <- function(x) {
  stringr::str_replace(x, "(-L|-R)$", "")
}

HERO_pick_region_col <- function(df, prefs = c("region_base", "Region", "region", "name")) {
  for (nm in prefs) if (nm %in% names(df)) return(nm)
  stop("None of the preferred region columns exist. Looked for: ",
       paste(prefs, collapse = ", "), ". Found columns: ", paste(names(df), collapse = ", "), call. = FALSE)
}

#' Build hierarchy table for a SmartAnalytics CSV

#' Internal helper used by HERO functions. Reads one SmartAnalytics CSV, detects
#' hemisphere from the original name/acronym, optionally merges left/right regions,
#' recalculates combined density, and rebuilds parent-child IDs in combined mode.

#' @param csv_file Path to one SmartAnalytics CSV.
#' @param laterality One of "combined", "left", or "right".
#' @param density_col Name of density column.
#' @param exclude_terms Character vector of terms to exclude from region names.
#' @return Data frame with id, parent_structure_id, depth, region_base, acronym_base, density.
HERO_build_hierarchy <- function(csv_file,
                                 laterality = c("combined", "left", "right"),
                                 density_col = "density (cells/mm^3)",
                                 exclude_terms = c("tract")) {
  laterality <- match.arg(laterality)
  HERO_check_packages(c("dplyr", "stringr"))

  df <- read.csv(csv_file, stringsAsFactors = FALSE, check.names = FALSE)
  req <- c("id", "name", "acronym", "parent_structure_id", "depth", "count", "volume (mm^3)", density_col)
  miss <- setdiff(req, names(df))
  if (length(miss) > 0) stop("Missing columns in ", basename(csv_file), ": ", paste(miss, collapse = ", "))

  df <- df[!is.na(df$depth) & df$depth >= 0, , drop = FALSE]
  if (length(exclude_terms) > 0) {
    pattern <- paste(exclude_terms, collapse = "|")
    df <- df[!grepl(pattern, df$name, ignore.case = TRUE), , drop = FALSE]
  }

  hemi_from_name <- rep(NA_character_, nrow(df))
  hemi_from_name[grepl("^left\\s+",  df$name, ignore.case = TRUE)] <- "left"
  hemi_from_name[grepl("^right\\s+", df$name, ignore.case = TRUE)] <- "right"

  hemi_from_acr <- rep(NA_character_, nrow(df))
  hemi_from_acr[grepl("-L$", df$acronym)] <- "left"
  hemi_from_acr[grepl("-R$", df$acronym)] <- "right"

  df$hemisphere <- dplyr::coalesce(hemi_from_name, hemi_from_acr, "combined")
  df$region_base  <- HERO_strip_lr(df$name)
  df$acronym_base <- HERO_strip_acr_lr(df$acronym)

  id_map <- df[, c("id", "region_base", "depth")]
  names(id_map) <- c("id", "parent_region_base_lookup", "parent_depth_lookup")
  df$parent_structure_id_chr <- as.character(df$parent_structure_id)
  id_map$id_chr <- as.character(id_map$id)

  df <- dplyr::left_join(
    df,
    id_map[, c("id_chr", "parent_region_base_lookup", "parent_depth_lookup")],
    by = c("parent_structure_id_chr" = "id_chr")
  )
  df$parent_base <- df$parent_region_base_lookup

  if (laterality %in% c("left", "right")) {
    return(df[df$hemisphere == laterality, , drop = FALSE] |>
             dplyr::transmute(
               id = as.integer(id),
               parent_structure_id = as.integer(parent_structure_id),
               depth = as.integer(depth),
               region_base = region_base,
               acronym_base = acronym_base,
               density = as.numeric(.data[[density_col]]),
               count = as.numeric(count),
               volume = as.numeric(`volume (mm^3)`)
             ))
  }

  df_lr <- df[df$hemisphere %in% c("left", "right"), , drop = FALSE]
  df_cb <- df[df$hemisphere == "combined", , drop = FALSE]
  src <- if (nrow(df_lr) > 0) df_lr else df_cb

  src$count <- as.numeric(src$count)
  src[["volume (mm^3)"]] <- as.numeric(src[["volume (mm^3)"]])

  merged <- src |>
    dplyr::mutate(parent_base = dplyr::coalesce(parent_base, NA_character_)) |>
    dplyr::group_by(depth, parent_base, region_base, acronym_base) |>
    dplyr::summarise(
      count = sum(count, na.rm = TRUE),
      volume = sum(`volume (mm^3)`, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(density = dplyr::if_else(volume > 0, count / volume, NA_real_)) |>
    dplyr::arrange(depth, parent_base, region_base) |>
    dplyr::mutate(id = dplyr::row_number())

  parent_lookup <- merged |>
    dplyr::select(id, depth, region_base) |>
    dplyr::rename(parent_id = id,
                  parent_depth = depth,
                  parent_region_base = region_base)

  merged |>
    dplyr::left_join(parent_lookup, by = c("parent_base" = "parent_region_base")) |>
    dplyr::mutate(
      parent_structure_id = dplyr::if_else(!is.na(parent_base) & parent_depth == depth - 1,
                                           parent_id, NA_integer_)
    ) |>
    dplyr::select(id, parent_structure_id, depth, region_base, acronym_base, density, count, volume)
}

#' Find descendants at a specified atlas depth
#'
#' @param hier Hierarchy table from HERO_build_hierarchy().
#' @param start_ids Starting parent ID(s).
#' @param want_depth Target atlas depth to return.
#' @return Integer vector of descendant IDs at want_depth.
HERO_descendants_at_depth <- function(hier, start_ids, want_depth) {
  start_ids <- unique(as.integer(start_ids))
  start_ids <- start_ids[!is.na(start_ids)]
  if (length(start_ids) == 0) return(integer(0))

  children_map <- split(hier$id, hier$parent_structure_id)
  current <- start_ids

  repeat {
    depths_now <- unique(hier$depth[hier$id %in% current])
    if (any(depths_now == want_depth, na.rm = TRUE)) break

    nxt <- unlist(children_map[as.character(current)], use.names = FALSE)
    nxt <- unique(as.integer(nxt))
    nxt <- nxt[!is.na(nxt)]
    if (length(nxt) == 0) break

    current <- nxt
    if (all(hier$depth[hier$id %in% current] >= want_depth, na.rm = TRUE)) break
  }

  unique(hier$id[hier$depth == want_depth & hier$id %in% current])
}

#' Rank regions by mean density
HERO_rank_table <- function(long_df, include_parent_cols = FALSE) {
  HERO_check_packages("dplyr")
  if (nrow(long_df) == 0) {
    if (include_parent_cols) {
      return(dplyr::tibble(Region=character(), Acronym=character(), Depth=integer(), MeanDensity=numeric(),
                           ParentRegion=character(), ParentAcronym=character(), ParentDepth=integer()))
    }
    return(dplyr::tibble(Region=character(), Acronym=character(), Depth=integer(), MeanDensity=numeric()))
  }
  if (include_parent_cols) {
    long_df |>
      dplyr::group_by(Region) |>
      dplyr::summarise(
        MeanDensity = mean(Density, na.rm = TRUE),
        Acronym = dplyr::first(stats::na.omit(Acronym)),
        Depth = dplyr::first(Depth),
        ParentRegion = dplyr::first(stats::na.omit(ParentRegion)),
        ParentAcronym = dplyr::first(stats::na.omit(ParentAcronym)),
        ParentDepth = dplyr::first(stats::na.omit(ParentDepth)),
        .groups = "drop"
      ) |>
      dplyr::arrange(dplyr::desc(MeanDensity))
  } else {
    long_df |>
      dplyr::group_by(Region) |>
      dplyr::summarise(
        MeanDensity = mean(Density, na.rm = TRUE),
        Acronym = dplyr::first(stats::na.omit(Acronym)),
        Depth = dplyr::first(Depth),
        .groups = "drop"
      ) |>
      dplyr::arrange(dplyr::desc(MeanDensity))
  }
}

# -----------------------------------------------------------------------------
# Function 1: FilterbyDepth
# -----------------------------------------------------------------------------
#' Filter SmartAnalytics CSVs by atlas depth and laterality
#'
#' Reads all CSV files in an input folder, cleans laterality, optionally merges
#' left and right hemisphere structures, removes excluded terms such as tracts,
#' and writes one filtered Excel file per sample.
#'
#' @param input_folder Folder containing raw SmartAnalytics CSV files.
#' @param output_folder Folder where filtered Excel files should be written.
#' @param depth_level Atlas depth to retain (e.g., 4). Shallow leaf regions are also retained.
#' @param laterality One of "combined", "left", or "right".
#' @param density_col Name of density column.
#' @param exclude_terms Terms to exclude from region names; default excludes "tract".
#' @return Invisibly returns paths to written Excel files.
FilterbyDepth <- function(input_folder,
                          output_folder,
                          depth_level = 4,
                          laterality = c("combined", "left", "right"),
                          density_col = "density (cells/mm^3)",
                          exclude_terms = c("tract")) {
  laterality <- match.arg(laterality)
  HERO_check_packages(c("writexl", "dplyr", "stringr"))
  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

  file_list <- list.files(input_folder, pattern = "\\.csv$", full.names = TRUE)
  if (length(file_list) == 0) stop("No .csv files found in input_folder: ", input_folder)

  out_paths <- character()
  for (filename in file_list) {
    hier <- HERO_build_hierarchy(filename, laterality = laterality, density_col = density_col, exclude_terms = exclude_terms)
    child_count <- table(as.character(hier$parent_structure_id))
    hier$has_children <- as.character(hier$id) %in% names(child_count)

    output_df <- hier |>
      dplyr::filter(depth == depth_level | (depth < depth_level & !has_children)) |>
      dplyr::arrange(depth, region_base)

    base_name <- tools::file_path_sans_ext(basename(filename))
    out_file <- file.path(output_folder, paste0(base_name, "_FilteredDepth", depth_level, "_", laterality, ".xlsx"))
    writexl::write_xlsx(output_df, out_file)
    out_paths <- c(out_paths, out_file)
  }
  invisible(out_paths)
}

# -----------------------------------------------------------------------------
# Function 2: RegionsbyDensity
# -----------------------------------------------------------------------------
#' Average region density across samples
#'
#' Reads filtered Excel files and produces AveragedRegions.xlsx with sample-level
#' densities and region-level means.
#'
#' @param filtered_folder Folder containing FilterbyDepth Excel files.
#' @param output_file Output Excel path or filename.
#' @param density_col Density column name.
#' @param region_col_preference Candidate region-name columns, checked in order.
#' @param id_from_filename If TRUE, sample IDs are derived from filenames.
#' @return Invisibly returns output path.
RegionsbyDensity <- function(filtered_folder,
                             output_file = "AveragedRegions.xlsx",
                             density_col = "density",
                             region_col_preference = c("region_base", "name"),
                             id_from_filename = TRUE) {
  HERO_check_packages(c("readxl", "writexl", "dplyr", "tidyr"))

  xlsx_files <- list.files(filtered_folder, pattern = "\\.(xlsx|xls)$", full.names = TRUE)
  if (length(xlsx_files) == 0) stop("No Excel files found in filtered_folder: ", filtered_folder)

  sample_long_list <- list()
  meta_long_list <- list()
  for (f in xlsx_files) {
    df <- readxl::read_excel(f)
    region_col <- HERO_pick_region_col(df, region_col_preference)
    if (!density_col %in% names(df)) stop("Density column '", density_col, "' not found in: ", basename(f))

    sample_id <- tools::file_path_sans_ext(basename(f))
    if (!id_from_filename && "sample_id" %in% names(df)) sample_id <- as.character(df$sample_id[[1]])

    tmp <- df |>
      dplyr::transmute(sample_id = sample_id,
                       region = as.character(.data[[region_col]]),
                       density = as.numeric(.data[[density_col]])) |>
      dplyr::filter(!is.na(region) & region != "") |>
      dplyr::group_by(sample_id, region) |>
      dplyr::summarise(density = mean(density, na.rm = TRUE), .groups = "drop")
    sample_long_list[[length(sample_long_list)+1]] <- tmp

    if (all(c("id", "depth") %in% names(df))) {
      meta_long_list[[length(meta_long_list)+1]] <- df |>
        dplyr::transmute(sample_id = sample_id,
                         region = as.character(.data[[region_col]]),
                         id = suppressWarnings(as.integer(id)),
                         depth = suppressWarnings(as.integer(depth))) |>
        dplyr::filter(!is.na(region) & region != "") |>
        dplyr::group_by(sample_id, region) |>
        dplyr::summarise(id = dplyr::first(id[!is.na(id)]),
                         depth = dplyr::first(depth[!is.na(depth)]),
                         .groups = "drop")
    }
  }

  long_all <- dplyr::bind_rows(sample_long_list)
  wide <- long_all |> tidyr::pivot_wider(names_from = region, values_from = density)
  region_means <- wide |>
    dplyr::select(-sample_id) |>
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ mean(.x, na.rm = TRUE)))
  means_vec <- as.numeric(region_means[1, ])
  names(means_vec) <- names(region_means)
  ordered_regions <- names(sort(means_vec, decreasing = TRUE))

  wide_ordered <- wide |> dplyr::select(sample_id, dplyr::all_of(ordered_regions))
  averages_df <- data.frame(region = ordered_regions,
                            mean_density = unname(means_vec[ordered_regions]),
                            stringsAsFactors = FALSE)

  if (length(meta_long_list) > 0) {
    meta_region <- dplyr::bind_rows(meta_long_list) |>
      dplyr::filter(!is.na(id) & !is.na(depth)) |>
      dplyr::group_by(region, id, depth) |>
      dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
      dplyr::arrange(region, dplyr::desc(n)) |>
      dplyr::group_by(region) |>
      dplyr::slice_head(n = 1) |>
      dplyr::ungroup() |>
      dplyr::select(region, id, depth)
    averages_df <- averages_df |> dplyr::left_join(meta_region, by = "region")
  }

  out_path <- if (grepl("\\.(xlsx|xls)$", output_file, ignore.case = TRUE)) output_file else file.path(filtered_folder, output_file)
  writexl::write_xlsx(list("Densities" = wide_ordered, "Averages" = averages_df), out_path)
  invisible(out_path)
}


# -----------------------------------------------------------------------------
# Function 3: Top30Cleaned
# -----------------------------------------------------------------------------
#' Select and clean top hierarchy-aware expression regions
#'
#' Starting from top parent regions in AveragedRegions.xlsx, this function traverses
#' the atlas hierarchy to a target depth, ranks candidate regions by mean density,
#' collapses overrepresented parents (> collapse_threshold children), refills from
#' additional descendants, and writes a multi-tab Excel workbook.
#'
#' @param averaged_file Path to AveragedRegions.xlsx from RegionsbyDensity().
#' @param original_csv_folder Folder containing original SmartAnalytics CSVs.
#' @param filtered_folder_for_id_map Folder containing FilterbyDepth outputs.
#' @param output_file Output Excel file path.
#' @param depth_start Starting atlas depth used for parent selection.
#' @param target_depth Target descendant depth, usually 6.
#' @param top_n_parents Number of starting parent regions to expand.
#' @param top_n_final Number of final regions to report.
#' @param collapse_threshold Collapse a parent if more than this number of children appear.
#' @param laterality One of "combined", "left", or "right".
#' @param density_col Name of density column in raw CSVs.
#' @param max_iter Maximum collapse/refill iterations.
#' @return Invisibly returns a list with final tables.
#Function 5: collapses “overrepresented” parents: if > 3 children from the same 
#parent appear in the top-30 list, replace those children with the parent (depth-5), 
#then backfill remaining slots with the next-best depth-6 regions (excluding 
#children of collapsed parents).
#Things to change/update:
  #Rank regions within subject density
  #Add median as well as averages

Top30Cleaned <- function(averaged_file,
                         original_csv_folder,
                         filtered_folder_for_id_map,
                         output_file = file.path(dirname(averaged_file), "Top30RegionsDepth6_CleanedGFP.xlsx"),
                         depth_start = 4,
                         target_depth = 6,
                         top_n_parents = 20,
                         top_n_final = 30,
                         collapse_threshold = 3,     # > 3 children => collapse
                         laterality = c("combined", "left", "right"),
                         density_col = "density (cells/mm^3)",
                         max_iter = 10) {
  
  laterality <- match.arg(laterality)
  
  if (!requireNamespace("readxl", quietly = TRUE)) stop("Install 'readxl'")
  if (!requireNamespace("writexl", quietly = TRUE)) stop("Install 'writexl'")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Install 'dplyr'")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Install 'tidyr'")
  if (!requireNamespace("stringr", quietly = TRUE)) stop("Install 'stringr'")
  
  strip_lr <- function(x) {
    x <- stringr::str_replace(x, stringr::regex("^left\\s+",  ignore_case = TRUE), "")
    x <- stringr::str_replace(x, stringr::regex("^right\\s+", ignore_case = TRUE), "")
    x
  }
  strip_acr_lr <- function(x) stringr::str_replace(x, "(-L|-R)$", "")
  
  build_hierarchy <- function(csv_file) {
    df <- read.csv(csv_file, stringsAsFactors = FALSE, check.names = FALSE)
    
    req <- c("id","name","acronym","parent_structure_id","depth","count","volume (mm^3)", density_col)
    miss <- setdiff(req, names(df))
    if (length(miss) > 0) stop("Missing columns in ", basename(csv_file), ": ", paste(miss, collapse=", "))
    
    df <- df[!is.na(df$depth) & df$depth >= 0, , drop = FALSE]
    
    hemi_from_name <- rep(NA_character_, nrow(df))
    hemi_from_name[grepl("^left\\s+",  df$name, ignore.case=TRUE)] <- "left"
    hemi_from_name[grepl("^right\\s+", df$name, ignore.case=TRUE)] <- "right"
    
    hemi_from_acr <- rep(NA_character_, nrow(df))
    hemi_from_acr[grepl("-L$", df$acronym)] <- "left"
    hemi_from_acr[grepl("-R$", df$acronym)] <- "right"
    
    df$hemisphere <- dplyr::coalesce(hemi_from_name, hemi_from_acr, "combined")
    df$region_base  <- strip_lr(df$name)
    df$acronym_base <- strip_acr_lr(df$acronym)
    
    # parent_base lookup
    id_map <- df[, c("id","region_base","depth")]
    names(id_map) <- c("id","parent_region_base_lookup","parent_depth_lookup")
    df$parent_structure_id_chr <- as.character(df$parent_structure_id)
    id_map$id_chr <- as.character(id_map$id)
    
    df <- dplyr::left_join(
      df,
      id_map[, c("id_chr","parent_region_base_lookup","parent_depth_lookup")],
      by = c("parent_structure_id_chr" = "id_chr")
    )
    df$parent_base <- df$parent_region_base_lookup
    
    if (laterality %in% c("left","right")) {
      return(df[df$hemisphere == laterality, , drop=FALSE] |>
               dplyr::transmute(
                 id = as.integer(id),
                 parent_structure_id = as.integer(parent_structure_id),
                 depth = as.integer(depth),
                 region_base = region_base,
                 acronym_base = acronym_base,
                 density = as.numeric(.data[[density_col]])
               ))
    }
    
    # Combined: merge L/R, recompute density, rebuild parent ids
    df_lr <- df[df$hemisphere %in% c("left","right"), , drop=FALSE]
    df_cb <- df[df$hemisphere == "combined", , drop=FALSE]
    src <- if (nrow(df_lr) > 0) df_lr else df_cb
    
    src$count <- as.numeric(src$count)
    src[["volume (mm^3)"]] <- as.numeric(src[["volume (mm^3)"]])
    
    merged <- src |>
      dplyr::mutate(parent_base = dplyr::coalesce(parent_base, NA_character_)) |>
      dplyr::group_by(depth, parent_base, region_base, acronym_base) |>
      dplyr::summarise(
        count = sum(count, na.rm=TRUE),
        volume = sum(`volume (mm^3)`, na.rm=TRUE),
        .groups = "drop"
      ) |>
      dplyr::mutate(density = dplyr::if_else(volume > 0, count/volume, NA_real_)) |>
      dplyr::arrange(depth, parent_base, region_base) |>
      dplyr::mutate(id = dplyr::row_number())
    
    parent_lookup <- merged |>
      dplyr::select(id, depth, region_base) |>
      dplyr::rename(parent_id = id,
                    parent_depth = depth,
                    parent_region_base = region_base)
    
    merged |>
      dplyr::left_join(parent_lookup, by = c("parent_base" = "parent_region_base")) |>
      dplyr::mutate(
        parent_structure_id = dplyr::if_else(!is.na(parent_base) & parent_depth == depth - 1,
                                             parent_id, NA_integer_)
      ) |>
      dplyr::select(id, parent_structure_id, depth, region_base, acronym_base, density)
  }
  
  descendants_at_depth <- function(hier, start_ids, want_depth) {
    start_ids <- unique(as.integer(start_ids))
    start_ids <- start_ids[!is.na(start_ids)]
    if (length(start_ids) == 0) return(integer(0))
    
    children_map <- split(hier$id, hier$parent_structure_id)
    current <- start_ids
    
    repeat {
      depths_now <- unique(hier$depth[hier$id %in% current])
      if (any(depths_now == want_depth, na.rm=TRUE)) break
      
      nxt <- unlist(children_map[as.character(current)], use.names=FALSE)
      nxt <- unique(as.integer(nxt))
      nxt <- nxt[!is.na(nxt)]
      if (length(nxt) == 0) break
      
      current <- nxt
      if (all(hier$depth[hier$id %in% current] >= want_depth, na.rm=TRUE)) break
    }
    
    unique(hier$id[hier$depth == want_depth & hier$id %in% current])
  }
  
  # ---- Determine top depth_start parents from AveragedRegions.xlsx ----
  avg_df <- readxl::read_excel(averaged_file, sheet = "Averages")
  if ("Region" %in% names(avg_df)) {
    avg_df <- avg_df |> dplyr::rename(region = Region, mean_density = MeanDensity)
  }
  top_parents <- avg_df |>
    dplyr::arrange(dplyr::desc(mean_density)) |>
    dplyr::slice_head(n = top_n_parents) |>
    dplyr::pull(region)
  
  # ---- For left/right, map top parents to IDs using a filtered file ----
  map_files <- list.files(filtered_folder_for_id_map, pattern="\\.xlsx$", full.names=TRUE)
  if (length(map_files) == 0) stop("No .xlsx files found in filtered_folder_for_id_map")
  map_df <- readxl::read_excel(map_files[1])
  region_col <- if ("region_base" %in% names(map_df)) "region_base" else if ("name" %in% names(map_df)) "name" else NA_character_
  if (is.na(region_col)) stop("Filtered file needs 'region_base' or 'name'.")
  if (!all(c("id","depth") %in% names(map_df))) stop("Filtered-by-depth files must include id + depth.")
  
  parent_id_map <- map_df |>
    dplyr::transmute(region = as.character(.data[[region_col]]),
                     id = as.integer(id),
                     depth = as.integer(depth)) |>
    dplyr::filter(depth == depth_start, region %in% top_parents)
  
  # ---- Read all CSVs and create per-sample pools at depth 6, 5, and depth_start ----
  csv_files <- list.files(original_csv_folder, pattern="\\.csv$", full.names=TRUE)
  if (length(csv_files) == 0) stop("No .csv files found in original_csv_folder")
  
  long_by_depth <- function(depth_want) {
    out <- list()
    for (csv in csv_files) {
      sample_id <- tools::file_path_sans_ext(basename(csv))
      hier <- build_hierarchy(csv)
      
      parent_ids <- if (laterality %in% c("left","right")) {
        parent_id_map$id
      } else {
        hier |> dplyr::filter(depth == depth_start, region_base %in% top_parents) |> dplyr::pull(id)
      }
      
      ids <- unique(unlist(lapply(parent_ids, function(pid) {
        descendants_at_depth(hier, pid, depth_want)
      }), use.names = FALSE))
      
      if (length(ids) == 0) next
      
      out[[length(out)+1]] <- hier |>
        dplyr::filter(id %in% ids, depth == depth_want) |>
        dplyr::transmute(SampleID = sample_id,
                         Region = region_base,
                         Acronym = acronym_base,
                         Depth = depth,
                         Density = density,
                         parent_id = parent_structure_id)
    }
    dplyr::bind_rows(out)
  }
  
  long6 <- long_by_depth(target_depth)
  long5 <- long_by_depth(target_depth - 1)
  long4 <- long_by_depth(depth_start) # last-resort fallback depth
  
  if (nrow(long6) == 0) stop("No depth-", target_depth, " candidates found under top parents.")
  
  # Parent map for depth6 parent resolution (depth 5)
  parent_lookup_list <- list()
  for (csv in csv_files) {
    hier <- build_hierarchy(csv)
    parent_lookup_list[[length(parent_lookup_list)+1]] <- hier |>
      dplyr::filter(depth == target_depth - 1) |>
      dplyr::transmute(parent_id = id,
                       ParentRegion = region_base,
                       ParentAcronym = acronym_base,
                       ParentDepth = depth)
  }
  parent_map <- dplyr::bind_rows(parent_lookup_list) |>
    dplyr::filter(!is.na(parent_id)) |>
    dplyr::group_by(parent_id, ParentRegion, ParentAcronym, ParentDepth) |>
    dplyr::summarise(n = dplyr::n(), .groups="drop") |>
    dplyr::arrange(parent_id, dplyr::desc(n)) |>
    dplyr::group_by(parent_id) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup()
  
  long6 <- long6 |> dplyr::left_join(parent_map, by = "parent_id")
  
  # Rank tables
  rank_table <- function(long_df, include_parent_cols = FALSE) {
    if (nrow(long_df) == 0) {
      return(dplyr::tibble(Region=character(), Acronym=character(), Depth=integer(), MeanDensity=numeric()))
    }
    if (include_parent_cols) {
      long_df |>
        dplyr::group_by(Region) |>
        dplyr::summarise(
          MeanDensity = mean(Density, na.rm=TRUE),
          Acronym = dplyr::first(na.omit(Acronym)),
          Depth = dplyr::first(Depth),
          ParentRegion = dplyr::first(na.omit(ParentRegion)),
          ParentAcronym = dplyr::first(na.omit(ParentAcronym)),
          ParentDepth = dplyr::first(na.omit(ParentDepth)),
          .groups="drop"
        ) |>
        dplyr::arrange(dplyr::desc(MeanDensity))
    } else {
      long_df |>
        dplyr::group_by(Region) |>
        dplyr::summarise(
          MeanDensity = mean(Density, na.rm=TRUE),
          Acronym = dplyr::first(na.omit(Acronym)),
          Depth = dplyr::first(Depth),
          .groups="drop"
        ) |>
        dplyr::arrange(dplyr::desc(MeanDensity))
    }
  }
  
  means6 <- rank_table(long6, include_parent_cols = TRUE)
  means5 <- rank_table(long5, include_parent_cols = FALSE)
  means4 <- rank_table(long4, include_parent_cols = FALSE)
  
  # Iterative collapse + refill
  collapsed_groups_all <- dplyr::tibble()
  collapsed_parents_all <- character()
  
  # Start with top30 at depth 6
  current <- means6 |> dplyr::slice_head(n = top_n_final) |>
    dplyr::select(Region, Acronym, Depth, MeanDensity, ParentRegion, ParentAcronym, ParentDepth)
  
  iter <- 1
  repeat {
    # Count children by parent among depth6 entries only
    over <- current |>
      dplyr::filter(Depth == target_depth, !is.na(ParentRegion)) |>
      dplyr::count(ParentRegion, name="n_children") |>
      dplyr::filter(n_children > collapse_threshold)
    
    if (nrow(over) == 0) break
    if (iter > max_iter) {
      warning("Reached max_iter=", max_iter, " before eliminating overrepresented parents.")
      break
    }
    
    collapse_parents <- over$ParentRegion
    collapsed_parents_all <- unique(c(collapsed_parents_all, collapse_parents))
    
    # Record which children are being collapsed this iteration (from current list)
    iter_groups <- current |>
      dplyr::filter(Depth == target_depth, ParentRegion %in% collapse_parents) |>
      dplyr::select(ParentRegion, ParentAcronym, ParentDepth,
                    ChildRegion = Region, ChildAcronym = Acronym, ChildDepth = Depth,
                    ChildMeanDensity = MeanDensity) |>
      dplyr::mutate(Iteration = iter) |>
      dplyr::arrange(ParentRegion, dplyr::desc(ChildMeanDensity))
    
    collapsed_groups_all <- dplyr::bind_rows(collapsed_groups_all, iter_groups)
    
    # Remove those children
    collapsed_children <- current |>
      dplyr::filter(Depth == target_depth, ParentRegion %in% collapse_parents) |>
      dplyr::pull(Region)
    
    kept <- current |>
      dplyr::filter(!(Region %in% collapsed_children)) |>
      dplyr::select(Region, Acronym, Depth, MeanDensity, ParentRegion, ParentAcronym, ParentDepth)
    
    # Add parents (depth 5) as replacements
    parent_rows <- means5 |>
      dplyr::filter(Region %in% collapse_parents) |>
      dplyr::transmute(Region, Acronym, Depth, MeanDensity,
                       ParentRegion = NA_character_,
                       ParentAcronym = NA_character_,
                       ParentDepth = NA_integer_)
    
    current <- dplyr::bind_rows(kept, parent_rows) |>
      dplyr::arrange(dplyr::desc(MeanDensity))
    
    # Exclude all depth6 children of collapsed parents from *any* refill
    excluded6 <- means6 |>
      dplyr::filter(ParentRegion %in% collapsed_parents_all) |>
      dplyr::pull(Region)
    
    # Helper: refill from ranked table until 30
    refill_from <- function(cur, ranked, exclude_regions) {
      need <- top_n_final - nrow(cur)
      if (need <= 0) return(cur |> dplyr::slice_head(n = top_n_final))
      
      already <- cur$Region
      add <- ranked |>
        dplyr::filter(!(Region %in% already),
                      !(Region %in% exclude_regions)) |>
        dplyr::slice_head(n = need)
      
      dplyr::bind_rows(cur, add) |>
        dplyr::arrange(dplyr::desc(MeanDensity)) |>
        dplyr::slice_head(n = top_n_final)
    }
    
    # Refill priority: depth6 -> depth5 -> depth4 (depth_start as last resort)
    current <- refill_from(current, means6 |>
                             dplyr::select(Region, Acronym, Depth, MeanDensity,
                                           ParentRegion, ParentAcronym, ParentDepth),
                           exclude_regions = excluded6)
    
    if (nrow(current) < top_n_final) {
      current <- refill_from(current, means5 |>
                               dplyr::select(Region, Acronym, Depth, MeanDensity) |>
                               dplyr::mutate(ParentRegion=NA_character_, ParentAcronym=NA_character_, ParentDepth=NA_integer_),
                             exclude_regions = character(0))
    }
    
    if (nrow(current) < top_n_final) {
      current <- refill_from(current, means4 |>
                               dplyr::select(Region, Acronym, Depth, MeanDensity) |>
                               dplyr::mutate(ParentRegion=NA_character_, ParentAcronym=NA_character_, ParentDepth=NA_integer_),
                             exclude_regions = character(0))
    }
    
    iter <- iter + 1
  }
  
  # Final cleaned list (drop parent columns from Averages output)
  cleaned_list <- current |>
    dplyr::select(Region, Acronym, Depth, MeanDensity)
  
  # Build final per-sample densities wide table using whichever depth the region ended up being
  final_regions <- cleaned_list$Region
  
  long_final <- dplyr::bind_rows(
    long6 |> dplyr::select(SampleID, Region, Density),
    long5 |> dplyr::select(SampleID, Region, Density),
    long4 |> dplyr::select(SampleID, Region, Density)
  ) |>
    dplyr::filter(Region %in% final_regions) |>
    dplyr::group_by(SampleID, Region) |>
    dplyr::summarise(Density = mean(Density, na.rm=TRUE), .groups="drop")
  
  wide <- long_final |>
    tidyr::pivot_wider(names_from = Region, values_from = Density)
  
  wide <- wide[, c("SampleID", final_regions), drop=FALSE]
  
  # If CollapsedGroups ended up empty, keep a header-only tibble
  if (nrow(collapsed_groups_all) == 0) {
    collapsed_groups_all <- dplyr::tibble(
      ParentRegion=character(), ParentAcronym=character(), ParentDepth=integer(),
      ChildRegion=character(), ChildAcronym=character(), ChildDepth=integer(),
      ChildMeanDensity=numeric(), Iteration=integer()
    )
  }
  
  # ---- Build AveragesAllRegions (full ranked pools under top parents) ----
  # Mark excluded depth6 children for transparency (based on FINAL collapsed parents)
  excluded6_final <- character(0)
  if (length(collapsed_parents_all) > 0) {
    excluded6_final <- means6 |>
      dplyr::filter(ParentRegion %in% collapsed_parents_all) |>
      dplyr::pull(Region)
  }
  
  avg_all_6 <- means6 |>
    dplyr::transmute(
      Region, Acronym, Depth, MeanDensity,
      ParentRegion, ParentAcronym, ParentDepth,
      Pool = paste0("Depth", target_depth),
      ExcludedBecauseCollapsed = Region %in% excluded6_final
    )
  
  avg_all_5 <- means5 |>
    dplyr::transmute(
      Region, Acronym, Depth, MeanDensity,
      ParentRegion = NA_character_,
      ParentAcronym = NA_character_,
      ParentDepth = NA_integer_,
      Pool = paste0("Depth", target_depth - 1),
      ExcludedBecauseCollapsed = FALSE
    )
  
  avg_all_4 <- means4 |>
    dplyr::transmute(
      Region, Acronym, Depth, MeanDensity,
      ParentRegion = NA_character_,
      ParentAcronym = NA_character_,
      ParentDepth = NA_integer_,
      Pool = paste0("Depth", depth_start, "_LastResort"),
      ExcludedBecauseCollapsed = FALSE
    )
  
  averages_all_regions <- dplyr::bind_rows(avg_all_6, avg_all_5, avg_all_4) |>
    dplyr::arrange(dplyr::desc(MeanDensity))

  # ---- Build RankedRegions (all candidate regions from AveragesAllRegions) ----
  # Use all candidate subject-by-region densities from the same pools that feed
  # AveragesAllRegions, then compute mean/median density per region and dense ranks
  # within each subject (highest density = rank 1; missing values remain NA).
  ranked_long <- dplyr::bind_rows(
    long6 |> dplyr::select(SampleID, Region, Density),
    long5 |> dplyr::select(SampleID, Region, Density),
    long4 |> dplyr::select(SampleID, Region, Density)
  ) |>
    dplyr::group_by(SampleID, Region) |>
    dplyr::summarise(Density = mean(Density, na.rm = TRUE), .groups = "drop")

  ranked_summary <- ranked_long |>
    dplyr::group_by(Region) |>
    dplyr::summarise(
      MeanDensity = mean(Density, na.rm = TRUE),
      MedianDensity = stats::median(Density, na.rm = TRUE),
      .groups = "drop"
    )

  ranked_subjects <- ranked_long |>
    dplyr::group_by(SampleID) |>
    dplyr::mutate(Rank = dplyr::dense_rank(dplyr::desc(Density))) |>
    dplyr::ungroup() |>
    dplyr::select(Region, SampleID, Rank) |>
    tidyr::pivot_wider(names_from = SampleID, values_from = Rank)

  ranked_regions <- ranked_summary |>
    dplyr::left_join(ranked_subjects, by = "Region") |>
    dplyr::arrange(dplyr::desc(MeanDensity), dplyr::desc(MedianDensity), Region)
  
  writexl::write_xlsx(
    list(
      "Densities" = wide,
      "Averages"  = cleaned_list,
      "AveragesAllRegions" = averages_all_regions,
      "RankedRegions" = ranked_regions,
      "CollapsedGroups" = collapsed_groups_all
    ),
    path = output_file
  )
  
  invisible(list(densities = wide,
                 averages = cleaned_list,
                 collapsed_groups = collapsed_groups_all,
                 collapsed_parents = collapsed_parents_all))
}

# -----------------------------------------------------------------------------
# Function 4: Top30FromRegionCleaned (focused parent-region analysis)
# -----------------------------------------------------------------------------
#Focused filtering for whole-brain data analysis

Top30FromRegionCleaned <- function(original_csv_folder,
                                   start_region,
                                   output_file = NULL,
                                   start_depth = NULL,          # optional: one value for all start regions, or named/vector values matching start_region
                                   target_depth = 6,
                                   top_n_final = 30,
                                   collapse_threshold = 3,      # > 3 depth-6 children from same depth-5 parent triggers collapse
                                   laterality = c("combined", "left", "right"),
                                   density_col = "density (cells/mm^3)",
                                   max_iter = 10,
                                   case_insensitive_match = TRUE) {
  
  laterality <- match.arg(laterality)
  
  if (!requireNamespace("writexl", quietly = TRUE)) stop("Install 'writexl'")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Install 'dplyr'")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Install 'tidyr'")
  if (!requireNamespace("stringr", quietly = TRUE)) stop("Install 'stringr'")
  
  library(dplyr)
  library(tidyr)
  library(stringr)
  
  # ---------------- Helpers ----------------
  strip_lr <- function(x) {
    x <- stringr::str_replace(x, stringr::regex("^left\\s+",  ignore_case = TRUE), "")
    x <- stringr::str_replace(x, stringr::regex("^right\\s+", ignore_case = TRUE), "")
    x
  }
  strip_acr_lr <- function(x) stringr::str_replace(x, "(-L|-R)$", "")
  
  # Build hierarchy for one CSV, with combined-mode merge support
  build_hierarchy <- function(csv_file) {
    df <- read.csv(csv_file, stringsAsFactors = FALSE, check.names = FALSE)
    
    req <- c("id","name","acronym","parent_structure_id","depth","count","volume (mm^3)", density_col)
    miss <- setdiff(req, names(df))
    if (length(miss) > 0) stop("Missing columns in ", basename(csv_file), ": ", paste(miss, collapse=", "))
    
    df <- df[!is.na(df$depth) & df$depth >= 0, , drop = FALSE]
    
    # hemisphere from ORIGINAL fields
    hemi_from_name <- rep(NA_character_, nrow(df))
    hemi_from_name[grepl("^left\\s+",  df$name, ignore.case = TRUE)] <- "left"
    hemi_from_name[grepl("^right\\s+", df$name, ignore.case = TRUE)] <- "right"
    
    hemi_from_acr <- rep(NA_character_, nrow(df))
    hemi_from_acr[grepl("-L$", df$acronym)] <- "left"
    hemi_from_acr[grepl("-R$", df$acronym)] <- "right"
    
    df$hemisphere <- dplyr::coalesce(hemi_from_name, hemi_from_acr, "combined")
    
    # base fields
    df$region_base  <- strip_lr(df$name)
    df$acronym_base <- strip_acr_lr(df$acronym)
    
    # parent_base lookup (to rebuild combined tree)
    id_map <- df[, c("id","region_base","depth")]
    names(id_map) <- c("id","parent_region_base_lookup","parent_depth_lookup")
    df$parent_structure_id_chr <- as.character(df$parent_structure_id)
    id_map$id_chr <- as.character(id_map$id)
    
    df <- dplyr::left_join(
      df,
      id_map[, c("id_chr","parent_region_base_lookup","parent_depth_lookup")],
      by = c("parent_structure_id_chr" = "id_chr")
    )
    df$parent_base <- df$parent_region_base_lookup
    
    # Left/right mode: keep only that hemisphere, keep original IDs
    if (laterality %in% c("left","right")) {
      out <- df[df$hemisphere == laterality, , drop=FALSE] |>
        dplyr::transmute(
          id = as.integer(id),
          parent_structure_id = as.integer(parent_structure_id),
          depth = as.integer(depth),
          region_base = as.character(region_base),
          acronym_base = as.character(acronym_base),
          density = as.numeric(.data[[density_col]])
        )
      return(out)
    }
    
    # Combined mode: prefer L/R rows if present; else use combined rows
    df_lr <- df[df$hemisphere %in% c("left","right"), , drop=FALSE]
    df_cb <- df[df$hemisphere == "combined", , drop=FALSE]
    src <- if (nrow(df_lr) > 0) df_lr else df_cb
    
    src$count <- as.numeric(src$count)
    src[["volume (mm^3)"]] <- as.numeric(src[["volume (mm^3)"]])
    
    merged <- src |>
      dplyr::mutate(parent_base = dplyr::coalesce(parent_base, NA_character_)) |>
      dplyr::group_by(depth, parent_base, region_base, acronym_base) |>
      dplyr::summarise(
        count = sum(count, na.rm = TRUE),
        volume = sum(`volume (mm^3)`, na.rm = TRUE),
        .groups = "drop"
      ) |>
      dplyr::mutate(density = dplyr::if_else(volume > 0, count / volume, NA_real_)) |>
      dplyr::arrange(depth, parent_base, region_base) |>
      dplyr::mutate(id = dplyr::row_number())
    
    parent_lookup <- merged |>
      dplyr::select(id, depth, region_base) |>
      dplyr::rename(parent_id = id,
                    parent_depth = depth,
                    parent_region_base = region_base)
    
    merged <- merged |>
      dplyr::left_join(parent_lookup, by = c("parent_base" = "parent_region_base")) |>
      dplyr::mutate(
        parent_structure_id = dplyr::if_else(!is.na(parent_base) & parent_depth == depth - 1,
                                             parent_id, NA_integer_)
      ) |>
      dplyr::select(id, parent_structure_id, depth, region_base, acronym_base, density)
    
    merged
  }
  
  # Descendants at a specific depth under a set of starting IDs
  descendants_at_depth <- function(hier, start_ids, want_depth) {
    start_ids <- unique(as.integer(start_ids))
    start_ids <- start_ids[!is.na(start_ids)]
    if (length(start_ids) == 0) return(integer(0))
    
    children_map <- split(hier$id, hier$parent_structure_id)
    current <- start_ids
    
    repeat {
      depths_now <- unique(hier$depth[hier$id %in% current])
      if (any(depths_now == want_depth, na.rm = TRUE)) break
      
      nxt <- unlist(children_map[as.character(current)], use.names = FALSE)
      nxt <- unique(as.integer(nxt))
      nxt <- nxt[!is.na(nxt)]
      if (length(nxt) == 0) break
      
      current <- nxt
      if (all(hier$depth[hier$id %in% current] >= want_depth, na.rm = TRUE)) break
    }
    
    unique(hier$id[hier$depth == want_depth & hier$id %in% current])
  }
  
  # Ranking helper
  rank_table <- function(long_df, include_parent_cols = FALSE) {
    if (nrow(long_df) == 0) {
      if (include_parent_cols) {
        return(tibble(Region=character(), Acronym=character(), Depth=integer(), MeanDensity=numeric(),
                      ParentRegion=character(), ParentAcronym=character(), ParentDepth=integer()))
      }
      return(tibble(Region=character(), Acronym=character(), Depth=integer(), MeanDensity=numeric()))
    }
    if (include_parent_cols) {
      long_df %>%
        group_by(Region) %>%
        summarise(
          MeanDensity = mean(Density, na.rm = TRUE),
          Acronym = first(na.omit(Acronym)),
          Depth = first(Depth),
          ParentRegion = first(na.omit(ParentRegion)),
          ParentAcronym = first(na.omit(ParentAcronym)),
          ParentDepth = first(na.omit(ParentDepth)),
          .groups = "drop"
        ) %>%
        arrange(desc(MeanDensity))
    } else {
      long_df %>%
        group_by(Region) %>%
        summarise(
          MeanDensity = mean(Density, na.rm = TRUE),
          Acronym = first(na.omit(Acronym)),
          Depth = first(Depth),
          .groups = "drop"
        ) %>%
        arrange(desc(MeanDensity))
    }
  }
  
  # ---------------- Main ----------------
  csv_files <- list.files(original_csv_folder, pattern = "\\.csv$", full.names = TRUE)
  if (length(csv_files) == 0) stop("No .csv files found in: ", original_csv_folder)
  
  # Pools
  long6_list <- list()
  long5_list <- list()
  longS_list <- list()  # last-resort: start_depth itself
  parent_lookup_list <- list()
  
  # Track each sample’s resolved start depth (if start_depth not provided)
  resolved_start_depths <- integer(0)
  
  # Match function for region names
  match_region <- function(x) {
    if (case_insensitive_match) {
      tolower(x) == tolower(start_region)
    } else {
      x == start_region
    }
  }
  
  for (csv in csv_files) {
    sample_id <- tools::file_path_sans_ext(basename(csv))
    hier <- build_hierarchy(csv)
    
    # Find the starting node(s) by region name (and optionally depth)
    candidates <- hier %>% filter(match_region(region_base))
    if (nrow(candidates) == 0) next
    
    if (!is.null(start_depth)) {
      candidates <- candidates %>% filter(depth == start_depth)
    } else {
      # Choose the broadest (minimum depth) match as the start node for that sample
      min_d <- min(candidates$depth, na.rm = TRUE)
      candidates <- candidates %>% filter(depth == min_d)
      resolved_start_depths <- c(resolved_start_depths, min_d)
    }
    
    start_ids <- candidates$id
    if (length(start_ids) == 0) next
    
    # Determine this sample's start depth for last-resort backfill
    this_start_depth <- if (!is.null(start_depth)) start_depth else min(candidates$depth, na.rm = TRUE)
    
    # Depth 6 descendants
    ids6 <- descendants_at_depth(hier, start_ids, target_depth)
    if (length(ids6) > 0) {
      long6_list[[length(long6_list)+1]] <- hier %>%
        filter(id %in% ids6, depth == target_depth) %>%
        transmute(
          SampleID = sample_id,
          Region = region_base,
          Acronym = acronym_base,
          Depth = depth,
          Density = density,
          parent_id = parent_structure_id
        )
    }
    
    # Depth 5 descendants
    ids5 <- descendants_at_depth(hier, start_ids, target_depth - 1)
    if (length(ids5) > 0) {
      long5_list[[length(long5_list)+1]] <- hier %>%
        filter(id %in% ids5, depth == (target_depth - 1)) %>%
        transmute(
          SampleID = sample_id,
          Region = region_base,
          Acronym = acronym_base,
          Depth = depth,
          Density = density
        )
    }
    
    # Last-resort pool: nodes at the start depth within this subtree
    idsS <- descendants_at_depth(hier, start_ids, this_start_depth)
    if (length(idsS) > 0) {
      longS_list[[length(longS_list)+1]] <- hier %>%
        filter(id %in% idsS, depth == this_start_depth) %>%
        transmute(
          SampleID = sample_id,
          Region = region_base,
          Acronym = acronym_base,
          Depth = depth,
          Density = density
        )
    }
    
    # Parent map for depth-6 parent resolution (depth 5)
    parent_lookup_list[[length(parent_lookup_list)+1]] <- hier %>%
      filter(depth == target_depth - 1) %>%
      transmute(parent_id = id,
                ParentRegion = region_base,
                ParentAcronym = acronym_base,
                ParentDepth = depth)
  }
  
  long6 <- bind_rows(long6_list)
  long5 <- bind_rows(long5_list)
  longS <- bind_rows(longS_list)
  
  if (nrow(long6) == 0) {
    stop("No depth-", target_depth, " descendants found for start_region='", start_region, "'. ",
         "Check spelling/case, laterality, and that the region exists in these CSVs.")
  }
  
  # Resolve parent mapping (most common)
  parent_map <- bind_rows(parent_lookup_list) %>%
    filter(!is.na(parent_id)) %>%
    group_by(parent_id, ParentRegion, ParentAcronym, ParentDepth) %>%
    summarise(n = n(), .groups = "drop") %>%
    arrange(parent_id, desc(n)) %>%
    group_by(parent_id) %>%
    slice_head(n = 1) %>%
    ungroup()
  
  long6 <- long6 %>% left_join(parent_map, by = "parent_id")
  
  # Ranked pools
  means6 <- rank_table(long6, include_parent_cols = TRUE)
  means5 <- rank_table(long5, include_parent_cols = FALSE)
  meansS <- rank_table(longS, include_parent_cols = FALSE)
  
  # ---------------- Iterative collapse + refill ----------------
  collapsed_groups_all <- tibble()
  collapsed_parents_all <- character()
  
  current <- means6 %>%
    slice_head(n = top_n_final) %>%
    select(Region, Acronym, Depth, MeanDensity, ParentRegion, ParentAcronym, ParentDepth)
  
  refill_from <- function(cur, ranked, exclude_regions) {
    need <- top_n_final - nrow(cur)
    if (need <= 0) return(cur %>% slice_head(n = top_n_final))
    already <- cur$Region
    add <- ranked %>%
      filter(!(Region %in% already),
             !(Region %in% exclude_regions)) %>%
      slice_head(n = need)
    bind_rows(cur, add) %>%
      arrange(desc(MeanDensity)) %>%
      slice_head(n = top_n_final)
  }
  
  iter <- 1
  repeat {
    # Overrepresented parents among depth-6 entries currently in list
    over <- current %>%
      filter(Depth == target_depth, !is.na(ParentRegion)) %>%
      count(ParentRegion, name = "n_children") %>%
      filter(n_children > collapse_threshold)
    
    if (nrow(over) == 0) break
    if (iter > max_iter) {
      warning("Reached max_iter=", max_iter, " before eliminating overrepresented parents.")
      break
    }
    
    collapse_parents <- over$ParentRegion
    collapsed_parents_all <- unique(c(collapsed_parents_all, collapse_parents))
    
    # Record which children are being collapsed this iteration
    iter_groups <- current %>%
      filter(Depth == target_depth, ParentRegion %in% collapse_parents) %>%
      select(ParentRegion, ParentAcronym, ParentDepth,
             ChildRegion = Region, ChildAcronym = Acronym, ChildDepth = Depth,
             ChildMeanDensity = MeanDensity) %>%
      mutate(Iteration = iter) %>%
      arrange(ParentRegion, desc(ChildMeanDensity))
    
    collapsed_groups_all <- bind_rows(collapsed_groups_all, iter_groups)
    
    # Remove collapsed children
    collapsed_children <- iter_groups$ChildRegion
    kept <- current %>%
      filter(!(Region %in% collapsed_children)) %>%
      select(Region, Acronym, Depth, MeanDensity, ParentRegion, ParentAcronym, ParentDepth)
    
    # Add collapsed parents (depth-5 rows) as replacements
    parent_rows <- means5 %>%
      filter(Region %in% collapse_parents) %>%
      transmute(Region, Acronym, Depth, MeanDensity,
                ParentRegion = NA_character_, ParentAcronym = NA_character_, ParentDepth = NA_integer_)
    
    current <- bind_rows(kept, parent_rows) %>%
      arrange(desc(MeanDensity))
    
    # Exclude all depth-6 children of any collapsed parent (across all iterations)
    excluded6 <- means6 %>%
      filter(ParentRegion %in% collapsed_parents_all) %>%
      pull(Region)
    
    # Refill priority: depth6 -> depth5 -> start-depth (last resort)
    current <- refill_from(current,
                           means6 %>% select(Region, Acronym, Depth, MeanDensity, ParentRegion, ParentAcronym, ParentDepth),
                           exclude_regions = excluded6)
    
    if (nrow(current) < top_n_final) {
      current <- refill_from(current,
                             means5 %>% transmute(Region, Acronym, Depth, MeanDensity,
                                                  ParentRegion = NA_character_, ParentAcronym = NA_character_, ParentDepth = NA_integer_),
                             exclude_regions = character(0))
    }
    
    if (nrow(current) < top_n_final) {
      current <- refill_from(current,
                             meansS %>% transmute(Region, Acronym, Depth, MeanDensity,
                                                  ParentRegion = NA_character_, ParentAcronym = NA_character_, ParentDepth = NA_integer_),
                             exclude_regions = character(0))
    }
    
    iter <- iter + 1
  }
  
  # Final cleaned list
  cleaned_list <- current %>% select(Region, Acronym, Depth, MeanDensity)
  
  # ---------------- AveragesAllRegions sheet ----------------
  excluded6_final <- if (length(collapsed_parents_all) > 0) {
    means6 %>% filter(ParentRegion %in% collapsed_parents_all) %>% pull(Region)
  } else character()
  
  avg_all_6 <- means6 %>%
    transmute(Region, Acronym, Depth, MeanDensity,
              ParentRegion, ParentAcronym, ParentDepth,
              Pool = paste0("Depth", target_depth),
              ExcludedBecauseCollapsed = Region %in% excluded6_final)
  
  avg_all_5 <- means5 %>%
    transmute(Region, Acronym, Depth, MeanDensity,
              ParentRegion = NA_character_, ParentAcronym = NA_character_, ParentDepth = NA_integer_,
              Pool = paste0("Depth", target_depth - 1),
              ExcludedBecauseCollapsed = FALSE)
  
  # start-depth pool (last resort) — may vary by sample if start_depth not provided, but the ranked list is still meaningful
  avg_all_S <- meansS %>%
    transmute(Region, Acronym, Depth, MeanDensity,
              ParentRegion = NA_character_, ParentAcronym = NA_character_, ParentDepth = NA_integer_,
              Pool = "StartDepth_LastResort",
              ExcludedBecauseCollapsed = FALSE)
  
  averages_all_regions <- bind_rows(avg_all_6, avg_all_5, avg_all_S) %>%
    arrange(desc(MeanDensity))
  
  # ---------------- Densities sheet (per-sample values) ----------------
  final_regions <- cleaned_list$Region
  
  long_final <- bind_rows(
    long6 %>% select(SampleID, Region, Density),
    long5 %>% select(SampleID, Region, Density),
    longS %>% select(SampleID, Region, Density)
  ) %>%
    filter(Region %in% final_regions) %>%
    group_by(SampleID, Region) %>%
    summarise(Density = mean(Density, na.rm = TRUE), .groups = "drop")
  
  wide <- long_final %>%
    tidyr::pivot_wider(names_from = Region, values_from = Density)
  
  wide <- wide[, c("SampleID", final_regions), drop = FALSE]
  
  # CollapsedGroups sheet (ensure it has headers even if empty)
  if (nrow(collapsed_groups_all) == 0) {
    collapsed_groups_all <- tibble(
      ParentRegion=character(), ParentAcronym=character(), ParentDepth=integer(),
      ChildRegion=character(), ChildAcronym=character(), ChildDepth=integer(),
      ChildMeanDensity=numeric(), Iteration=integer()
    )
  }
  
  writexl::write_xlsx(
    list(
      "Densities" = wide,
      "Averages" = cleaned_list,
      "AveragesAllRegions" = averages_all_regions,
      "CollapsedGroups" = collapsed_groups_all
    ),
    path = output_file
  )
  
  invisible(list(
    densities = wide,
    averages = cleaned_list,
    averages_all_regions = averages_all_regions,
    collapsed_groups = collapsed_groups_all,
    collapsed_parents = collapsed_parents_all
  ))
}

# Example usage:
# Top30FromRegionCleaned(
#   original_csv_folder = "~/Library/CloudStorage/Box-Box/Centanni/Ava/Sparta Collab/Co-expression/Unfiltered/Impov",
#   start_region = c("Pallidum", "Hypothalamus"),  # must match region_base naming
#   start_depth = 4,                                # one value for all start regions
#   laterality = "combined",
#   output_file = "~/Library/CloudStorage/Box-Box/Centanni/Ava/Sparta Collab/Focused analysis/Impov/CoexpressionImpovFocused_Top30SpartaCleaned.xlsx",
#   collapse_threshold = 30
# )
#
# You can also provide region-specific depths:
# start_depth = c(Pallidum = 4, Hypothalamus = 3)

# -----------------------------------------------------------------------------
# Function 5: make_density_chord_parentgroups (group comparison chord diagram)
# -----------------------------------------------------------------------------
## Chord diagram analysis to compare region overlap between groups
## Cody by: Ava Shipman, Centanni Lab
## Cleaned up + legend-safe layout (no legend overlap with circos)

make_density_chord_parentgroups <- function(group_files,
                                            group_names = NULL,
                                            sheet_avg = "Averages",
                                            sheet_all = "AveragesAllRegions",
                                            region_col = "Region",
                                            value_col = "MeanDensity",
                                            parent_col = "ParentRegion",
                                            acronym_col = "Acronym",
                                            # Display:
                                            top_regions = 30,
                                            min_density = NULL,
                                            use_acronym_labels = TRUE,
                                            show_parent_labels = TRUE,
                                            # Styling:
                                            shared_fill = "red",
                                            unique_fill = "blue",
                                            group_fill = "grey25",
                                            shared_link = shared_fill,
                                            unique_link = unique_fill,
                                            link_alpha = 0.35,
                                            # Tracks / labels:
                                            label_track_height = 0.12,
                                            parent_track_height = 0.08,
                                            group_label_cex = 0.9,
                                            region_label_cex = 0.55,
                                            parent_label_cex = 0.60,
                                            parent_label_facing = c("radial", "bending.inside"),
                                            # Legends (drawn in a separate right-hand panel):
                                            show_legends = TRUE,
                                            legend_panel_width = 2.2,
                                            density_legend_n = 3,       # number of example bars for the density scale
                                            start_degree = 90) {
  
  parent_label_facing <- match.arg(parent_label_facing)
  
  if (!requireNamespace("readxl", quietly = TRUE)) stop("Install readxl")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Install dplyr")
  if (!requireNamespace("circlize", quietly = TRUE)) stop("Install circlize")
  
  library(readxl)
  library(dplyr)
  library(circlize)
  
  if (length(group_files) < 2) stop("Provide at least 2 group_files.")
  if (is.null(group_names)) group_names <- tools::file_path_sans_ext(basename(group_files))
  if (length(group_names) != length(group_files)) stop("group_names must match group_files length.")
  
  # ---------- Read Averages (for chord links) ----------
  read_group_avg <- function(path, group_name) {
    df <- readxl::read_excel(path, sheet = sheet_avg)
    if (!(region_col %in% names(df))) stop("Missing ", region_col, " in ", basename(path))
    if (!(value_col %in% names(df))) stop("Missing ", value_col, " in ", basename(path))
    
    has_acr <- use_acronym_labels && (acronym_col %in% names(df))
    
    df %>%
      transmute(
        Region = as.character(.data[[region_col]]),
        Label  = if (has_acr) as.character(.data[[acronym_col]]) else as.character(.data[[region_col]]),
        Group  = group_name,
        Density = as.numeric(.data[[value_col]])
      ) %>%
      filter(!is.na(Region), Region != "", !is.na(Density)) %>%
      group_by(Group, Region, Label) %>%
      summarise(Density = mean(Density, na.rm = TRUE), .groups = "drop")
  }
  
  links_long <- bind_rows(Map(read_group_avg, group_files, group_names))
  if (!is.null(min_density)) links_long <- links_long %>% filter(Density >= min_density)
  
  # ---------- Parent map from AveragesAllRegions ----------
  read_parent_map <- function(path) {
    if (!(sheet_all %in% readxl::excel_sheets(path))) return(NULL)
    df <- readxl::read_excel(path, sheet = sheet_all)
    if (!(region_col %in% names(df))) return(NULL)
    if (!(parent_col %in% names(df))) return(NULL)
    
    df %>%
      transmute(
        Region = as.character(.data[[region_col]]),
        ParentRegion = as.character(.data[[parent_col]])
      ) %>%
      filter(!is.na(Region), Region != "",
             !is.na(ParentRegion), ParentRegion != "") %>%
      distinct()
  }
  
  parent_raw <- bind_rows(lapply(group_files, read_parent_map))
  if (!is.null(parent_raw) && nrow(parent_raw) > 0) {
    parent_map <- parent_raw %>%
      count(Region, ParentRegion, name = "n") %>%
      arrange(Region, desc(n)) %>%
      group_by(Region) %>%
      slice_head(n = 1) %>%
      ungroup() %>%
      select(Region, ParentRegion)
  } else {
    parent_map <- tibble(Region = unique(links_long$Region), ParentRegion = "Unmapped")
  }
  
  # ---------- Shared vs unique ----------
  presence <- links_long %>%
    distinct(Region, Group) %>%
    count(Region, name = "n_groups") %>%
    mutate(category = ifelse(n_groups == length(group_names), "shared", "unique"))
  
  # ---------- Totals for ordering + legend range ----------
  region_totals <- links_long %>%
    group_by(Region) %>%
    summarise(total = sum(Density, na.rm = TRUE), .groups = "drop") %>%
    left_join(parent_map, by = "Region") %>%
    left_join(presence %>% select(Region, category), by = "Region") %>%
    mutate(
      ParentRegion = ifelse(is.na(ParentRegion), "Unmapped", ParentRegion),
      category = ifelse(is.na(category), "unique", category)
    ) %>%
    arrange(desc(total))
  
  if (!is.null(top_regions)) {
    keep <- head(region_totals$Region, top_regions)
    links_long <- links_long %>% filter(Region %in% keep)
    presence <- presence %>% filter(Region %in% keep)
    region_totals <- region_totals %>% filter(Region %in% keep)
  }
  
  # order by ParentRegion, then density within parent
  region_totals <- region_totals %>% arrange(ParentRegion, desc(total))
  region_order <- region_totals$Region
  parent_blocks <- region_totals$ParentRegion
  
  # gaps between ParentRegion blocks
  region_gaps <- rep(1, length(region_order))
  if (length(region_order) > 1) {
    boundaries <- which(parent_blocks[-1] != parent_blocks[-length(parent_blocks)])
    region_gaps[boundaries] <- 6
  }
  if (length(region_gaps) > 0) region_gaps[length(region_gaps)] <- 8
  
  sector_order <- c(group_names, region_order)
  
  # sector colors
  grid.col <- c(
    setNames(rep(group_fill, length(group_names)), group_names),
    setNames(ifelse(region_totals$category == "shared", shared_fill, unique_fill), region_order)
  )
  
  # chord input + per-link colors (match region category)
  chord_df <- links_long %>%
    transmute(from = Region, to = Group, value = Density) %>%
    left_join(region_totals %>% select(Region, category), by = c("from" = "Region")) %>%
    mutate(link_col = ifelse(category == "shared", shared_link, unique_link))
  
  # ---------- Plot ----------
  # Create a separate legend panel so legends can never overlap/cut off.
  if (isTRUE(show_legends)) {
    layout(matrix(c(1, 2), nrow = 1), widths = c(4.5, legend_panel_width))
  } else {
    layout(1)
  }
  
  circos.clear()
  circos.par(
    start.degree = start_degree,
    gap.after = c(rep(12, length(group_names)), region_gaps)
  )
  
  pre_tracks <- if (show_parent_labels) {
    list(
      list(track.height = parent_track_height),  # parent labels
      list(track.height = label_track_height)    # labels
    )
  } else {
    list(
      list(track.height = label_track_height)    # labels
    )
  }
  
  chordDiagram(
    x = chord_df %>% select(from, to, value),
    order = sector_order,
    grid.col = grid.col,
    col = chord_df$link_col,
    transparency = link_alpha,
    directional = 0,
    annotationTrack = "grid",
    preAllocateTracks = pre_tracks
  )
  
  # Track indices depending on whether parent labels are shown
  parent_track_index <- if (show_parent_labels) 1 else NA_integer_
  label_track_index  <- if (show_parent_labels) 2 else 1
  
  # ---------- ParentRegion labels ----------
  if (show_parent_labels) {
    parent_ranges <- region_totals %>%
      group_by(ParentRegion) %>%
      summarise(first_region = first(Region), .groups = "drop")
    
    circos.trackPlotRegion(
      track.index = parent_track_index,
      bg.border = NA,
      panel.fun = function(x, y) {
        sector <- get.cell.meta.data("sector.index")
        if (sector %in% group_names) return(invisible(NULL))
        
        pr <- parent_ranges %>% filter(first_region == sector)
        if (nrow(pr) == 0) return(invisible(NULL))
        
        xcenter <- get.cell.meta.data("xcenter")
        
        circos.text(
          x = xcenter, y = 0.5,
          labels = pr$ParentRegion[1],
          cex = parent_label_cex,
          facing = if (parent_label_facing == "radial") "outside" else "bending.inside",
          niceFacing = TRUE,
          adj = c(0.5, 0.5)
        )
      }
    )
  }
  
  # ---------- Labels track (groups curved; regions tangential) ----------
  label_lookup <- links_long %>%
    select(Region, Label) %>%
    distinct() %>%
    group_by(Region) %>%
    slice_head(n = 1) %>%
    ungroup()
  
  get_label <- function(region) {
    lab <- label_lookup$Label[match(region, label_lookup$Region)]
    ifelse(is.na(lab) | lab == "", region, lab)
  }
  
  circos.trackPlotRegion(
    track.index = label_track_index,
    bg.border = NA,
    panel.fun = function(x, y) {
      sector <- get.cell.meta.data("sector.index")
      xcenter <- get.cell.meta.data("xcenter")
      
      if (sector %in% group_names) {
        circos.text(
          x = xcenter, y = 0.5,
          labels = sector,
          cex = group_label_cex,
          facing = "bending.inside",
          niceFacing = TRUE,
          adj = c(0.5, 0.5)
        )
      } else {
        circos.text(
          x = xcenter, y = 0.5,
          labels = get_label(sector),
          cex = region_label_cex,
          facing = "clockwise",
          niceFacing = TRUE,
          adj = c(0, 0.5)
        )
      }
    }
  )
  
  # ---- Legend panel (panel 2) ----
  if (isTRUE(show_legends)) {
    # Legend panel uses normalized coordinates (0..1), so nothing clips.
    par(mar = c(2, 2, 2, 2))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))

    # Shared vs unique legend
    legend(
      x = 0.05, y = 0.95,
      xjust = 0, yjust = 1,
      bty = "n",
      legend = c("Shared region", "Unique region"),
      fill = c(shared_fill, unique_fill),
      cex = 0.95
    )

    # Density legend: example bar heights proportional to total mean density
    min_total <- suppressWarnings(min(region_totals$total, na.rm = TRUE))
    max_total <- suppressWarnings(max(region_totals$total, na.rm = TRUE))
    if (!is.finite(min_total) || !is.finite(max_total) || max_total <= 0) {
      min_total <- 0
      max_total <- 1
    }
    mid_total <- (min_total + max_total) / 2

    text(0.05, 0.62, "Total mean density\n(bar height scale)", adj = c(0, 1), cex = 0.9)

    vals <- c(min_total, mid_total, max_total)
    if (is.numeric(density_legend_n) && density_legend_n >= 2) {
      vals <- seq(min_total, max_total, length.out = density_legend_n)
    }

    norm <- function(v) max(0, min(1, v / max_total))
    x0 <- 0.12
    dx <- 0.18
    y_base <- 0.18
    h_max <- 0.32

    for (i in seq_along(vals)) {
      xi <- x0 + (i - 1) * dx
      h <- norm(vals[i]) * h_max
      rect(xi, y_base, xi + 0.08, y_base + h, border = NA, col = "grey60")
      text(xi + 0.04, y_base - 0.06, format(round(vals[i], 2), nsmall = 0), cex = 0.85)
    }
  }
  
  invisible(list(
    links_long = links_long,
    presence = presence,
    region_totals = region_totals,
    parent_map_used = parent_map
  ))
  # Reset layout to prevent issues in future plots
  layout(1)
}

## Example usage (uncomment to run)
# make_density_chord_parentgroups(
#   group_files = c(
#     "~/path/to/GroupA/Top30RegionsDepth6_Cleaned.xlsx",
#     "~/path/to/GroupB/Top30RegionsDepth6_Cleaned.xlsx"
#   ),
#   group_names = c("Conditioning", "Fear Recall"),
#   top_regions = 30,
#   show_parent_labels = FALSE,
#   show_legends = TRUE
# )
