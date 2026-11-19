# HERO: Hierarchy-aware Expression Region Organization
# Functions for hierarchy-aware analysis of SmartAnalytics/LifeCanvas whole-brain light-sheet expression CSVs.
#
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
#'
#' Internal helper used by HERO functions. Reads one SmartAnalytics CSV, detects
#' hemisphere from the original name/acronym, optionally merges left/right regions,
#' recalculates combined density, and rebuilds parent-child IDs in combined mode.
#'
#' @param csv_file Path to one SmartAnalytics CSV.
#' @param laterality One of "combined", "left", or "right".
#' @param density_col Name of density column.
#' @param exclude_terms Character vector of terms to exclude from region names.
#' @return Data frame with id, parent_structure_id, depth, region_base, acronym_base, density.
#' Build a sample hierarchy with complete ID-based parent mapping.
#' The original atlas IDs are used to resolve parents *before* excluding
#' regions; in combined mode we reassign unique IDs after L/R aggregation.
HERO_build_hierarchy <- function(csv_file,
                                 laterality = c("combined", "left", "right"),
                                 density_col = "density (cells/mm^3)",
                                 exclude_terms = c("tract")) {
  laterality <- match.arg(laterality)
  HERO_check_packages(c("dplyr", "stringr"))
  raw <- utils::read.csv(csv_file, stringsAsFactors = FALSE, check.names = FALSE)
  req <- c("id", "name", "acronym", "parent_structure_id", "depth",
           "count", "volume (mm^3)", density_col)
  missing <- setdiff(req, names(raw))
  if (length(missing)) {
    stop("Missing columns in ", basename(csv_file), ": ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
  raw$region_base <- HERO_strip_lr(as.character(raw$name))
  raw$acronym_base <- HERO_strip_acr_lr(as.character(raw$acronym))
  raw$original_id <- as.character(raw$id)
  raw$original_parent_id <- as.character(raw$parent_structure_id)
  raw$depth <- suppressWarnings(as.integer(raw$depth))
  raw <- raw[!is.na(raw$depth) & raw$depth >= 0L, , drop = FALSE]

  # Unlike a join on raw IDs, this lookup never multiplies counts when
  # left/right records carry the same atlas ID.
  by_id <- split(raw$region_base[!is.na(raw$original_id)],
                 raw$original_id[!is.na(raw$original_id)])
  conflicted <- names(by_id)[vapply(by_id, function(z)
    length(unique(z[!is.na(z)])) > 1L, logical(1))]
  if (length(conflicted)) {
    stop("Conflicting structure names for atlas IDs in ", basename(csv_file),
         ": ", paste(head(conflicted, 10), collapse = ", "), call. = FALSE)
  }
  id_to_name <- vapply(by_id, function(z) {
    v <- unique(z[!is.na(z)]); if (length(v)) v[[1]] else NA_character_
  }, character(1))
  raw$parent_base <- unname(id_to_name[raw$original_parent_id])

  nm <- as.character(raw$name)
  ac <- as.character(raw$acronym)
  hemi_name <- ifelse(grepl("^left\\s+", nm, ignore.case = TRUE), "left",
                     ifelse(grepl("^right\\s+", nm, ignore.case = TRUE), "right", NA_character_))
  hemi_ac <- ifelse(grepl("-L$", ac), "left",
                    ifelse(grepl("-R$", ac), "right", "combined"))
  raw$hemisphere <- ifelse(is.na(hemi_name), hemi_ac, hemi_name)

  excluded <- HERO_region_is_excluded(raw$region_base, exclude_terms)
  raw <- raw[!excluded, , drop = FALSE]
  if (laterality != "combined") {
    chosen <- raw[raw$hemisphere == laterality, , drop = FALSE]
    return(dplyr::transmute(
      chosen,
      id = original_id,
      parent_structure_id = original_parent_id,
      depth = depth,
      region_base = region_base,
      acronym_base = acronym_base,
      density = suppressWarnings(as.numeric(.data[[density_col]])),
      count = suppressWarnings(as.numeric(count)),
      volume = suppressWarnings(as.numeric(`volume (mm^3)`)),
      parent_base = parent_base
    ))
  }

  # Prefer actual L/R measurements over precomputed combined rows. Retain
  # combined-only root or intermediate nodes when their L/R equivalents do
  # not exist; they are needed for valid ancestry across atlas depths.
  lr <- raw[raw$hemisphere %in% c("left", "right"), , drop = FALSE]
  cb <- raw[raw$hemisphere == "combined", , drop = FALSE]
  if (nrow(lr)) {
    cb <- cb[!cb$region_base %in% lr$region_base, , drop = FALSE]
    src <- dplyr::bind_rows(lr, cb)
  } else {
    src <- cb
  }
  if (!nrow(src)) {
    stop("No atlas rows for laterality='combined' in ", basename(csv_file),
         call. = FALSE)
  }
  src$count <- suppressWarnings(as.numeric(src$count))
  src[["volume (mm^3)"]] <- suppressWarnings(as.numeric(src[["volume (mm^3)"]]))
  merged <- src |>
    dplyr::group_by(depth, parent_base, region_base, acronym_base) |>
    dplyr::summarise(
      count = HERO_safe_sum(count),
      volume = HERO_safe_sum(`volume (mm^3)`),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      density = dplyr::if_else(is.finite(volume) & volume > 0 &
                                 is.finite(count), count / volume, NA_real_)
    ) |>
    dplyr::arrange(depth, parent_base, region_base) |>
    dplyr::mutate(id = as.character(dplyr::row_number()))
  # Region is the public workbook key; an ambiguous atlas name must not
  # be mapped to an arbitrary parent or duplicate a subject's measurement.
  dup <- unique(merged$region_base[duplicated(merged$region_base)])
  if (length(dup)) {
    stop("Atlas structure name is ambiguous across depths/parents in ",
         basename(csv_file), ": ", paste(head(dup, 10), collapse = ", "),
         call. = FALSE)
  }
  merged$parent_structure_id <- unname(stats::setNames(merged$id, merged$region_base)[merged$parent_base])
  dplyr::select(merged, id, parent_structure_id, depth, region_base,
                acronym_base, density, count, volume, parent_base)
}

# Fixed-string, case-insensitive exclusion shared by all analysis scopes.
HERO_region_is_excluded <- function(regions, exclude_terms) {
  exclude_terms <- as.character(exclude_terms)
  exclude_terms <- exclude_terms[!is.na(exclude_terms) & nzchar(exclude_terms)]
  if (!length(exclude_terms)) return(rep(FALSE, length(regions)))
  result <- rep(FALSE, length(regions))
  for (term in exclude_terms) {
    result <- result | grepl(term, regions, fixed = TRUE, ignore.case = TRUE)
  }
  result[is.na(result)] <- FALSE
  result
}

HERO_safe_sum <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (!any(is.finite(x))) return(NA_real_)
  sum(x[is.finite(x)])
}


#' Find descendants at a specified atlas depth
#'
#' @param hier Hierarchy table from HERO_build_hierarchy().
#' @param start_ids Starting parent ID(s).
#' @param want_depth Target atlas depth to return.
#' @return Integer vector of descendant IDs at want_depth.
HERO_descendants_at_depth <- function(hier, start_ids, want_depth) {
  ids <- as.character(hier$id)
  parents <- as.character(hier$parent_structure_id)
  depths <- suppressWarnings(as.integer(hier$depth))
  starts <- unique(as.character(start_ids))
  starts <- starts[!is.na(starts) & nzchar(starts)]
  if (!length(starts)) return(character(0))
  if (anyDuplicated(ids)) {
    stop("Duplicate node IDs in hierarchy; cannot safely traverse atlas.", call. = FALSE)
  }
  children <- split(ids[!is.na(parents)], parents[!is.na(parents)])
  queue <- starts
  visited <- character(0)
  result <- character(0)
  while (length(queue)) {
    node <- queue[[1L]]
    queue <- queue[-1L]
    if (node %in% visited) next
    visited <- c(visited, node)
    k <- match(node, ids)
    if (is.na(k) || is.na(depths[[k]])) next
    if (depths[[k]] == want_depth) {
      result <- c(result, node)
    } else if (depths[[k]] < want_depth) {
      next_children <- children[[node]]
      if (length(next_children)) queue <- c(queue, next_children)
    }
  }
  unique(result)
}


#' Rank regions by mean density
HERO_first_nonmissing <- function(x) {
  if (is.character(x)) {
    x <- x[!is.na(x) & nzchar(x)]
  } else {
    x <- x[!is.na(x)]
  }
  if (length(x)) x[[1L]] else if (is.character(x)) NA_character_ else NA_integer_
}

HERO_rank_table <- function(long_df, include_parent_cols = FALSE) {
  HERO_check_packages("dplyr")
  base <- c("Region", "Acronym", "Depth", "MeanDensity")
  parent <- c("ParentRegion", "ParentAcronym", "ParentDepth")
  if (!nrow(long_df)) {
    empty <- dplyr::tibble(Region = character(), Acronym = character(),
                           Depth = integer(), MeanDensity = numeric())
    if (include_parent_cols) {
      empty$ParentRegion <- character()
      empty$ParentAcronym <- character()
      empty$ParentDepth <- integer()
    }
    return(empty)
  }
  if (include_parent_cols) {
    long_df |>
      dplyr::group_by(Region) |>
      dplyr::summarise(
        MeanDensity = HERO_safe_mean(Density),
        Acronym = HERO_first_nonmissing(Acronym),
        Depth = HERO_first_nonmissing(Depth),
        ParentRegion = HERO_first_nonmissing(ParentRegion),
        ParentAcronym = HERO_first_nonmissing(ParentAcronym),
        ParentDepth = HERO_first_nonmissing(ParentDepth),
        .groups = "drop"
      ) |>
      dplyr::arrange(dplyr::desc(MeanDensity), Region)
  } else {
    long_df |>
      dplyr::group_by(Region) |>
      dplyr::summarise(
        MeanDensity = HERO_safe_mean(Density),
        Acronym = HERO_first_nonmissing(Acronym),
        Depth = HERO_first_nonmissing(Depth),
        .groups = "drop"
      ) |>
      dplyr::arrange(dplyr::desc(MeanDensity), Region)
  }
}


# -----------------------------------------------------------------------------
# HERO: audit helpers (atlas-wide overlap, dense subject ranks, selection report)
# -----------------------------------------------------------------------------
HERO_safe_mean <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else mean(x)
}

HERO_safe_median <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else stats::median(x)
}

# One rank column per subject. Pool rows are aggregated once per subject/region;
# absent or invalid measurements remain NA rather than taking the lowest rank.
HERO_make_ranked_regions <- function(long_df) {
  HERO_check_packages(c("dplyr", "tidyr"))
  needed <- c("SampleID", "Region", "Density")
  if (!all(needed %in% names(long_df))) {
    stop("HERO_make_ranked_regions needs columns: ", paste(needed, collapse = ", "),
         call. = FALSE)
  }
  if (nrow(long_df) == 0L) {
    return(dplyr::tibble(Region = character(), MeanDensity = numeric(),
                         MedianDensity = numeric()))
  }
  observed <- long_df |>
    dplyr::transmute(SampleID = as.character(SampleID),
                     Region = as.character(Region),
                     Density = suppressWarnings(as.numeric(Density))) |>
    dplyr::filter(!is.na(SampleID), nzchar(SampleID),
                  !is.na(Region), nzchar(Region)) |>
    dplyr::mutate(Density = dplyr::if_else(is.finite(Density), Density, NA_real_)) |>
    dplyr::group_by(SampleID, Region) |>
    dplyr::summarise(Density = HERO_safe_mean(Density), .groups = "drop")

  summary <- observed |>
    dplyr::group_by(Region) |>
    dplyr::summarise(MeanDensity = HERO_safe_mean(Density),
                     MedianDensity = HERO_safe_median(Density),
                     .groups = "drop")
  rank_wide <- observed |>
    dplyr::group_by(SampleID) |>
    dplyr::mutate(Rank = dplyr::dense_rank(dplyr::desc(Density))) |>
    dplyr::ungroup() |>
    dplyr::select(Region, SampleID, Rank) |>
    tidyr::pivot_wider(names_from = SampleID, values_from = Rank)

  summary |>
    dplyr::left_join(rank_wide, by = "Region") |>
    dplyr::arrange(dplyr::desc(MeanDensity), dplyr::desc(MedianDensity), Region)
}

# Read the COMPLETE, unfiltered atlas links; excluded candidate regions must
# still be allowed to occur as intermediate ancestors. Original IDs are used
# here (not just the immediately annotated parent in an analysis depth pool).
# For combined laterality, left/right nodes are unified by base region name.
HERO_atlas_parent_lookup <- function(csv_files, laterality = "combined") {
  HERO_check_packages("stringr")
  if (length(csv_files) == 0L) stop("No CSVs supplied for atlas validation.", call. = FALSE)
  rows <- list()
  for (csv in csv_files) {
    raw <- utils::read.csv(csv, check.names = FALSE, stringsAsFactors = FALSE)
    required <- c("id", "parent_structure_id", "name", "acronym", "depth")
    missing <- setdiff(required, names(raw))
    if (length(missing)) stop("Cannot validate atlas in ", basename(csv),
                              "; missing ", paste(missing, collapse = ", "), call. = FALSE)
    x <- data.frame(id = as.character(raw$id),
                    parent_id = as.character(raw$parent_structure_id),
                    Region = HERO_strip_lr(as.character(raw$name)),
                    Depth = suppressWarnings(as.integer(raw$depth)),
                    stringsAsFactors = FALSE)
    hemi_name <- ifelse(grepl("^left\\s+", raw$name, ignore.case = TRUE), "left",
                        ifelse(grepl("^right\\s+", raw$name, ignore.case = TRUE), "right", NA_character_))
    hemi_acr <- ifelse(grepl("-L$", raw$acronym), "left",
                       ifelse(grepl("-R$", raw$acronym), "right", "combined"))
    x$hemisphere <- ifelse(is.na(hemi_name), hemi_acr, hemi_name)
    # Some atlas exports reuse the same ontology ID for L/R rows. Resolve IDs
    # to hemisphere-stripped names, and reject only true ID/name conflicts.
    id_region <- split(x$Region[!is.na(x$id)], x$id[!is.na(x$id)])
    ambiguous_id <- names(id_region)[vapply(id_region, function(z) {
      length(unique(z[!is.na(z)])) > 1L
    }, logical(1))]
    if (length(ambiguous_id)) {
      stop("Atlas IDs map to different structures in ", basename(csv), ": ",
           paste(ambiguous_id, collapse = ", "), call. = FALSE)
    }
    id_lookup <- vapply(id_region, function(z) {
      nonmissing <- unique(z[!is.na(z)])
      if (!length(nonmissing)) NA_character_ else nonmissing[[1L]]
    }, character(1))
    # An unresolved intermediate parent would make a no-overlap result unsafe.
    # Root/sentinel links are allowed; missing non-root parents are not.
    unresolved <- !is.na(x$parent_id) & nzchar(x$parent_id) &
      !(x$parent_id %in% c("0", "-1")) & !(x$parent_id %in% names(id_lookup))
    minimum_depth <- suppressWarnings(min(x$Depth, na.rm = TRUE))
    if (any(unresolved & !is.na(x$Depth) & x$Depth > minimum_depth)) {
      stop("Incomplete atlas ancestry in ", basename(csv),
           ": a non-root region references a parent ID not present in the CSV. ",
           "Overlap validation requires the full atlas hierarchy.", call. = FALSE)
    }
    x$ParentRegion <- unname(id_lookup[x$parent_id])
    selected <- if (laterality == "combined") {
      x
    } else {
      x[x$hemisphere == laterality, , drop = FALSE]
    }
    selected <- selected[!is.na(selected$Region) & nzchar(selected$Region), , drop = FALSE]
    rows[[length(rows) + 1L]] <- selected[, c("Region", "ParentRegion", "Depth")]
  }
  all <- do.call(rbind, rows)
  if (is.null(all) || nrow(all) == 0L) {
    stop("No atlas regions for requested laterality: ", laterality, call. = FALSE)
  }
  all$ParentRegion[!is.na(all$ParentRegion) & all$ParentRegion == ""] <- NA_character_
  by_region <- split(all, all$Region)
  # Analysis uses 'Region' as its unique key. A genuinely ambiguous name would
  # invalidate both mean densities and overlap checks; do not silently guess.
  for (nm in names(by_region)) {
    d <- by_region[[nm]]
    parent_names <- unique(d$ParentRegion[!is.na(d$ParentRegion)])
    if (length(parent_names) > 1L) {
      stop("Ambiguous atlas ancestry for region '", nm, "': ",
           paste(parent_names, collapse = ", "), ". Check atlas nomenclature/IDs.",
           call. = FALSE)
    }
  }
  parent_map <- vapply(by_region, function(d) {
    p <- unique(d$ParentRegion[!is.na(d$ParentRegion)])
    if (length(p) == 0L) NA_character_ else p[[1L]]
  }, character(1))
  # Check all chains, including ancestors that will never be output candidates.
  for (nm in names(parent_map)) {
    seen <- character()
    cursor <- nm
    while (!is.na(cursor) && nzchar(cursor)) {
      if (cursor %in% seen) {
        stop("Cycle detected in complete atlas at '", cursor, "'.", call. = FALSE)
      }
      seen <- c(seen, cursor)
      if (!(cursor %in% names(parent_map))) break
      cursor <- unname(parent_map[[cursor]])
    }
  }
  parent_map
}

HERO_ancestors <- function(region, parent_map) {
  if (!(region %in% names(parent_map))) {
    stop("Region missing from complete atlas hierarchy: '", region,
         "'. Cannot safely validate overlaps.", call. = FALSE)
  }
  ancestors <- character()
  cursor <- unname(parent_map[[region]])
  while (!is.na(cursor) && nzchar(cursor)) {
    if (cursor %in% ancestors || identical(cursor, region)) {
      stop("Cycle in atlas ancestor chain for region '", region, "'.", call. = FALSE)
    }
    ancestors <- c(ancestors, cursor)
    if (!(cursor %in% names(parent_map))) break
    cursor <- unname(parent_map[[cursor]])
  }
  ancestors
}

HERO_overlap_pairs <- function(regions, parent_map) {
  regions <- as.character(regions)
  if (anyNA(regions) || any(!nzchar(regions))) {
    stop("Final selection contains an unnamed region.", call. = FALSE)
  }
  if (anyDuplicated(regions)) {
    stop("Final selection contains duplicate regions: ",
         paste(unique(regions[duplicated(regions)]), collapse = ", "), call. = FALSE)
  }
  if (!length(regions)) {
    return(data.frame(Ancestor = character(), Descendant = character()))
  }
  pairs <- lapply(regions, function(child) {
    parents <- intersect(HERO_ancestors(child, parent_map), regions)
    if (!length(parents)) return(NULL)
    data.frame(Ancestor = parents, Descendant = rep(child, length(parents)))
  })
  out <- do.call(rbind, pairs)
  if (is.null(out)) data.frame(Ancestor = character(), Descendant = character()) else out
}

HERO_assert_nonoverlap <- function(selected_regions, parent_map, stage) {
  overlaps <- HERO_overlap_pairs(selected_regions, parent_map)
  if (nrow(overlaps)) {
    stop("HERO anatomical overlap at ", stage, ": ",
         paste0(overlaps$Ancestor, " > ", overlaps$Descendant, collapse = "; "),
         ". Selection stopped to prevent parent–child double counting.", call. = FALSE)
  }
  invisible(TRUE)
}

# Preserves the original depth6 -> depth5 -> start-depth refill precedence.
# All candidate additions (including replacement parents) are tested against
# every already selected region's entire ancestor chain.
HERO_select_cleaned_regions <- function(means6, means5, means_start,
                                        csv_files, laterality, target_depth,
                                        top_n_final, collapse_threshold, max_iter) {
  HERO_check_packages("dplyr")
  if (length(top_n_final) != 1L || is.na(top_n_final) ||
      top_n_final < 1 || top_n_final != as.integer(top_n_final)) {
    stop("top_n_final must be a positive integer.", call. = FALSE)
  }
  cols <- c("Region", "Acronym", "Depth", "MeanDensity",
            "ParentRegion", "ParentAcronym", "ParentDepth")
  for (pool in list(means6, means5, means_start)) {
    if (!all(cols %in% names(pool))) {
      stop("A HERO selection pool is missing required columns.", call. = FALSE)
    }
  }
  top_n_final <- as.integer(top_n_final)
  parent_map <- HERO_atlas_parent_lookup(csv_files, laterality)
  current <- means6[0, cols, drop = FALSE]
  total_rejected <- 0L

  # Adds candidates one-by-one, never truncating a good prior selection in
  # order to add a new higher-depth region that covers an existing region.
  try_add <- function(cur, ranked, excluded = character()) {
    rejected <- 0L
    for (i in seq_len(nrow(ranked))) {
      if (nrow(cur) >= top_n_final) break
      candidate <- as.character(ranked$Region[[i]])
      if (!is.finite(as.numeric(ranked$MeanDensity[[i]]))) next
      if (is.na(candidate) || !nzchar(candidate) ||
          candidate %in% excluded || candidate %in% cur$Region) next
      if (!(candidate %in% names(parent_map))) {
        stop("Candidate '", candidate, "' absent from complete atlas.", call. = FALSE)
      }
      ancestor_set <- HERO_ancestors(candidate, parent_map)
      if (any(cur$Region %in% ancestor_set) ||
          any(vapply(cur$Region, function(s) candidate %in% HERO_ancestors(s, parent_map),
                     logical(1)))) {
        rejected <- rejected + 1L
        next
      }
      cur <- dplyr::bind_rows(cur, ranked[i, cols, drop = FALSE])
    }
    list(current = cur, rejected = rejected)
  }
  collect <- function(cur, ranked, excluded = character()) {
    one <- try_add(cur, ranked, excluded)
    total_rejected <<- total_rejected + one$rejected
    one
  }
  refill <- function(cur, excluded = character()) {
    iteration_rejected <- 0L
    for (pool in list(means6, means5, means_start)) {
      one <- collect(cur, pool, excluded)
      cur <- one$current
      iteration_rejected <- iteration_rejected + one$rejected
      if (nrow(cur) >= top_n_final) break
    }
    list(current = cur, rejected = iteration_rejected)
  }

  initial <- refill(current)
  current <- initial$current
  HERO_assert_nonoverlap(current$Region, parent_map, "initial top-N selection")
  collapsed_parents <- character()
  group_rows <- list()
  iter <- 1L
  repeat {
    over <- current |>
      dplyr::filter(Depth == target_depth, !is.na(ParentRegion)) |>
      dplyr::count(ParentRegion, name = "n_children") |>
      dplyr::filter(n_children > collapse_threshold)
    if (nrow(over) == 0L) break
    if (iter > max_iter) {
      warning("HERO reached max_iter=", max_iter,
              " with overrepresented parents remaining; inspect the selection.",
              call. = FALSE)
      break
    }
    to_collapse <- as.character(over$ParentRegion)
    collapsed_parents <- unique(c(collapsed_parents, to_collapse))
    children <- current |>
      dplyr::filter(Depth == target_depth, ParentRegion %in% to_collapse)
    child_rank <- match(children$Region, means6$Region)
    pd <- means5$MeanDensity[match(children$ParentRegion, means5$Region)]
    cd <- as.numeric(children$MeanDensity)
    valid <- is.finite(cd) & cd > 0 & is.finite(pd) & pd >= 0
    attenuation <- rep(NA_real_, length(cd))
    attenuation[valid] <- 100 * (1 - pd[valid] / cd[valid])
    retention <- ifelse(!valid, "NotEvaluable",
                        ifelse(attenuation > 50, "flagged", "pass"))
    details <- data.frame(
      RecordType = rep("CollapsedChild", nrow(children)),
      ParentRegion = children$ParentRegion,
      ParentAcronym = children$ParentAcronym,
      ParentDepth = children$ParentDepth,
      ChildRegion = children$Region,
      ChildAcronym = children$Acronym,
      ChildDepth = children$Depth,
      ChildMeanDensity = cd,
      ChildOriginalRank = as.integer(child_rank),
      ParentMeanDensity = as.numeric(pd),
      DensityAttenuation = attenuation,
      SignalRetention = retention,
      Iteration = rep(iter, nrow(children)),
      stringsAsFactors = FALSE
    )
    kept <- current[!(current$Region %in% children$Region), cols, drop = FALSE]
    parent_candidates <- means5[means5$Region %in% to_collapse, cols, drop = FALSE]
    missing_parents <- setdiff(to_collapse, parent_candidates$Region)
    if (length(missing_parents)) {
      stop("Cannot collapse selected children: parent density is unavailable for ",
           paste(missing_parents, collapse = ", "),
           ". Preserve these children and inspect the atlas/filters.", call. = FALSE)
    }
    added_parents <- collect(kept, parent_candidates)
    current <- added_parents$current

    # Do not refill with any descendants of a parent that was collapsed, even
    # if an intervening hierarchy level has been filtered out by exclude_terms.
    all_candidate_names <- unique(c(means6$Region, means5$Region,
                                    means_start$Region))
    excluded_descendants <- all_candidate_names[vapply(all_candidate_names, function(candidate) {
      any(collapsed_parents %in% HERO_ancestors(candidate, parent_map))
    }, logical(1))]
    filled <- refill(current, excluded = excluded_descendants)
    current <- filled$current
    iter_rejected <- added_parents$rejected + filled$rejected
    details$OverlapOccurred <- rep(iter_rejected > 0L, nrow(details))
    details$OverlapCandidatesRejected <- rep(as.integer(iter_rejected), nrow(details))
    group_rows[[length(group_rows) + 1L]] <- details
    HERO_assert_nonoverlap(current$Region, parent_map,
                          paste0("collapse/refill iteration ", iter))
    iter <- iter + 1L
  }

  HERO_assert_nonoverlap(current$Region, parent_map, "final selection")
  # A final validation record exists independently of whether a collapse occurred.
  # Iteration=0 denotes the final overall audit, not a collapse iteration.
  validation <- data.frame(
    RecordType = "SelectionValidation", ParentRegion = NA_character_,
    ParentAcronym = NA_character_, ParentDepth = NA_integer_,
    ChildRegion = NA_character_, ChildAcronym = NA_character_,
    ChildDepth = NA_integer_, ChildMeanDensity = NA_real_,
    ChildOriginalRank = NA_integer_, ParentMeanDensity = NA_real_,
    DensityAttenuation = NA_real_, SignalRetention = NA_character_,
    Iteration = 0L, OverlapOccurred = total_rejected > 0L,
    OverlapCandidatesRejected = as.integer(total_rejected),
    stringsAsFactors = FALSE
  )
  collapsed_groups <- do.call(dplyr::bind_rows, c(group_rows, list(validation)))
  selected_n <- nrow(current)
  missing_n <- max(0L, top_n_final - selected_n)
  if (missing_n > 0L) {
    warning("HERO top-N incomplete: requested ", top_n_final, ", selected ",
            selected_n, "; ", missing_n,
            " position(s) unfilled because no additional anatomically non-overlapping " ,
            "eligible regions were available.", call. = FALSE)
  }
  selection_report <- data.frame(
    RequestedTopN = top_n_final, SelectedN = selected_n,
    MissingN = missing_n,
    Status = if (missing_n == 0L) "Complete" else "Incomplete",
    OverlapValidation = "Pass",
    OverlapCandidatesRejected = as.integer(total_rejected),
    stringsAsFactors = FALSE
  )
  list(current = current |> dplyr::arrange(dplyr::desc(MeanDensity)),
       collapsed_groups = collapsed_groups,
       collapsed_parents = collapsed_parents,
       selection_report = selection_report)
}


HERO_empty_candidate_long <- function() {
  dplyr::tibble(SampleID = character(), Region = character(),
                Acronym = character(), Depth = integer(),
                Density = numeric(), parent_id = character(),
                ParentRegion = character(), ParentAcronym = character(),
                ParentDepth = integer())
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
    ranked_regions <- HERO_make_ranked_regions(data.frame(
      SampleID = base_name,
      Region = output_df$region_base,
      Density = output_df$density,
      stringsAsFactors = FALSE
    ))
    # First tab stays the filtered source table for existing read_excel(f) calls.
    writexl::write_xlsx(list("FilteredRegions" = output_df,
                             "RankedRegions" = ranked_regions), out_file)
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
  
  xlsx_files <- list.files(filtered_folder, pattern = "\\.(xlsx|xls)$",
                           full.names = TRUE, ignore.case = TRUE)
  canonical <- grepl("_FilteredDepth[0-9]+_(combined|left|right)\\.xlsx$",
                     basename(xlsx_files), ignore.case = TRUE)
  if (any(canonical)) xlsx_files <- xlsx_files[canonical]
  if (!length(xlsx_files)) {
    stop("No filtered Excel files found in filtered_folder: ", filtered_folder)
  }
  
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
      dplyr::summarise(density = HERO_safe_mean(density), .groups = "drop")
    sample_long_list[[length(sample_long_list)+1]] <- tmp
    
    if (all(c("id", "depth") %in% names(df))) {
      meta_long_list[[length(meta_long_list)+1]] <- df |>
        dplyr::transmute(sample_id = sample_id,
                         region = as.character(.data[[region_col]]),
                         id = as.character(id),
                         depth = suppressWarnings(as.integer(depth))) |>
        dplyr::filter(!is.na(region) & region != "") |>
        dplyr::group_by(sample_id, region) |>
        dplyr::summarise(id = HERO_first_nonmissing(id),
                         depth = HERO_first_nonmissing(depth),
                         .groups = "drop")
    }
  }
  
  long_all <- dplyr::bind_rows(sample_long_list)
  wide <- long_all |> tidyr::pivot_wider(names_from = region, values_from = density)
  region_means <- wide |>
    dplyr::select(-sample_id) |>
    dplyr::summarise(dplyr::across(dplyr::everything(), HERO_safe_mean))
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
  ranked_regions <- HERO_make_ranked_regions(long_all |>
    dplyr::transmute(SampleID = sample_id, Region = region, Density = density))
  writexl::write_xlsx(list("Densities" = wide_ordered,
                           "Averages" = averages_df,
                           "RankedRegions" = ranked_regions), out_path)
  invisible(out_path)
}



# -----------------------------------------------------------------------------
# Function 3: Top30Regions
# -----------------------------------------------------------------------------
#' Select uncleaned top hierarchy-aware expression regions
#'
#' Starting from the highest-density parent regions in AveragedRegions.xlsx, this
#' function traverses the atlas hierarchy to a target depth and writes the
#' uncleaned top-region output before any parent-collapse correction is applied.
#' This output is useful when users want to inspect the raw top target-depth
#' regions selected by density.
#'
#' @param averaged_file Path to AveragedRegions.xlsx from RegionsbyDensity().
#' @param original_csv_folder Folder containing original SmartAnalytics CSVs.
#' @param filtered_folder_for_id_map Folder containing FilterbyDepth outputs.
#' @param output_file Output Excel file path, usually Top30RegionsDepth6.xlsx.
#' @param depth_start Starting atlas depth used for parent selection.
#' @param target_depth Target descendant depth, usually 6.
#' @param top_n_parents Number of starting parent regions to expand.
#' @param top_n_final Number of uncleaned top regions to report.
#' @param laterality One of "combined", "left", or "right".
#' @param density_col Name of density column in raw CSVs.
#' @param exclude_terms Terms to exclude from region names; default excludes "tract".
#' @return Invisibly returns a list containing densities, averages, and the full candidate table.
Top30Regions <- function(averaged_file,
                         original_csv_folder,
                         filtered_folder_for_id_map,
                         output_file = file.path(dirname(averaged_file), "Top30RegionsDepth6.xlsx"),
                         depth_start = 4,
                         target_depth = 6,
                         top_n_parents = 20,
                         top_n_final = 30,
                         laterality = c("combined", "left", "right"),
                         density_col = "density (cells/mm^3)",
                         exclude_terms = c("tract")) {
  laterality <- match.arg(laterality)
  HERO_check_packages(c("readxl", "writexl", "dplyr", "tidyr", "stringr"))
  
  # ---- Step 1: top parent regions from AveragedRegions.xlsx ----
  avg_df <- readxl::read_excel(averaged_file, sheet = "Averages")
  if ("Region" %in% names(avg_df)) {
    avg_df <- avg_df |>
      dplyr::rename(region = Region, mean_density = MeanDensity)
  }
  if (!all(c("region", "mean_density") %in% names(avg_df))) {
    stop("Sheet 'Averages' must include columns named region/mean_density or Region/MeanDensity.", call. = FALSE)
  }
  avg_df <- avg_df |>
    dplyr::filter(!HERO_region_is_excluded(region, exclude_terms))

  top_parents <- avg_df |>
    dplyr::arrange(dplyr::desc(mean_density)) |>
    dplyr::slice_head(n = top_n_parents) |>
    dplyr::pull(region)
  
  # ---- Step 2: map top parent names to IDs for left/right mode ----
  map_files <- list.files(filtered_folder_for_id_map, pattern = "\\.(xlsx|xls)$", full.names = TRUE)
  filtered_candidates <- grepl("_FilteredDepth[0-9]+_(combined|left|right)\\.xlsx$",
                               basename(map_files), ignore.case = TRUE)
  if (any(filtered_candidates)) map_files <- map_files[filtered_candidates]
  if (length(map_files) == 0) stop("No .xlsx/.xls files found in filtered_folder_for_id_map: ", filtered_folder_for_id_map, call. = FALSE)
  
  map_df <- readxl::read_excel(map_files[1])
  region_col <- HERO_pick_region_col(map_df, c("region_base", "Region", "region", "name"))
  if (!all(c("id", "depth") %in% names(map_df))) {
    stop("Filtered-by-depth files must include id and depth columns. Re-run FilterbyDepth().", call. = FALSE)
  }
  
  parent_id_map <- map_df |>
    dplyr::transmute(region = as.character(.data[[region_col]]),
                     id = as.character(id),
                     depth = as.integer(depth)) |>
    dplyr::filter(depth == depth_start, region %in% top_parents)
  
  if (laterality %in% c("left", "right") && nrow(parent_id_map) == 0) {
    stop("Could not map top parent regions to IDs at depth_start=", depth_start,
         ". Check that FilterbyDepth used the same depth and laterality.", call. = FALSE)
  }
  
  # ---- Step 3: collect target-depth descendants under the top parents for each sample ----
  csv_files <- list.files(original_csv_folder, pattern = "\\.csv$", full.names = TRUE)
  if (length(csv_files) == 0) stop("No .csv files found in original_csv_folder: ", original_csv_folder, call. = FALSE)
  
  long_list <- list()
  parent_lookup_list <- list()
  
  for (csv in csv_files) {
    sample_id <- tools::file_path_sans_ext(basename(csv))
    hier <- HERO_build_hierarchy(csv,
                                 laterality = laterality,
                                 density_col = density_col,
                                 exclude_terms = exclude_terms)
    
    parent_ids <- if (laterality %in% c("left", "right")) {
      parent_id_map$id
    } else {
      hier |>
        dplyr::filter(depth == depth_start, region_base %in% top_parents) |>
        dplyr::pull(id)
    }
    
    target_ids <- unique(unlist(lapply(parent_ids, function(pid) {
      HERO_descendants_at_depth(hier, pid, target_depth)
    }), use.names = FALSE))
    
    if (length(target_ids) == 0) next
    
    long_list[[length(long_list) + 1]] <- hier |>
      dplyr::filter(id %in% target_ids, depth == target_depth) |>
      dplyr::transmute(
        SampleID = sample_id,
        Region = region_base,
        Acronym = acronym_base,
        Depth = depth,
        Density = density
      )
  }
  
  long_all <- dplyr::bind_rows(long_list)
  if (nrow(long_all) == 0) {
    stop("No target-depth descendants found. Check depth_start, target_depth, laterality, and input folders.", call. = FALSE)
  }
  
  # ---- Step 4: build a depth-general parent lookup table ----
  # ParentRegion should mean the immediate anatomical parent of each region,
  # regardless of whether that region is at target_depth, target_depth - 1,
  # depth_start, or any other atlas depth. This prevents AveragesAllRegions
  # from having parent columns only for Depth 6 rows.
  for (csv in csv_files) {
    hier <- HERO_build_hierarchy(csv,
                                 laterality = laterality,
                                 density_col = density_col,
                                 exclude_terms = exclude_terms)
    
    parent_lookup_list[[length(parent_lookup_list) + 1]] <- hier |>
      dplyr::select(
        Region = region_base,
        Acronym = acronym_base,
        Depth = depth,
        id,
        parent_structure_id
      ) |>
      dplyr::left_join(
        hier |>
          dplyr::select(
            parent_structure_id = id,
            ParentRegion = region_base,
            ParentAcronym = acronym_base,
            ParentDepth = depth
          ),
        by = "parent_structure_id"
      ) |>
      dplyr::select(Region, Acronym, Depth, ParentRegion, ParentAcronym, ParentDepth) |>
      dplyr::filter(!is.na(Region))
  }
  
  parent_lookup_all <- dplyr::bind_rows(parent_lookup_list) |>
    dplyr::group_by(Region, Acronym, Depth, ParentRegion, ParentAcronym, ParentDepth) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
    dplyr::arrange(Region, Acronym, Depth, dplyr::desc(n)) |>
    dplyr::group_by(Region, Acronym, Depth) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup() |>
    dplyr::select(-n)
  
  # ---- Step 5: rank target-depth candidates and keep uncleaned top N ----
  region_means_all <- long_all |>
    dplyr::group_by(Region) |>
    dplyr::summarise(
      MeanDensity = mean(Density, na.rm = TRUE),
      Acronym = dplyr::first(stats::na.omit(Acronym)),
      Depth = dplyr::first(Depth),
      .groups = "drop"
    ) |>
    dplyr::left_join(parent_lookup_all, by = c("Region", "Acronym", "Depth")) |>
    dplyr::arrange(dplyr::desc(MeanDensity))
  
  region_means <- region_means_all |>
    dplyr::slice_head(n = top_n_final)
  
  top_regions <- region_means$Region
  
  # ---- Step 5: preserve subject-level values for selected uncleaned regions ----
  wide <- long_all |>
    dplyr::filter(Region %in% top_regions) |>
    dplyr::select(SampleID, Region, Density) |>
    dplyr::group_by(SampleID, Region) |>
    dplyr::summarise(Density = mean(Density, na.rm = TRUE), .groups = "drop") |>
    tidyr::pivot_wider(names_from = Region, values_from = Density)
  
  wide <- wide[, c("SampleID", top_regions), drop = FALSE]
  
  writexl::write_xlsx(
    list(
      "Densities" = wide,
      "Averages" = region_means,
      "AveragesAllRegions" = region_means_all,
      "RankedRegions" = HERO_make_ranked_regions(long_all)
    ),
    path = output_file
  )
  
  invisible(list(densities = wide, averages = region_means, averages_all_regions = region_means_all, long = long_all))
}

# -----------------------------------------------------------------------------
# Function 4: Top30Cleaned
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
#' @param exclude_terms Character vector of terms to exclude from region names; default excludes "tract".
#' @param max_iter Maximum collapse/refill iterations.
#' @return Invisibly returns a list with final tables.
# The cleaned selection has an independent atlas-overlap audit, density
# attenuation audit, and requested-top-N completeness report.

Top30Cleaned <- function(averaged_file,
                         original_csv_folder,
                         filtered_folder_for_id_map,
                         output_file = file.path(dirname(averaged_file), "Top30RegionsDepth6_Cleaned.xlsx"),
                         depth_start = 4,
                         target_depth = 6,
                         top_n_parents = 20,
                         top_n_final = 30,
                         collapse_threshold = 3,     # > 3 children => collapse
                         laterality = c("combined", "left", "right"),
                         density_col = "density (cells/mm^3)",
                         exclude_terms = c("tract"),
                         max_iter = 10
)
{
  laterality <- match.arg(laterality)
  if (depth_start >= target_depth) {
    stop("Top30Cleaned requires depth_start < target_depth. ",
         "For example use depth_start = 4 and target_depth = 6. ",
         "Equal depths cannot provide the parent/child collapse pools.",
         call. = FALSE)
  }
  
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
  
  # Helper: identify regions to exclude from all Top30Cleaned candidate pools.
  # Terms are treated as plain text and matched case-insensitively.
  exclude_terms <- exclude_terms[!is.na(exclude_terms) & nzchar(exclude_terms)]
  region_is_excluded <- function(x) {
    if (length(exclude_terms) == 0) return(rep(FALSE, length(x)))
    hits <- lapply(exclude_terms, function(term) {
      grepl(term, x, ignore.case = TRUE, fixed = TRUE)
    })
    Reduce(`|`, hits)
  }
  
  
  build_hierarchy <- function(csv_file) {
    HERO_build_hierarchy(csv_file, laterality = laterality,
                         density_col = density_col,
                         exclude_terms = exclude_terms)
  }
  descendants_at_depth <- HERO_descendants_at_depth

  # ---- Determine top depth_start parents from AveragedRegions.xlsx ----
  avg_df <- readxl::read_excel(averaged_file, sheet = "Averages")
  if ("Region" %in% names(avg_df)) {
    avg_df <- avg_df |> dplyr::rename(region = Region, mean_density = MeanDensity)
  }
  avg_df <- avg_df |>
    dplyr::filter(!region_is_excluded(region))
  top_parents <- avg_df |>
    dplyr::arrange(dplyr::desc(mean_density)) |>
    dplyr::slice_head(n = top_n_parents) |>
    dplyr::pull(region)
  
  # ---- For left/right, map top parents to IDs using a filtered file ----
  map_files <- list.files(filtered_folder_for_id_map, pattern="\\.xlsx$", full.names=TRUE)
  filtered_candidates <- grepl("_FilteredDepth[0-9]+_(combined|left|right)\\.xlsx$",
                               basename(map_files), ignore.case = TRUE)
  if (any(filtered_candidates)) map_files <- map_files[filtered_candidates]
  if (length(map_files) == 0) stop("No .xlsx files found in filtered_folder_for_id_map")
  map_df <- readxl::read_excel(map_files[1])
  region_col <- if ("region_base" %in% names(map_df)) "region_base" else if ("name" %in% names(map_df)) "name" else NA_character_
  if (is.na(region_col)) stop("Filtered file needs 'region_base' or 'name'.")
  if (!all(c("id","depth") %in% names(map_df))) stop("Filtered-by-depth files must include id + depth.")
  
  parent_id_map <- map_df |>
    dplyr::transmute(region = as.character(.data[[region_col]]),
                     id = as.character(id),
                     depth = as.integer(depth)) |>
    dplyr::filter(depth == depth_start,
                  region %in% top_parents,
                  !region_is_excluded(region))
  
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
                         parent_id = parent_structure_id) |>
        dplyr::filter(!region_is_excluded(Region))
    }
    if (!length(out)) {
      return(HERO_empty_candidate_long() |>
               dplyr::select(SampleID, Region, Acronym, Depth, Density, parent_id))
    }
    dplyr::bind_rows(out)
  }

  long6 <- long_by_depth(target_depth)
  long5 <- long_by_depth(target_depth - 1)
  long4 <- long_by_depth(depth_start) # last-resort fallback depth
  parent_lookup_list <- list()
  
  if (nrow(long6) == 0) stop("No depth-", target_depth, " candidates found under top parents.")
  
  # Parent map for all candidate depths.
  # ParentRegion should mean the immediate anatomical parent of each region,
  # regardless of whether that row is Depth 6, Depth 5, depth_start, or another
  # atlas depth. This keeps AveragesAllRegions interpretable for mixed-depth
  # candidate pools created by collapse/refill logic.
  for (csv in csv_files) {
    hier <- build_hierarchy(csv)
    
    parent_lookup_list[[length(parent_lookup_list) + 1]] <- hier |>
      dplyr::select(
        Region = region_base,
        Acronym = acronym_base,
        Depth = depth,
        id,
        parent_structure_id
      ) |>
      dplyr::left_join(
        hier |>
          dplyr::select(
            parent_structure_id = id,
            ParentRegion = region_base,
            ParentAcronym = acronym_base,
            ParentDepth = depth
          ),
        by = "parent_structure_id"
      ) |>
      dplyr::select(Region, Acronym, Depth, ParentRegion, ParentAcronym, ParentDepth) |>
      dplyr::filter(!is.na(Region))
  }
  
  parent_lookup_all <- dplyr::bind_rows(parent_lookup_list) |>
    dplyr::group_by(Region, Acronym, Depth, ParentRegion, ParentAcronym, ParentDepth) |>
    dplyr::summarise(n = dplyr::n(), .groups = "drop") |>
    dplyr::arrange(Region, Acronym, Depth, dplyr::desc(n)) |>
    dplyr::group_by(Region, Acronym, Depth) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup() |>
    dplyr::select(-n)
  
  annotate_parent_info <- function(x) {
    x |>
      dplyr::select(-dplyr::any_of(c("ParentRegion", "ParentAcronym", "ParentDepth"))) |>
      dplyr::left_join(parent_lookup_all, by = c("Region", "Acronym", "Depth"))
  }
  
  long6 <- annotate_parent_info(long6)
  long5 <- annotate_parent_info(long5)
  long4 <- annotate_parent_info(long4)
  
  # Shared ranking for all candidate depth pools.
  rank_table <- HERO_rank_table

  means6 <- rank_table(long6, include_parent_cols = TRUE) |>
    dplyr::filter(!region_is_excluded(Region))
  means5 <- rank_table(long5, include_parent_cols = TRUE) |>
    dplyr::filter(!region_is_excluded(Region))
  means4 <- rank_table(long4, include_parent_cols = TRUE) |>
    dplyr::filter(!region_is_excluded(Region))
  
  # Atlas-wide, overlap-safe selection with collapse/refill audit.
  cleaned_selection <- HERO_select_cleaned_regions(
    means6 = means6, means5 = means5, means_start = means4,
    csv_files = csv_files, laterality = laterality,
    target_depth = target_depth, top_n_final = top_n_final,
    collapse_threshold = collapse_threshold, max_iter = max_iter
  )
  current <- cleaned_selection$current
  collapsed_groups_all <- cleaned_selection$collapsed_groups
  collapsed_parents_all <- cleaned_selection$collapsed_parents
  selection_report <- cleaned_selection$selection_report

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
      ParentRegion, ParentAcronym, ParentDepth,
      Pool = paste0("Depth", target_depth - 1),
      ExcludedBecauseCollapsed = FALSE
    )
  
  avg_all_4 <- means4 |>
    dplyr::transmute(
      Region, Acronym, Depth, MeanDensity,
      ParentRegion, ParentAcronym, ParentDepth,
      Pool = paste0("Depth", depth_start, "_LastResort"),
      ExcludedBecauseCollapsed = FALSE
    )
  
  averages_all_regions <- dplyr::bind_rows(avg_all_6, avg_all_5, avg_all_4) |>
    dplyr::arrange(dplyr::desc(MeanDensity))
  
  ranked_regions <- HERO_make_ranked_regions(dplyr::bind_rows(
    long6 |> dplyr::select(SampleID, Region, Density),
    long5 |> dplyr::select(SampleID, Region, Density),
    long4 |> dplyr::select(SampleID, Region, Density)
  ))

  writexl::write_xlsx(
    list(
      "Densities" = wide,
      "Averages"  = cleaned_list,
      "AveragesAllRegions" = averages_all_regions,
      "RankedRegions" = ranked_regions,
      "CollapsedGroups" = collapsed_groups_all,
      "SelectionReport" = selection_report
    ),
    path = output_file
  )
  
  invisible(list(densities = wide,
                 averages = cleaned_list,
                 collapsed_groups = collapsed_groups_all,
                 collapsed_parents = collapsed_parents_all,
                 selection_report = selection_report))
}

# -----------------------------------------------------------------------------
# Function 5: Top30FocusedCleaned (focused parent-region analysis)
# -----------------------------------------------------------------------------
# Function to select top 30 child regions from specified starting parent
# regions defined by the user. Then uses the same process as Function 4.
# -----------------------------------------------------------------------------
#Focused filtering for whole-brain data analysis
Top30FocusedCleaned <- function(original_csv_folder,
                                start_region,
                                output_file = NULL,
                                start_depth = NULL,          # NULL, scalar, named vector/list, or same length as start_region
                                target_depth = 6,
                                top_n_final = 30,
                                collapse_threshold = 3,      # > 3 depth-6 children from same depth-5 parent triggers collapse
                                laterality = c("combined", "left", "right"),
                                density_col = "density (cells/mm^3)",
                                exclude_terms = c("tract"),
                                max_iter = 10,
                                case_insensitive_match = TRUE) {
  
  laterality <- match.arg(laterality)
  
  # Allow one or multiple starting regions.
  # Examples:
  #   start_region = "Pallidum"
  #   start_region = c("Pallidum", "Hypothalamus")
  start_region <- as.character(start_region)
  if (length(start_region) == 0 || any(is.na(start_region)) || any(start_region == "")) {
    stop("start_region must contain at least one non-empty region name.")
  }
  start_region_clean <- unique(start_region)
  
  # If no output_file is provided, create one from all starting regions.
  if (is.null(output_file)) {
    safe_region_name <- paste(gsub("[^A-Za-z0-9_]+", "_", start_region_clean), collapse = "_")
    output_file <- file.path(original_csv_folder, paste0("Top30_From_", safe_region_name, "_Cleaned.xlsx"))
  }
  
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
  
  # Helper: identify regions to exclude from all Top30FocusedCleaned candidate pools.
  # Terms are treated as plain text and matched case-insensitively.
  exclude_terms <- exclude_terms[!is.na(exclude_terms) & nzchar(exclude_terms)]
  region_is_excluded <- function(x) {
    if (length(exclude_terms) == 0) return(rep(FALSE, length(x)))
    hits <- lapply(exclude_terms, function(term) {
      grepl(term, x, ignore.case = TRUE, fixed = TRUE)
    })
    Reduce(`|`, hits)
  }
  
  match_region <- function(x, regions = start_region_clean) {
    if (case_insensitive_match) {
      tolower(x) %in% tolower(regions)
    } else {
      x %in% regions
    }
  }
  
  # Optional region-specific start depth.
  # Accepted forms:
  #   NULL: automatically choose broadest/minimum depth per matched region
  #   scalar: use same depth for all start regions
  #   vector same length as start_region: paired in order
  #   named vector/list: names are region names
  get_depth_for_region <- function(region_name) {
    if (is.null(start_depth)) return(NULL)
    
    if (!is.null(names(start_depth)) && region_name %in% names(start_depth)) {
      return(as.integer(start_depth[[region_name]]))
    }
    
    if (length(start_depth) == 1) return(as.integer(start_depth[[1]]))
    
    if (length(start_depth) == length(start_region_clean)) {
      idx <- match(region_name, start_region_clean)
      return(as.integer(start_depth[[idx]]))
    }
    
    stop("start_depth must be NULL, length 1, a named vector/list, or the same length as start_region.")
  }
  
  # Use the same full hierarchy traversal as every other HERO analysis.
  build_hierarchy <- function(csv_file) {
    HERO_build_hierarchy(csv_file, laterality = laterality,
                         density_col = density_col,
                         exclude_terms = exclude_terms)
  }
  descendants_at_depth <- HERO_descendants_at_depth


  # Shared ranking for all candidate depth pools.
  rank_table <- HERO_rank_table

  # ---------------- Main ----------------
  csv_files <- list.files(original_csv_folder, pattern = "\\.csv$", full.names = TRUE)
  if (length(csv_files) == 0) stop("No .csv files found in: ", original_csv_folder)
  
  long6_list <- list()
  long5_list <- list()
  longS_list <- list()  # last-resort: start depth(s) themselves
  parent_lookup_list <- list()
  resolved_start_depths <- integer(0)
  
  for (csv in csv_files) {
    sample_id <- tools::file_path_sans_ext(basename(csv))
    hier <- build_hierarchy(csv)
    
    # Build an immediate-parent lookup for this sample's hierarchy.
    # This is depth-general, so Depth 6, Depth 5, and start-depth fallback
    # regions all receive their anatomical parent in the output.
    parent_lookup_this <- hier |>
      transmute(
        parent_id = id,
        ParentRegion = region_base,
        ParentAcronym = acronym_base,
        ParentDepth = depth
      )
    
    # Find starting node(s) for each requested start region.
    # If multiple regions are supplied, descendants from all are pooled before ranking.
    candidate_list <- list()
    
    for (sr in start_region_clean) {
      sr_candidates <- hier |> filter(match_region(region_base, sr))
      if (nrow(sr_candidates) == 0) next
      
      sr_depth <- get_depth_for_region(sr)
      if (!is.null(sr_depth)) {
        sr_candidates <- sr_candidates |> filter(depth == sr_depth)
      } else {
        # Choose the broadest/minimum depth match for this specific start region
        min_d <- min(sr_candidates$depth, na.rm = TRUE)
        sr_candidates <- sr_candidates |> filter(depth == min_d)
        resolved_start_depths <- c(resolved_start_depths, min_d)
      }
      
      if (nrow(sr_candidates) > 0) {
        sr_candidates$StartRegion <- sr
        candidate_list[[length(candidate_list) + 1]] <- sr_candidates
      }
    }
    
    if (length(candidate_list) == 0) next
    candidates <- dplyr::bind_rows(candidate_list)
    
    start_ids <- unique(candidates$id)
    if (length(start_ids) == 0) next
    
    # Multiple start regions can have different depths, so this can be a vector.
    this_start_depths <- sort(unique(candidates$depth))
    
    # Depth 6 descendants
    ids6 <- descendants_at_depth(hier, start_ids, target_depth)
    if (length(ids6) > 0) {
      long6_list[[length(long6_list) + 1]] <- hier |>
        filter(id %in% ids6, depth == target_depth) |>
        transmute(
          SampleID = sample_id,
          Region = region_base,
          Acronym = acronym_base,
          Depth = depth,
          Density = density,
          parent_id = parent_structure_id
        ) |>
        left_join(parent_lookup_this, by = "parent_id")
    }
    
    # Depth 5 descendants
    ids5 <- descendants_at_depth(hier, start_ids, target_depth - 1)
    if (length(ids5) > 0) {
      long5_list[[length(long5_list) + 1]] <- hier |>
        filter(id %in% ids5, depth == (target_depth - 1)) |>
        transmute(
          SampleID = sample_id,
          Region = region_base,
          Acronym = acronym_base,
          Depth = depth,
          Density = density,
          parent_id = parent_structure_id
        ) |>
        left_join(parent_lookup_this, by = "parent_id")
    }
    
    # Last-resort pool: nodes at the starting depth(s) within this subtree.
    idsS <- unique(unlist(lapply(this_start_depths, function(sd) {
      descendants_at_depth(hier, start_ids, sd)
    }), use.names = FALSE))
    
    if (length(idsS) > 0) {
      longS_list[[length(longS_list) + 1]] <- hier |>
        filter(id %in% idsS, depth %in% this_start_depths) |>
        transmute(
          SampleID = sample_id,
          Region = region_base,
          Acronym = acronym_base,
          Depth = depth,
          Density = density,
          parent_id = parent_structure_id
        ) |>
        left_join(parent_lookup_this, by = "parent_id")
    }
  }
  
  empty_long <- HERO_empty_candidate_long()
  long6 <- if (length(long6_list)) bind_rows(long6_list) else empty_long
  long5 <- if (length(long5_list)) bind_rows(long5_list) else empty_long
  longS <- if (length(longS_list)) bind_rows(longS_list) else empty_long
  
  if (nrow(long6) == 0) {
    stop("No depth-", target_depth, " descendants found for start_region(s): ",
         paste(start_region_clean, collapse = ", "), ". ",
         "Check spelling/case, laterality, start_depth, and that the region(s) exist in these CSVs.")
  }
  
  # Ranked pools. Parent columns are included for all pools so that
  # AveragesAllRegions contains immediate parent annotations for every depth.
  means6 <- rank_table(long6, include_parent_cols = TRUE)
  means5 <- rank_table(long5, include_parent_cols = TRUE)
  meansS <- rank_table(longS, include_parent_cols = TRUE)
  
  # Atlas-wide, overlap-safe selection with collapse/refill audit.
  cleaned_selection <- HERO_select_cleaned_regions(
    means6 = means6, means5 = means5, means_start = meansS,
    csv_files = csv_files, laterality = laterality,
    target_depth = target_depth, top_n_final = top_n_final,
    collapse_threshold = collapse_threshold, max_iter = max_iter
  )
  current <- cleaned_selection$current
  collapsed_groups_all <- cleaned_selection$collapsed_groups
  collapsed_parents_all <- cleaned_selection$collapsed_parents
  selection_report <- cleaned_selection$selection_report

  # Final cleaned list
  cleaned_list <- current |> select(Region, Acronym, Depth, MeanDensity)
  
  # ---------------- AveragesAllRegions sheet ----------------
  excluded6_final <- if (length(collapsed_parents_all) > 0) {
    means6 |> filter(ParentRegion %in% collapsed_parents_all) |> pull(Region)
  } else character()
  
  avg_all_6 <- means6 |>
    transmute(Region, Acronym, Depth, MeanDensity,
              ParentRegion, ParentAcronym, ParentDepth,
              Pool = paste0("Depth", target_depth),
              ExcludedBecauseCollapsed = Region %in% excluded6_final)
  
  avg_all_5 <- means5 |>
    transmute(Region, Acronym, Depth, MeanDensity,
              ParentRegion, ParentAcronym, ParentDepth,
              Pool = paste0("Depth", target_depth - 1),
              ExcludedBecauseCollapsed = FALSE)
  
  # Start-depth pool (last resort); can include multiple depths if multiple start regions were used.
  avg_all_S <- meansS |>
    transmute(Region, Acronym, Depth, MeanDensity,
              ParentRegion, ParentAcronym, ParentDepth,
              Pool = "StartDepth_LastResort",
              ExcludedBecauseCollapsed = FALSE)
  
  averages_all_regions <- bind_rows(avg_all_6, avg_all_5, avg_all_S) |>
    arrange(desc(MeanDensity))
  
  # ---------------- Densities sheet (per-sample values) ----------------
  final_regions <- cleaned_list$Region
  
  long_final <- bind_rows(
    long6 |> select(SampleID, Region, Density),
    long5 |> select(SampleID, Region, Density),
    longS |> select(SampleID, Region, Density)
  ) |>
    filter(Region %in% final_regions) |>
    group_by(SampleID, Region) |>
    summarise(Density = mean(Density, na.rm = TRUE), .groups = "drop")
  
  wide <- long_final |>
    tidyr::pivot_wider(names_from = Region, values_from = Density)
  
  wide <- wide[, c("SampleID", final_regions), drop = FALSE]
  
  # CollapsedGroups sheet (ensure it has headers even if empty)
  if (nrow(collapsed_groups_all) == 0) {
    collapsed_groups_all <- tibble(
      ParentRegion = character(), ParentAcronym = character(), ParentDepth = integer(),
      ChildRegion = character(), ChildAcronym = character(), ChildDepth = integer(),
      ChildMeanDensity = numeric(), Iteration = integer()
    )
  }
  
  start_regions_sheet <- data.frame(
    StartRegion = start_region_clean,
    StartDepth = vapply(start_region_clean, function(sr) {
      sd <- get_depth_for_region(sr)
      if (is.null(sd)) NA_integer_ else as.integer(sd)
    }, integer(1)),
    stringsAsFactors = FALSE
  )
  
  writexl::write_xlsx(
    list(
      "Densities" = wide,
      "Averages" = cleaned_list,
      "AveragesAllRegions" = averages_all_regions,
      "RankedRegions" = HERO_make_ranked_regions(dplyr::bind_rows(
        long6 |> dplyr::select(SampleID, Region, Density),
        long5 |> dplyr::select(SampleID, Region, Density),
        longS |> dplyr::select(SampleID, Region, Density)
      )),
      "CollapsedGroups" = collapsed_groups_all,
      "SelectionReport" = selection_report,
      "StartRegions" = start_regions_sheet
    ),
    path = output_file
  )
  
  invisible(list(
    densities = wide,
    averages = cleaned_list,
    averages_all_regions = averages_all_regions,
    collapsed_groups = collapsed_groups_all,
    collapsed_parents = collapsed_parents_all,
    selection_report = selection_report,
    start_regions = start_regions_sheet
  ))
}

# Example usage:
# Top30FocusedCleaned(
#   original_csv_folder = "Pathtofolderwithrawcelldetectioncsvs",
#   start_region = c("Pallidum", "Hypothalamus"),  # must match region_base naming
#   start_depth = 4,                                # one value for all start regions
#   laterality = "combined",
#   output_file = "filepathandfilename",
#   collapse_threshold = 30
# )
#
# You can also provide region-specific depths:
# start_depth = c(Pallidum = 4, Hypothalamus = 3)

# -----------------------------------------------------------------------------
# Function 6: make_density_chord_parentgroups (group comparison chord diagram)
# -----------------------------------------------------------------------------
## Chord diagram analysis to compare region overlap between groups
## Cody by: Ava Shipman, Centanni Lab
## Cleaned up + legend-safe layout (no legend overlap with circos)

## Updated: selects top N Group-Region pairs per group by default

make_density_chord_parentgroups <- function(group_files,
                                            group_names = NULL,
                                            sheet_avg = "Averages",
                                            sheet_all = "AveragesAllRegions",
                                            region_col = "Region",
                                            value_col = "MeanDensity",
                                            parent_col = "ParentRegion",
                                            acronym_col = "Acronym",
                                            # Display / region selection:
                                            # Default behavior: select the top N Group-Region pairs WITHIN EACH GROUP.
                                            # This prevents one group from inheriting extra links from regions
                                            # selected by another group.
                                            top_regions_per_group = 15,
                                            # Legacy option: if top_regions_per_group = NULL, select top N
                                            # regions overall by summed density across groups.
                                            top_regions = NULL,
                                            min_density = NULL,
                                            use_acronym_labels = TRUE,
                                            show_parent_labels = TRUE,
                                            # Styling:
                                            shared_fill = "red",
                                            unique_fill = "blue",
                                            unique_group_fills = NULL,
                                            group_fill = unique_group_fills,
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
                                            show_legends = FALSE,
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
  
  # ---------- Colors for unique regions ----------
  # Shared regions remain shared_fill/shared_link.
  # Unique regions are colored by the group in which they appear so viewers can
  # immediately identify which experimental group contributed each unique region.
  if (is.null(unique_group_fills)) {
    unique_group_fills <- grDevices::hcl.colors(length(group_names), palette = "Dark 3")
    names(unique_group_fills) <- group_names
  } else {
    if (is.null(names(unique_group_fills))) {
      if (length(unique_group_fills) != length(group_names)) {
        stop("If unique_group_fills is unnamed, it must have one color per group_name.", call. = FALSE)
      }
      names(unique_group_fills) <- group_names
    }
    missing_cols <- setdiff(group_names, names(unique_group_fills))
    if (length(missing_cols) > 0) {
      stop("unique_group_fills is missing colors for: ", paste(missing_cols, collapse = ", "), call. = FALSE)
    }
    unique_group_fills <- unique_group_fills[group_names]
  }
  
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
  
  # ---------- Select regions to plot ----------
  # Primary behavior: choose the top N regions WITHIN EACH GROUP and keep only
  # those selected Group-Region pairs.
  #
  # Important: do NOT filter by the union of selected regions alone. If we keep
  # every group's density for the union, then a group can display many more than
  # top_regions_per_group links because it will inherit links for regions selected
  # by other groups. semi_join() below prevents that by preserving only the pairs
  # selected for each group.
  selected_pairs <- NULL
  if (!is.null(top_regions_per_group)) {
    selected_pairs <- links_long %>%
      group_by(Group) %>%
      arrange(desc(Density), .by_group = TRUE) %>%
      slice_head(n = top_regions_per_group) %>%
      ungroup() %>%
      distinct(Group, Region)
    
    links_long <- links_long %>%
      semi_join(selected_pairs, by = c("Group", "Region"))
  }
  
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
  # A region is "shared" if it appears in more than one group.
  # A region is "unique" if it appears in only one group; UniqueGroup stores
  # the group identity used for group-specific unique-region ribbon colors.
  presence <- links_long %>%
    distinct(Region, Group) %>%
    group_by(Region) %>%
    summarise(
      n_groups = dplyr::n(),
      UniqueGroup = ifelse(n_groups == 1, dplyr::first(Group), NA_character_),
      .groups = "drop"
    ) %>%
    mutate(category = ifelse(n_groups > 1, "shared", "unique"))
  
  # ---------- Totals for ordering + legend range ----------
  region_totals <- links_long %>%
    group_by(Region) %>%
    summarise(total = sum(Density, na.rm = TRUE), .groups = "drop") %>%
    left_join(parent_map, by = "Region") %>%
    left_join(presence %>% select(Region, category, UniqueGroup), by = "Region") %>%
    mutate(
      ParentRegion = ifelse(is.na(ParentRegion), "Unmapped", ParentRegion),
      category = ifelse(is.na(category), "unique", category)
    ) %>%
    arrange(desc(total))
  
  # Optional legacy behavior: only select top N overall when per-group selection is off.
  if (is.null(top_regions_per_group) && !is.null(top_regions)) {
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
  # Group sectors remain group_fill by default. Region sectors are red if shared
  # and group-colored if unique.
  unique_region_cols <- unique_group_fills[region_totals$UniqueGroup]
  unique_region_cols[is.na(unique_region_cols)] <- unique_fill
  region_sector_cols <- ifelse(region_totals$category == "shared", shared_fill, unique_region_cols)
  
  grid.col <- c(
    setNames(rep(group_fill, length(group_names)), group_names),
    setNames(region_sector_cols, region_order)
  )
  
  # chord input + per-link colors
  # Shared ribbons are red. Unique ribbons use the color assigned to the group
  # that contributed that unique region.
  chord_df <- links_long %>%
    transmute(from = Region, to = Group, value = Density) %>%
    left_join(region_totals %>% select(Region, category, UniqueGroup), by = c("from" = "Region")) %>%
    mutate(
      link_col = ifelse(
        category == "shared",
        shared_link,
        unique_group_fills[as.character(to)]
      ),
      link_col = ifelse(is.na(link_col), unique_link, link_col)
    )
  
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
    
    # Shared + unique-by-group legend
    legend(
      x = 0.05, y = 0.95,
      xjust = 0, yjust = 1,
      bty = "n",
      legend = c("Shared region", paste0("Unique region: ", names(unique_group_fills))),
      fill = c(shared_fill, unique_group_fills),
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
      rect(xi, y_base, xi + 0.08, y_base + h, border = NA, col = unique_group_fills)
      text(xi + 0.04, y_base - 0.06, format(round(vals[i], 2), nsmall = 0), cex = 0.85)
    }
  }
  
  invisible(list(
    links_long = links_long,
    presence = presence,
    region_totals = region_totals,
    parent_map_used = parent_map,
    selected_regions = region_order,
    selected_pairs = selected_pairs,
    selection_summary = links_long %>% dplyr::count(Group, name = "n_displayed_links"),
    selection_settings = list(
      top_regions_per_group = top_regions_per_group,
      top_regions_overall = top_regions,
      shared_color = shared_fill,
      unique_group_fills = unique_group_fills
    )
  ))
  # Reset layout to prevent issues in future plots
  layout(1)
}



