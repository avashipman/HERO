# HERO
HERO (Hierarchy-aware Expression Region Organization): A hierarchy-aware analysis pipeline based in R for interpreting rodent whole-brain 3D datasets. HERO is designed to perform hierarchy-aware selection, refinement, ranking, and visualization. 

**HERO** is designed to help users identify, organize, rank, and visualize brain regions with high expression while preserving the anatomical hierarchy of brain atlases. HERO was created with the Allen Brain Atlas-style ontology used in LifeCanvas' **SmartAnalytics** outputs but can be applied to other atlas-registered cell detection datasets. Whole-brain cell detection datasets contain hundreds to thousands of atlas-registered regions, organized as parent-child anatomical structures across multiple depths. This hierarchy is powerful, but it also creates challenges: overrepresentation of subregions, laterality consideration, and biasm. HERO addresses these challenges by using hiararchy-aware filtering, density averaging, overrepresentation correction, and group comparison visualization to bridge the gab to facilitate biological meaningful interpretations and results from whole-brain cell detection datasets .

---

HERO performs five major steps:
1. **Depth Filtering**
    - Import whole-brain cell detection CSV files
    - Clean region names and acronyms
    - Handle laterality: left, right, or combined
    - Filter regions by atlas depth
    - Remove non-target structures (e.g. 'tracts') if necessary
2. **Group-Level Averaging**
    - Combine sample-level filtered outputs
    - Calculate mean density per region across samples
    - Preserve individual subject-level density values
3. **Hierarchical-Guided Region Selection**
    - Identify highly expressing parent regions
    - Hierarchy mapping using parent-child IDs
    - Expand selected regions to a target depth (e.g. Depth 6)
4. **Overrepresentation Correction**
    - Detect overrepresentation from multiple children of a single parent region
    - Collapse overrepresented children back to their parent structure
    - Refill the top-region list with the next best candidate regions
    - Preserve the final selected region list and all candidate regions
5. **Visualization**
    - Generate chord diagrams comparing expression patterns between groups
    - Display shared and group-specific regions
    - Handle laterality: left, right, or combined
    - Represent density using ribbon thickness

  ---
  ## Repository Structure

  This repository is organized around two main R scripts:

  ```text 
  #Contains all reusable HERO functions. Users generally do not need to edit
  #this file unless they want to modify the pipeline itself
  functions_HERO.R

  #The actual script file that the user will input their data and set the parameters
  #needed for their analysis. This script calls functions from the functions_HERO.R
  #and runs the full HERO workflow.
  analysis_HERO.R
```

---
## Input Data

HERO is designed for region-level CSV files exported from SmartAnalytics LifeCanvas. 
Each CSV file should contain, at minimum, the following columns:

| Column                 | Description             |
| ---------------------- | ----------------------- |
| `id`                   | Atlas region ID         |
| `name`                 | Region name             |
| `acronym`              | Region acronym          |
| `parent_structure_id`  | ID of the parent region |
| `depth`                | Atlas hierarchy depth   |
| `count`                | Expression/cell count   |
| `volume (mm^3)`        | Region volume           |
| `density (cells/mm^3)` | Expression density      |

The 'parent_structure_id' comlumn is especially important because HERO uses it to reconstruct the atlas hierarchy and identify child regions.

In order to get sample individual data (e.g. 'M1'), each csv file name needs to put the sample ID in the front like so, 'M1_'. The code is designed to pull the sample ID from the beginning of the file name until the underscore.

---
## Required R Packages

```text
install.packages(c(
  "readxl",
  "writexl",
  "dplyr",
  "tidyr",
  "stringr",
  "circlize"
))

#For visualizations
install.packages(c("ComplexHeatmap", "RColorBrewer"))
```
---
## Quick Start

1. Clone or download the repository
```text
git clone https://github.com/your-username/HERO.git
```
2. Edit analysis_HERO.R to match the paths, files, and parameters for your project
   - Recommended folder organization
     ```
       HERO_project/
       │
       ├── functions_HERO.R
       ├── analysis_HERO.R
       │
       ├── data/
       │   ├── raw_csv/
       │   └── filtered/
       │
       ├── results/
       │   ├── group_A/
       │   ├── group_B/
       │   └── figures/
       │
       └── README.md
     ```
    - Enter your paths and project files    
      ```text
      source("functions_HERO.R")
      
      input_folder <- "path/to/SmartAnalytics_CSVs"
      filtered_folder <- "path/to/filtered_outputs"
      output_folder <- "path/to/final_outputs"
      ```
    - Set your core analysis parameters:
      ```text
      depth_level <- 4
      target_depth <- 6
      laterality <- "combined"
      top_n_parents <- 10
      top_n_final <- 30
      collapse_threshold <- 3
      ```
3. Run through the commands of analysis_HERO.R
   Will give you three main files:
  - Filtered sample files
    - Each file contains the filtered regions for one sample after depth filtering and laterality handling 
  - AveragedRegions.xlsx

    | Tab         | Description                                      |
    | ----------- | ------------------------------------------------ |
    | `Densities` | Individual sample density values for each region |
    | `Averages`  | Mean density across samples for each region      |

  - Top30RegionsDepth6.xlsx
      - This file contains top target-depth regions identified by utilizing the hierarchy from high-density parent regions.

    | Tab                  | Description                                                                 |
    | -------------------- | --------------------------------------------------------------------------- |
    | `Densities`          | Subject-level density values for the final selected regions                 |
    | `Averages`           | Final cleaned top-region list with mean density, acronym, and depth         |
    | `RankedRegions`      | Region-level mean density, median density, and subject-specific dense ranks |
    | `AveragesAllRegions` | Full candidate pool considered during selection and refill                  |
    | `CollapsedGroups`    | Regions collapsed due to parent overrepresentation                          |
 
4. Optional functions in analysis_HERO.R
    - Focused Analysis: start from one or more user-defined parent regions
    - Chord Diagram Visualization: compares region-level expression for 2 or more groups
      - The cord diagram uses cleaned HERO outputs from each group. Links, or ribbons, connect brain regions to groups, and link thickness represents mean density.
      - Regions are categorized by Shared and Unique
        - Shared: Region appears in the top N regions in multiple groups
        - Unique: Region appears in the top N regions in only one group 

---
## Citation
If you use HERO in your work, please cite the associated manuscript:

---
## Contact
For questions, issues, or suggestions, please open a GitHub issue or contact Ava Shipman (ava.shipman@wfusm.edu).

