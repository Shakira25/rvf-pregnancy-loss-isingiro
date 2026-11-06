# =============================================================================
# FILE:    download_MODIS_NDVI_rainfall_v2.1_FDA_17May2026.R
# VERSION: 3.0  (12 May 2026)
#
# PURPOSE: Download MODIS NDVI (via NASA AppEEARS REST API) for each
#          sub-county of Isingiro District, Uganda, using actual sub-county
#          centroids derived from the geoBoundaries ADM4 shapefile.
#
# APPROACH — SUB-COUNTY CENTROIDS:
#   geoBoundaries-UGA-ADM4 contains 1,521 Uganda sub-county polygons. We:
#     1. Read Isingiro district boundary from ADM3 shapefile
#     2. Spatially intersect ADM4 polygons with the Isingiro boundary
#     3. Compute centroids of the matched sub-county polygons (~30 sub-counties)
#     4. Submit sub-county centroids to AppEEARS in ONE task
#     5. After download, cluster the sub-county NDVI profiles using Functional
#        Data Analysis (FDA) — Fourier basis + functional PCA + k-means on fPCA
#        scores — to identify livestock vs. crop vegetation zones by named sub-county
#
# CHANGES FROM v2 (12 May 2026 — grid approach):
#   1. Sample points now = actual sub-county centroids from ADM4, not a
#      regular grid. Each point is labelled with its sub-county name.
#   2. ADM4 shapefile added: geoBoundaries-UGA-ADM4.shp (Admin 4 folder).
#   3. DTW cluster map now labels each point by sub-county name.
#   4. AppEEARS task renamed: Isingiro_RVF_NDVI_SubCounty_v3.
#   5. Task ID saved to appeears_task_id_v3.txt (separate from v2 grid run).
#   6. DTW-PAM clustering replaced with Functional Data Analysis (FDA): Fourier
#      basis smoothing, fPCA, k-means on scores — better suited for phenological
#      cycle analysis.
#
# CHANGES FROM v1 (Abel Wilson Walekhwa feedback, 27 Apr 2026):
#   1. All Isingiro sub-counties covered (Abel: "run for each sub-county").
#   2. Extended period: 1995-01-01 to 2026-03-31 (MODIS from 2000 onward).
#   3. CHIRPS removed — rainfall embedded in NDVI signal.
#   4. Cyclic NDVI climatology (12 monthly means) replaces single-value mean.
#   5. NA imputation deferred to after FDA clustering (within-cluster strategy).
#   6. Per-sub-county lambda_env multiplier for spatial heterogeneity in ODEs.
#
# STUDY:   Estimating the Hidden Burden: Modelling RVF-Attributable Excess
#          Pregnancy Loss Among Women in Isingiro District, Uganda
# AUTHOR:  Shakira Babirye
# SUPERVISORS: Prof. Betty Kivumbi & Abel Wilson Walekhwa
#
# HOW TO RUN:
#   Open in RStudio → Ctrl+Shift+Enter. Run time: 45–90 min (AppEEARS wait).
#   All outputs saved to BASE_DIR/data/.
#
# OUTPUT FILES:
#   NDVI_MODIS/                            — raw AppEEARS CSV
#   isingiro_subcounty_points.csv          — sub-county centroids used
#   isingiro_subcounty_map.png             — map of sub-county centroids
#   ndvi_subcounty_16day_raw.csv           — cleaned 16-day NDVI, all sub-counties
#   ndvi_subcounty_monthly.csv             — monthly mean NDVI per sub-county
#   ndvi_subcounty_monthly_imputed.csv     — after NA imputation
#   ndvi_climatology_monthly.csv           — 12-month climatological mean per sub-county
#   ndvi_timeseries_plot.png               — multi-panel time series by cluster
#
#   FDA-specific outputs (saved to FDA_OUT_DIR = outputs_FDA_17May2026/):
#   fda_cluster_assignments.csv            — FDA cluster for each sub-county
#   fda_cluster_map.png                    — spatial map with sub-county names
#   ndvi_cluster_profiles.png             — mean NDVI profile per FDA cluster
#   lambda_env_fda_cluster.csv             — lambda_env multiplier per cluster × month
#   lambda_env_plot.png                    — seasonal multiplier per FDA cluster
#   fda_fpca_variance.png                  — fPCA variance explained plot
# =============================================================================


# -----------------------------------------------------------------------------
# SECTION 0: INSTALL AND LOAD PACKAGES
# -----------------------------------------------------------------------------
required <- c(
  "sf",          # read shapefile, generate spatial grid
  "httr",        # HTTP requests to AppEEARS API
  "jsonlite",    # JSON encoding
  "tidyverse",   # data wrangling + ggplot2
  "lubridate",   # date arithmetic
  "terra",       # raster handling (dependency)
  "curl",        # HTTP dependency
  "fda",         # Functional Data Analysis (Ramsay & Silverman)
  "cluster",     # silhouette() for cluster quality
  "patchwork",   # combine ggplot panels
  "viridis"      # colour-blind-safe palettes
)

for (pkg in required) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# fda masks dplyr::select() and dplyr::filter() — re-pin to dplyr after loading
select <- dplyr::select
filter <- dplyr::filter


# -----------------------------------------------------------------------------
# SECTION 1: CONFIGURATION
# -----------------------------------------------------------------------------

# [07Oct2026] Paths re-pointed to the self-contained RVF v2.1 folder (copy of this script).
# The original copy in code/ still uses the old locations. Old lines kept for record:
# BASE_DIR    <- "D:/Shakira Tasks/UNIPH Modelling course/Abel new revisions"
# DATA_DIR    <- file.path(BASE_DIR, "data")
# [05 Nov 2026] Local E: drive path replaced so the script runs from the repository root on any computer.
# BASE_DIR    <- "E:/Shakira Tasks/UNIPH Modelling course/Abel new revisions/code/RVF_v2.1_revised_07Oct2026"
BASE_DIR    <- normalizePath(".", winslash = "/")
DATA_DIR    <- file.path(BASE_DIR, "data")
NDVI_DIR    <- file.path(DATA_DIR, "NDVI_MODIS")
# ADM3: district-level boundaries — used to extract the Isingiro polygon
# SHAPE_DIR   <- file.path(BASE_DIR, "Admin 3 shape files uganda")
SHAPE_DIR   <- file.path(DATA_DIR, "shapefiles", "Admin 3 shape files uganda")
SHAPE_FILE  <- file.path(SHAPE_DIR, "geoBoundaries-UGA-ADM3.shp")
# ADM4: sub-county boundaries — 1,521 Uganda sub-counties; used to get centroids
# SHAPE_DIR4  <- file.path(BASE_DIR, "Admin 4 shape files uganda")
SHAPE_DIR4  <- file.path(DATA_DIR, "shapefiles", "Admin 4 shape files uganda")
SHAPE_FILE4 <- file.path(SHAPE_DIR4, "geoBoundaries-UGA-ADM4.shp")

dir.create(NDVI_DIR, recursive = TRUE, showWarnings = FALSE)

# FDA_OUT_DIR <- file.path(BASE_DIR, "outputs_FDA_17May2026")
FDA_OUT_DIR <- file.path(DATA_DIR, "ndvi_fda_outputs")   # read by RVF_Models_Isingiro_Uganda_v2.1_07Oct2026.Rmd
dir.create(FDA_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# [05 Nov 2026] Earthdata credentials removed from the public copy of this script (they were typed in here).
# Before re-running, set them as environment variables, e.g. in ~/.Renviron:
#   EARTHDATA_USER=your_username
#   EARTHDATA_PASS=your_password
EARTHDATA_USER <- Sys.getenv("EARTHDATA_USER")
EARTHDATA_PASS <- Sys.getenv("EARTHDATA_PASS")

DATE_START <- "1995-01-01"   # MODIS data available from 2000-02-18 onwards
DATE_END   <- "2026-03-31"

APPEEARS_BASE <- "https://appeears.earthdatacloud.nasa.gov/api"



# =============================================================================
# SECTION 2: EXTRACT ISINGIRO BOUNDARY AND SUB-COUNTY CENTROIDS
# =============================================================================

cat("── Section 2: Extracting Isingiro sub-county centroids from ADM4...\n")

# Read Isingiro district boundary from ADM3
uga_adm3 <- sf::st_read(SHAPE_FILE, quiet = TRUE)
isingiro_poly <- uga_adm3 |> filter(shapeName == "Isingiro")
if (nrow(isingiro_poly) == 0)
  stop("'Isingiro' not found in ADM3 shapefile.")
cat(sprintf("   Isingiro district polygon found (ADM3).\n"))

# Read ADM4 sub-county boundaries
uga_adm4 <- sf::st_read(SHAPE_FILE4, quiet = TRUE)
uga_adm4 <- sf::st_transform(uga_adm4, sf::st_crs(isingiro_poly))
cat(sprintf("   ADM4 loaded: %d Uganda sub-county polygons.\n", nrow(uga_adm4)))

# Compute sub-county centroids, keep only those whose centroid falls
# strictly inside the Isingiro district polygon (excludes neighbouring
# districts whose polygons touch the boundary, e.g. Mbarara City divisions)
adm4_centroids <- sf::st_centroid(uga_adm4)
within_idx     <- sf::st_within(adm4_centroids,
                                sf::st_union(isingiro_poly),
                                sparse = FALSE)[, 1]

isingiro_subs  <- uga_adm4[within_idx, ]
sub_centroids  <- adm4_centroids[within_idx, ]

cat(sprintf("   Isingiro sub-counties found: %d\n", nrow(isingiro_subs)))
cat(paste0("     ", sort(isingiro_subs$shapeName), collapse = "\n"), "\n\n")

# Build SAMPLE_POINTS with sub-county names
sub_coords <- sf::st_coordinates(sub_centroids)
SAMPLE_POINTS <- data.frame(
  id        = sprintf("ISI_%03d", seq_len(nrow(sub_coords))),
  subcounty = isingiro_subs$shapeName,
  longitude = round(sub_coords[, "X"], 6),
  latitude  = round(sub_coords[, "Y"], 6),
  stringsAsFactors = FALSE
)

cat(sprintf("   %d sample points created (one per sub-county centroid).\n",
            nrow(SAMPLE_POINTS)))
cat(sprintf("   lon range: [%.4f, %.4f]  lat range: [%.4f, %.4f]\n\n",
            min(SAMPLE_POINTS$longitude), max(SAMPLE_POINTS$longitude),
            min(SAMPLE_POINTS$latitude),  max(SAMPLE_POINTS$latitude)))

write_csv(SAMPLE_POINTS, file.path(DATA_DIR, "isingiro_subcounty_points.csv"))

# ── Map of sub-county centroids on district boundary ─────────────────────────
p_subcounty_map <- ggplot() +
  geom_sf(data = isingiro_poly, fill = "#f0f4f0", colour = "grey40",
          linewidth = 0.8) +
  geom_point(data = SAMPLE_POINTS,
             aes(longitude, latitude), colour = "#e74c3c",
             size = 2.5, alpha = 0.85) +
  geom_text(data = SAMPLE_POINTS,
            aes(longitude, latitude, label = subcounty),
            size = 2.0, vjust = -0.6, colour = "grey20") +
  labs(
    title    = "MODIS NDVI Sub-County Centroids — Isingiro District, Uganda",
    subtitle = sprintf("%d sub-county centroids | ADM4 from geoBoundaries | %s to %s",
                       nrow(SAMPLE_POINTS), DATE_START, DATE_END),
    x = "Longitude", y = "Latitude",
    caption  = "One NDVI time series extracted per sub-county; FDA phenology clustering (Fourier basis + fPCA) identifies seasonal vegetation zones"
  ) +
  theme_minimal(base_size = 11) +
  coord_sf()

ggsave(file.path(DATA_DIR, "isingiro_subcounty_map.png"),
       p_subcounty_map, width = 12, height = 10, dpi = 150)
cat("   Saved: isingiro_subcounty_map.png\n\n")


# =============================================================================
# SECTION 3: HELPER FUNCTIONS — RESILIENT HTTP + AUTH
#
# WHY: The AppEEARS task for 25-30 points over 1995-2026 takes 45-90 min to
# process. A single connection can time out or reset during that wait. These
# helpers add:
#   - Automatic retry (up to MAX_RETRIES attempts with exponential back-off)
#   - Re-authentication when a token expires mid-run
#   - Task ID saved to disk immediately after submission (so the script can
#     be restarted without resubmitting the task)
# =============================================================================

TASK_ID_FILE <- file.path(DATA_DIR, "appeears_task_id_v3.txt")
MAX_RETRIES  <- 5       # retry count for transient network errors
POLL_SLEEP   <- 90      # seconds between status polls (90 s is fine; reduces noise)

# ── Authenticate and return a fresh auth header ───────────────────────────────
appeears_auth <- function() {
  cat("   Authenticating with NASA AppEEARS...\n")
  for (attempt in seq_len(MAX_RETRIES)) {
    resp <- tryCatch(
      httr::POST(
        url    = paste0(APPEEARS_BASE, "/login"),
        httr::authenticate(EARTHDATA_USER, EARTHDATA_PASS),
        httr::timeout(120)
      ),
      error = function(e) {
        cat(sprintf("   Auth attempt %d failed (%s). Retrying in 30 s...\n",
                    attempt, conditionMessage(e)))
        Sys.sleep(30); NULL
      }
    )
    if (!is.null(resp) && httr::status_code(resp) == 200) {
      tok <- httr::content(resp)$token
      cat("   Authentication successful.\n")
      return(httr::add_headers(Authorization = paste("Bearer", tok)))
    }
  }
  stop("Authentication failed after ", MAX_RETRIES, " attempts.")
}

# ── Resilient GET: retries on timeout / connection reset ─────────────────────
# The ... argument allows passing extra httr configs such as write_disk() or
# progress() for file downloads while still getting retry behaviour.
safe_GET <- function(url, auth_hdr, timeout_s = 120, ...) {
  for (attempt in seq_len(MAX_RETRIES)) {
    resp <- tryCatch(
      httr::GET(url, config = auth_hdr, httr::timeout(timeout_s), ...),
      error = function(e) {
        wait <- 30 * attempt
        cat(sprintf("   GET failed (attempt %d/%d): %s\n   Waiting %d s...\n",
                    attempt, MAX_RETRIES, conditionMessage(e), wait))
        Sys.sleep(wait); NULL
      }
    )
    if (!is.null(resp)) return(resp)
  }
  stop("GET ", url, " failed after ", MAX_RETRIES, " attempts.")
}


# =============================================================================
# SECTION 4: AUTHENTICATE
# =============================================================================

cat("── Section 4: Authenticating...\n")
auth_header <- appeears_auth()
cat("\n")


# =============================================================================
# SECTION 5: SUBMIT TASK — OR RESUME FROM SAVED TASK ID
#
# If appeears_task_id.txt already exists the task was already submitted
# (possibly from a previous run that timed out). We skip re-submission and
# go straight to polling. This prevents duplicate tasks on AppEEARS.
# =============================================================================

task_id_file_exists <- file.exists(TASK_ID_FILE)

if (task_id_file_exists) {
  task_id <- trimws(readLines(TASK_ID_FILE, n = 1))
  cat(sprintf("── Section 5: Resuming — found saved task ID: %s\n", task_id))
  cat("   Skipping submission. Polling for completion.\n\n")

} else {
  cat(sprintf("── Section 5: Submitting NDVI task for %d sub-county centroids...\n",
              nrow(SAMPLE_POINTS)))

  coords_list <- lapply(seq_len(nrow(SAMPLE_POINTS)), function(i) {
    list(
      id        = SAMPLE_POINTS$id[i],
      longitude = SAMPLE_POINTS$longitude[i],
      latitude  = SAMPLE_POINTS$latitude[i],
      category  = SAMPLE_POINTS$subcounty[i]
    )
  })

  task_body <- list(
    task_type = "point",
    task_name = "Isingiro_RVF_NDVI_SubCounty_v3",
    params = list(
      dates = list(
        list(
          startDate = format(as.Date(DATE_START), "%m-%d-%Y"),
          endDate   = format(as.Date(DATE_END),   "%m-%d-%Y")
        )
      ),
      layers = list(
        list(product = "MOD13Q1.061", layer = "_250m_16_days_NDVI"),
        list(product = "MOD13Q1.061", layer = "_250m_16_days_EVI"),
        list(product = "MOD13Q1.061", layer = "_250m_16_days_pixel_reliability")
      ),
      coordinates = coords_list,
      output = list(format = list(type = "csv"), projection = "geographic")
    )
  )

  task_resp <- tryCatch(
    httr::POST(
      url    = paste0(APPEEARS_BASE, "/task"),
      body   = jsonlite::toJSON(task_body, auto_unbox = TRUE),
      config = auth_header,
      httr::content_type_json(),
      httr::timeout(120)
    ),
    error = function(e) stop("Task submission error: ", conditionMessage(e))
  )

  if (httr::status_code(task_resp) != 202) {
    cat("Task submission response:\n")
    print(httr::content(task_resp))
    stop("Task submission failed (HTTP ", httr::status_code(task_resp), ").")
  }

  task_id <- httr::content(task_resp)$task_id

  # Save immediately — if the script crashes during polling this file lets us resume
  writeLines(task_id, TASK_ID_FILE)
  cat(sprintf("   Task submitted. ID: %s\n", task_id))
  cat(sprintf("   Task ID saved to: %s\n", TASK_ID_FILE))
  cat("   Monitor at: https://appeears.earthdatacloud.nasa.gov/task\n\n")
}


# =============================================================================
# SECTION 6: POLL UNTIL TASK COMPLETE
#
# Re-authenticates every ~90 min (token lifespan) to keep the session alive.
# Catches and retries on connection resets / timeouts instead of crashing.
# =============================================================================

cat("── Section 6: Polling for task completion...\n")
cat(sprintf("   Task ID: %s\n", task_id))
cat("   Checking every 90 seconds. This may take 45–90 min for a full grid.\n\n")

auth_refresh_counter <- 0    # re-auth every 60 polls (~90 min)

repeat {
  # Re-authenticate periodically to avoid token expiry
  auth_refresh_counter <- auth_refresh_counter + 1
  if (auth_refresh_counter %% 60 == 0) {
    cat("   Refreshing authentication token...\n")
    auth_header <- appeears_auth()
  }

  resp <- tryCatch(
    safe_GET(paste0(APPEEARS_BASE, "/task/", task_id), auth_header, timeout_s = 120),
    error = function(e) {
      cat(sprintf("   Poll error: %s — will retry next cycle.\n",
                  conditionMessage(e)))
      NULL
    }
  )

  if (!is.null(resp)) {
    task_status <- httr::content(resp)$status
    cat(sprintf("   [%s] Status: %s\n", format(Sys.time(), "%H:%M:%S"),
                task_status))

    if (task_status == "done") {
      cat("   Task complete!\n\n")
      break
    }
    if (task_status %in% c("error", "deleted")) {
      stop("AppEEARS task failed with status: ", task_status,
           "\nCheck https://appeears.earthdatacloud.nasa.gov/task for details.")
    }
  }

  Sys.sleep(POLL_SLEEP)
}


# =============================================================================
# SECTION 7: DOWNLOAD RESULT FILES
# =============================================================================

cat("── Section 7: Downloading NDVI result files...\n")

files_resp <- safe_GET(paste0(APPEEARS_BASE, "/bundle/", task_id),
                        auth_header, timeout_s = 120)
files_list <- httr::content(files_resp)$files

for (f in files_list) {
  file_name <- basename(f$file_name)
  dest_path <- file.path(NDVI_DIR, file_name)
  if (file.exists(dest_path)) {
    cat(sprintf("   Skipping (exists): %s\n", file_name)); next
  }
  cat(sprintf("   Downloading: %s\n", file_name))
  safe_GET(
    url       = paste0(APPEEARS_BASE, "/bundle/", task_id, "/", f$file_id),
    auth_hdr  = auth_header,
    timeout_s = 300,
    httr::write_disk(dest_path, overwrite = TRUE),
    httr::progress()
  )
}

httr::POST(paste0(APPEEARS_BASE, "/logout"), config = auth_header,
           httr::timeout(30))
cat("   Logged out.\n\n")


# =============================================================================
# SECTION 8: CLEAN AND PROCESS NDVI DATA
# =============================================================================

cat("── Section 8: Processing NDVI CSV...\n")

ndvi_csv <- list.files(NDVI_DIR, pattern = "MOD13Q1.*\\.csv$",
                        full.names = TRUE)

if (length(ndvi_csv) == 0)
  stop("No NDVI CSV found in ", NDVI_DIR, ". Check that the task completed.")

ndvi_raw <- read_csv(ndvi_csv[1], show_col_types = FALSE)
cat(sprintf("   Raw CSV: %d rows, %d columns\n", nrow(ndvi_raw), ncol(ndvi_raw)))

# Standardise and quality-filter
ndvi_clean <- ndvi_raw |>
  rename_with(tolower) |>
  filter(id %in% SAMPLE_POINTS$id) |>
  mutate(
    date     = as.Date(date, format = "%Y-%m-%d"),
    point_id = id,
    ndvi_raw = `mod13q1_061__250m_16_days_ndvi`,
    evi      = `mod13q1_061__250m_16_days_evi`,
    qc_flag  = `mod13q1_061__250m_16_days_pixel_reliability`,
    # Quality mask: 0 = good, 1 = marginal (keep); 2 = snow, 3 = cloud (remove)
    ndvi_qc  = if_else(qc_flag %in% c(0, 1), ndvi_raw, NA_real_),
    month    = floor_date(date, "month"),
    year     = year(date),
    month_n  = month(date)
  ) |>
  select(date, month, year, month_n, point_id,
         latitude, longitude, ndvi_raw, ndvi_qc, evi, qc_flag) |>
  arrange(point_id, date)

write_csv(ndvi_clean, file.path(DATA_DIR, "ndvi_subcounty_16day_raw.csv"))
cat(sprintf("   Saved: ndvi_subcounty_16day_raw.csv  (%d rows, %d sub-counties)\n",
            nrow(ndvi_clean), n_distinct(ndvi_clean$point_id)))

# Monthly aggregation
ndvi_monthly <- ndvi_clean |>
  group_by(point_id, latitude, longitude, month, year, month_n) |>
  summarise(
    ndvi_mean  = mean(ndvi_qc, na.rm = TRUE),
    ndvi_n_obs = sum(!is.na(ndvi_qc)),
    ndvi_sd    = sd(ndvi_qc, na.rm = TRUE),
    evi_mean   = mean(evi, na.rm = TRUE),
    .groups    = "drop"
  ) |>
  mutate(
    ndvi_flag = if_else(ndvi_n_obs < 1 | is.na(ndvi_mean), "MISSING", "OK"),
    source    = "MODIS_MOD13Q1.061_250m"
  )

write_csv(ndvi_monthly, file.path(DATA_DIR, "ndvi_subcounty_monthly.csv"))
cat(sprintf("   Saved: ndvi_subcounty_monthly.csv  (%d point × month rows)\n\n",
            nrow(ndvi_monthly)))


# =============================================================================
# SECTION 9: NA IMPUTATION (WITHIN-CLUSTER STRATEGY)
#
# Strategy:
#   1. LOCF (Last Observation Carried Forward) as an immediate placeholder.
#   2. After FDA clustering (Section 10), replace cluster-member imputations
#      with the within-cluster mean for that calendar month — this is more
#      appropriate than LOCF for the long 1995-2026 series.
#   3. For isolated missing periods (>3 months), use the climatological mean
#      for that point and month.
#
# Abel's guidance: do NOT use the district-wide mean. Use cluster-neighbours.
# =============================================================================

cat("── Section 9: NA imputation (LOCF placeholder; update after clustering)...\n")

ndvi_monthly_imp <- ndvi_monthly |>
  arrange(point_id, month) |>
  group_by(point_id) |>
  mutate(
    ndvi_imputed = ndvi_mean,
    ndvi_imputed = if_else(
      ndvi_flag == "MISSING",
      lag(ndvi_mean, 1, default = mean(ndvi_mean, na.rm = TRUE)),
      ndvi_imputed
    ),
    imputed = (ndvi_flag == "MISSING")
  ) |>
  ungroup()

pct_missing <- mean(ndvi_monthly_imp$imputed) * 100
cat(sprintf("   Overall NA rate: %.1f%%  (will refine with cluster-neighbour imputation)\n\n",
            pct_missing))

write_csv(ndvi_monthly_imp, file.path(DATA_DIR, "ndvi_subcounty_monthly_imputed.csv"))


# =============================================================================
# SECTION 10: FDA CLUSTERING — IDENTIFY VEGETATION PHENOLOGY ZONES
#
# Goal: find k groups of sub-counties with similar seasonal NDVI functional shapes.
# These groups correspond to the livestock vs. crop vegetation zones that Abel
# described — differentiated by both amplitude AND timing of the annual NDVI cycle.
#
# Method: Functional Data Analysis (Ramsay & Silverman 2005)
#   1. Compute 12-month climatological NDVI per sub-county (from LOCF-imputed data)
#   2. Represent each profile as a functional data object using a Fourier basis
#      (natural basis for periodic/cyclic data; period = 12 months)
#   3. Apply functional PCA (fPCA) to decompose between-sub-county variation into
#      orthogonal modes (PC1 ≈ amplitude; PC2 ≈ phase/timing shifts)
#   4. Cluster sub-counties with k-means on fPCA scores (k = 2, 3, 4 tested)
#   5. Select k by maximising average silhouette width
#
# Advantage over DTW-PAM: FDA explicitly models the continuous curve shape,
# captures both amplitude and timing variation, and produces smooth functional
# mean curves (not a single observed medoid) for each cluster.
#
# Reference: Ramsay & Silverman (2005) Functional Data Analysis, Springer.
# =============================================================================

cat("── Section 10: FDA clustering of NDVI seasonal profiles...\n")

# Step 1: Per-sub-county 12-month climatological NDVI (from LOCF-imputed data)
sub_clim <- ndvi_monthly_imp |>
  group_by(point_id, month_n) |>
  summarise(ndvi_clim = mean(ndvi_imputed, na.rm = TRUE), .groups = "drop")

sub_clim_wide <- sub_clim |>
  pivot_wider(names_from = point_id, values_from = ndvi_clim) |>
  arrange(month_n) |>
  select(-month_n)

# Matrix: rows = months 1:12, columns = sub-county point IDs
clim_mat      <- as.matrix(sub_clim_wide)   # 12 × n_sub
pts_ordered   <- colnames(clim_mat)

cat(sprintf("   Climatology matrix: %d months × %d sub-counties\n",
            nrow(clim_mat), ncol(clim_mat)))

# Replace any residual NAs (isolated points with all-NA months) with column mean
for (j in seq_len(ncol(clim_mat))) {
  na_idx <- is.na(clim_mat[, j])
  if (any(na_idx)) clim_mat[na_idx, j] <- mean(clim_mat[, j], na.rm = TRUE)
}

# Step 2: Fourier basis — natural for cyclic annual data (period = 12 months)
# nbasis = 5 → 1 constant + 2 sin/cos pairs (fundamental + 1st harmonic)
# Captures Uganda's bimodal rainy season pattern without overfitting
n_basis <- 5
fbasis  <- fda::create.fourier.basis(rangeval = c(1, 12), nbasis = n_basis)
fd_smooth <- fda::smooth.basis(
  argvals  = 1:12,
  y        = clim_mat,
  fdParobj = fbasis
)$fd

cat(sprintf("   Fourier basis: %d basis functions | period = 12 months\n", n_basis))

# Step 3: Functional PCA — extract dominant modes of seasonal variation
n_harm   <- min(3, ncol(clim_mat) - 1)
fpca_out <- fda::pca.fd(fd_smooth, nharm = n_harm)

pve     <- fpca_out$varprop * 100
cum_pve <- cumsum(pve)
cat("   fPCA variance explained per harmonic:\n")
for (h in seq_along(pve))
  cat(sprintf("     FPC%d: %.1f%%  (cumulative: %.1f%%)\n", h, pve[h], cum_pve[h]))

# fPCA scores: n_sub × n_harm
fpc_scores <- fpca_out$scores
rownames(fpc_scores) <- pts_ordered

# fPCA variance-explained bar plot — saved to FDA_OUT_DIR
pve_df <- data.frame(
  Component = paste0("FPC", seq_along(pve)),
  VarExp    = pve,
  CumVarExp = cum_pve
)
p_fpca_var <- ggplot(pve_df, aes(Component, VarExp)) +
  geom_col(fill = "#27ae60", alpha = 0.85, width = 0.6) +
  geom_line(aes(y = CumVarExp, group = 1), colour = "#e74c3c",
            linewidth = 1.2) +
  geom_point(aes(y = CumVarExp), colour = "#e74c3c", size = 3) +
  geom_text(aes(y = VarExp + 1.5,
                label = sprintf("%.1f%%", VarExp)),
            size = 4) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 20)) +
  labs(
    title    = "FDA fPCA — Variance Explained by Functional Principal Components",
    subtitle = "Fourier basis (5 fns) | Isingiro District sub-county NDVI climatology",
    x = NULL, y = "Variance explained (%)",
    caption  = "Red line = cumulative variance"
  ) +
  theme_minimal(base_size = 13)
ggsave(file.path(FDA_OUT_DIR, "fda_fpca_variance.png"),
       p_fpca_var, width = 8, height = 6, dpi = 150)
cat("   Saved: fda_fpca_variance.png\n")

# Step 4: K-means on fPCA scores (k = 2, 3, 4)
cat("   Testing k = 2, 3, 4 clusters (silhouette on fPCA scores):\n")

best_k   <- 2
best_sil <- -Inf
cluster_results_fda <- list()

for (k in 2:min(4, nrow(fpc_scores) - 1)) {
  set.seed(2026)
  km  <- kmeans(fpc_scores, centers = k, nstart = 25, iter.max = 100)
  sil <- tryCatch({
    sil_obj <- cluster::silhouette(km$cluster, dist(fpc_scores))
    mean(sil_obj[, "sil_width"])
  }, error = function(e) NA_real_)
  cat(sprintf("     k = %d: silhouette = %.3f\n", k, sil))
  cluster_results_fda[[as.character(k)]] <- km
  if (!is.na(sil) && sil > best_sil) {
    best_sil <- sil
    best_k   <- k
  }
}

cat(sprintf("   Best k = %d (highest silhouette = %.3f)\n\n", best_k, best_sil))

best_km <- cluster_results_fda[[as.character(best_k)]]

cluster_assignments <- data.frame(
  point_id    = pts_ordered,
  fda_cluster = best_km$cluster,
  stringsAsFactors = FALSE
) |>
  left_join(SAMPLE_POINTS, by = c("point_id" = "id"))

write_csv(cluster_assignments, file.path(FDA_OUT_DIR, "fda_cluster_assignments.csv"))
cat("   Saved: fda_cluster_assignments.csv\n\n")

# Cluster summary
cluster_summary <- cluster_assignments |>
  count(fda_cluster) |>
  mutate(
    lat_mean = sapply(fda_cluster, function(cl)
                 mean(cluster_assignments$latitude[cluster_assignments$fda_cluster == cl])),
    lon_mean = sapply(fda_cluster, function(cl)
                 mean(cluster_assignments$longitude[cluster_assignments$fda_cluster == cl]))
  )
cat("   FDA cluster summary:\n")
print(cluster_summary)
cat("\n")


# ── Update NA imputation with within-cluster calendar-month mean ─────────────
cat("── Section 10b: Refining NA imputation using cluster neighbours...\n")

# For each missing value, substitute the mean NDVI of all OTHER cluster members
# for the same (year, month_n)
ndvi_with_cluster <- ndvi_monthly_imp |>
  left_join(cluster_assignments |> select(point_id, fda_cluster),
            by = "point_id")

cluster_month_means <- ndvi_with_cluster |>
  filter(!imputed) |>
  group_by(fda_cluster, year, month_n) |>
  summarise(cluster_ndvi_mean = mean(ndvi_mean, na.rm = TRUE),
            .groups = "drop")

ndvi_monthly_imp <- ndvi_with_cluster |>
  left_join(cluster_month_means, by = c("fda_cluster", "year", "month_n")) |>
  mutate(
    ndvi_imputed = if_else(imputed & !is.na(cluster_ndvi_mean),
                           cluster_ndvi_mean,
                           ndvi_imputed),
    imputation_method = case_when(
      !imputed                         ~ "observed",
      !is.na(cluster_ndvi_mean)        ~ "cluster_neighbour",
      TRUE                             ~ "LOCF"
    )
  ) |>
  select(-cluster_ndvi_mean)

write_csv(ndvi_monthly_imp, file.path(DATA_DIR, "ndvi_subcounty_monthly_imputed.csv"))

imp_tbl <- ndvi_monthly_imp |> count(imputation_method)
cat("   Final imputation summary:\n")
print(imp_tbl)
cat("\n")


# =============================================================================
# SECTION 11: CYCLIC NDVI CLIMATOLOGY (12-MONTH MEANS PER FDA CLUSTER)
#
# Abel's guidance: do not use a single mean across years. Compute a 12-month
# climatological mean (one value per calendar month) per cluster. This gives
# a stable cyclic seasonal forcing function for the long ODE runs.
# =============================================================================

cat("── Section 11: Computing 12-month cyclic climatology per FDA cluster...\n")

ndvi_climatology <- ndvi_monthly_imp |>
  group_by(fda_cluster, month_n) |>
  summarise(
    ndvi_clim_mean = mean(ndvi_imputed, na.rm = TRUE),
    ndvi_clim_sd   = sd(ndvi_imputed,   na.rm = TRUE),
    n_obs          = sum(!imputed),
    .groups        = "drop"
  ) |>
  mutate(
    month_abbr = factor(month.abb[month_n], levels = month.abb),
    cluster_label = sprintf("Cluster %d", fda_cluster)
  )

write_csv(ndvi_climatology, file.path(DATA_DIR, "ndvi_climatology_monthly.csv"))
cat(sprintf("   Saved: ndvi_climatology_monthly.csv  (%d rows)\n\n",
            nrow(ndvi_climatology)))


# =============================================================================
# SECTION 12: PER-CLUSTER ENVIRONMENTAL LAMBDA MULTIPLIER
#
# lambda_env(cluster, month) = 0.2 + 0.8 × normalised_NDVI_clim
#
# Livestock-zone clusters (higher, more variable NDVI) will have a higher
# peak lambda_env than crop-zone clusters — reflecting elevated vector density
# in areas with short shrubs, flooding, and livestock-related habitat.
# =============================================================================

cat("── Section 12: Computing per-cluster lambda_env multiplier...\n")

lambda_env_cluster <- ndvi_climatology |>
  group_by(fda_cluster) |>
  mutate(
    ndvi_norm  = (ndvi_clim_mean - min(ndvi_clim_mean, na.rm = TRUE)) /
                 (max(ndvi_clim_mean, na.rm = TRUE) -
                  min(ndvi_clim_mean, na.rm = TRUE)),
    lambda_env = 0.2 + 0.8 * ndvi_norm
  ) |>
  ungroup()

write_csv(lambda_env_cluster, file.path(FDA_OUT_DIR, "lambda_env_fda_cluster.csv"))
cat(sprintf("   Saved: lambda_env_fda_cluster.csv\n"))

lambda_env_cluster |>
  group_by(fda_cluster, cluster_label) |>
  summarise(min_lam = min(lambda_env), max_lam = max(lambda_env),
            mean_lam = mean(lambda_env), .groups = "drop") |>
  print()
cat("\n")


# =============================================================================
# SECTION 13: PLOTS
# =============================================================================

cat("── Section 13: Generating plots...\n")

# ── Cluster spatial map ────────────────────────────────────────────────────────
cluster_sf <- cluster_assignments |>
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

p_cluster_map <- ggplot() +
  geom_sf(data = isingiro_poly, fill = "grey95", colour = "grey50",
          linewidth = 0.8) +
  geom_sf(data = cluster_sf,
          aes(colour = factor(fda_cluster), shape = factor(fda_cluster)),
          size = 4, alpha = 0.9) +
  geom_text(data = cluster_assignments,
            aes(longitude, latitude, label = subcounty),
            size = 2.2, vjust = -0.8, colour = "grey20") +
  scale_colour_viridis_d(option = "D", name = "FDA Cluster",
                          labels = sprintf("Cluster %d", sort(unique(cluster_assignments$fda_cluster)))) +
  scale_shape_manual(values = c(16, 17, 15, 18)[1:best_k],
                     name = "FDA Cluster",
                     labels = sprintf("Cluster %d", sort(unique(cluster_assignments$fda_cluster)))) +
  labs(
    title    = "FDA Phenology Clusters by Sub-County — Isingiro District, Uganda",
    subtitle = sprintf("k = %d clusters | %d sub-county centroids (ADM4) | %s to %s",
                       best_k, nrow(SAMPLE_POINTS), DATE_START, DATE_END),
    x = "Longitude", y = "Latitude",
    caption  = "FDA clusters (Fourier basis + fPCA + k-means) reflect distinct seasonal NDVI phenology — amplitude and timing differentiation"
  ) +
  theme_minimal(base_size = 12) +
  coord_sf()

ggsave(file.path(FDA_OUT_DIR, "fda_cluster_map.png"),
       p_cluster_map, width = 10, height = 9, dpi = 150)
cat("   Saved: fda_cluster_map.png\n")

# ── Cluster NDVI profiles ──────────────────────────────────────────────────────
p_profiles <- ggplot(ndvi_climatology,
                     aes(month_abbr, ndvi_clim_mean,
                         colour = cluster_label, group = cluster_label)) +
  geom_ribbon(aes(ymin = ndvi_clim_mean - ndvi_clim_sd,
                  ymax = ndvi_clim_mean + ndvi_clim_sd,
                  fill = cluster_label),
              alpha = 0.15, colour = NA) +
  geom_line(linewidth = 1.3) +
  geom_point(size = 2.5) +
  scale_colour_viridis_d(option = "D", name = "Cluster") +
  scale_fill_viridis_d(option = "D", name = "Cluster") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    title    = "Mean Functional NDVI Profiles by FDA Cluster",
    subtitle = "Functional mean curves (Fourier basis) | 12-month climatology ± SD | Isingiro District",
    x = "Month", y = "Mean NDVI (climatological)"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(file.path(FDA_OUT_DIR, "ndvi_cluster_profiles.png"),
       p_profiles, width = 12, height = 7, dpi = 150)
cat("   Saved: ndvi_cluster_profiles.png\n")

# ── Lambda_env per cluster ────────────────────────────────────────────────────
p_lambda <- ggplot(lambda_env_cluster,
                   aes(month_abbr, lambda_env,
                       colour = cluster_label, group = cluster_label)) +
  geom_hline(yintercept = c(0.2, 1.0), linetype = "dashed",
             colour = "grey70", linewidth = 0.6) +
  geom_line(linewidth = 1.3) +
  geom_point(size = 2.5) +
  scale_colour_viridis_d(option = "D", name = "Cluster") +
  scale_y_continuous(limits = c(0, 1.15), breaks = seq(0.2, 1.0, 0.2)) +
  annotate("text", x = 0.7, y = 1.04, label = "Peak (1.0)",
           hjust = 0, colour = "grey50", size = 3.5) +
  annotate("text", x = 0.7, y = 0.14, label = "Minimum (0.2)",
           hjust = 0, colour = "grey50", size = 3.5) +
  labs(
    title    = "Seasonal λ_env Multiplier per FDA Cluster — Isingiro District",
    subtitle = "Higher peak = more vector-favourable seasonal conditions (expected in livestock clusters)",
    x = "Month", y = "λ_env  [0.2 – 1.0]"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom")

ggsave(file.path(FDA_OUT_DIR, "lambda_env_plot.png"),
       p_lambda, width = 12, height = 7, dpi = 150)
cat("   Saved: lambda_env_plot.png\n")

# ── Full NDVI time series (all sub-counties, coloured by cluster) ────────────────
ndvi_monthly_plot <- ndvi_monthly_imp |>
  mutate(cluster_label = sprintf("Cluster %d", fda_cluster))

p_ts <- ggplot(ndvi_monthly_plot,
               aes(month, ndvi_imputed, group = point_id,
                   colour = cluster_label)) +
  geom_line(linewidth = 0.5, alpha = 0.6) +
  scale_colour_viridis_d(option = "D", name = "Cluster") +
  scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
  facet_wrap(~cluster_label, ncol = 1) +
  labs(
    title    = "NDVI Time Series by FDA Cluster — Isingiro District",
    subtitle = sprintf("All %d sub-counties | MODIS MOD13Q1.061 | %s to %s",
                       nrow(SAMPLE_POINTS), DATE_START, DATE_END),
    x = NULL, y = "NDVI (quality-masked, imputed)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        strip.text = element_text(face = "bold"))

ggsave(file.path(FDA_OUT_DIR, "ndvi_timeseries_plot.png"),
       p_ts, width = 16, height = 4 * best_k, dpi = 150)
cat("   Saved: ndvi_timeseries_plot.png\n\n")


# =============================================================================
# SECTION 14: HOW TO USE CLUSTER lambda_env IN THE ODE MODEL
# =============================================================================

cat("── Section 14: How to use FDA cluster lambda_env in the ODE models ──────────\n")
cat("
  Load the FDA cluster forcing data:
    lambda_data <- read_csv('data/ndvi_fda_outputs/lambda_env_fda_cluster.csv')

  Build one interpolation function per FDA cluster:
    make_lambda_fn <- function(cl_data, lambda_base) {
      approxfun(
        x      = c(0, cl_data$month_n, 13),
        y      = c(cl_data$lambda_env[12],
                   cl_data$lambda_env,
                   cl_data$lambda_env[1]) * lambda_base,
        method = 'linear', rule = 2
      )
    }

    lambda_fns <- lapply(
      split(lambda_data, lambda_data$fda_cluster),
      make_lambda_fn, lambda_base = lambda_base_central
    )

  In the ODE, select the FDA cluster for the sub-county you are modelling:
    lambda_t <- lambda_fns[['1']](t_to_month_n(t))

  Run Model 2 separately for each FDA cluster to compare:
    - Livestock-zone cluster: higher lambda_env → more infections, more losses
    - Crop-zone cluster: lower lambda_env → fewer infections, fewer losses
    FDA clustering captures both amplitude and timing (phenological) differences,
    providing better differentiation between seasonal vegetation zones.
\n")


# =============================================================================
# SECTION 15: OUTPUT MANIFEST
# =============================================================================

cat("── Output manifest ──────────────────────────────────────────────────────\n")

outputs <- c(
  file.path(DATA_DIR,     "isingiro_subcounty_points.csv"),
  file.path(DATA_DIR,     "isingiro_subcounty_map.png"),
  file.path(DATA_DIR,     "ndvi_subcounty_16day_raw.csv"),
  file.path(DATA_DIR,     "ndvi_subcounty_monthly.csv"),
  file.path(DATA_DIR,     "ndvi_subcounty_monthly_imputed.csv"),
  file.path(DATA_DIR,     "ndvi_climatology_monthly.csv"),
  file.path(FDA_OUT_DIR,  "fda_cluster_assignments.csv"),
  file.path(FDA_OUT_DIR,  "fda_cluster_map.png"),
  file.path(FDA_OUT_DIR,  "ndvi_cluster_profiles.png"),
  file.path(FDA_OUT_DIR,  "lambda_env_fda_cluster.csv"),
  file.path(FDA_OUT_DIR,  "lambda_env_plot.png"),
  file.path(FDA_OUT_DIR,  "ndvi_timeseries_plot.png"),
  file.path(FDA_OUT_DIR,  "fda_fpca_variance.png")
)

for (f in outputs) {
  status <- if (file.exists(f)) "✓" else "✗ MISSING"
  size   <- if (file.exists(f)) sprintf("(%.1f KB)", file.size(f)/1024) else ""
  cat(sprintf("   %s  %-50s %s\n", status, basename(f), size))
}

cat("\n── Complete. Next steps:\n")
cat("   1. Open fda_cluster_map.png and fda_fpca_variance.png — verify clusters align with livestock/crop zones\n")
cat("   2. If k is wrong, change best_k and re-run from Section 10 onwards\n")
cat("   3. Run RVF_Models_Isingiro_Uganda_v1.3_FDA_17May2026.Rmd — loads lambda_env_fda_cluster.csv\n")
