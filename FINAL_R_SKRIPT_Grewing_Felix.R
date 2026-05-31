# FINAL_R_SKRIPT_Grewing_Felix.R
#
# Author:      Felix Grewing, University of Oldenburg
# Description: Replication script for OLS analysis of EIU Democracy Index dimensions
#              and international life satisfaction (cross-section, N = 102 countries).
# R version:   4.5.1
#
# Usage: Set BASE_DIR to the folder containing all data files, then run top to bottom.
#        Results are written to out_analysis_<timestamp>/out/.
#        On first run, EXPECTED_N_ANALYSIS_SAMPLE.txt is created to lock the analysis N.

options(stringsAsFactors = FALSE)
set.seed(123)

# USER SETTINGS ----
BASE_DIR <- file.path(path.expand("~"), "Desktop", "Datenteile Data")  # <- adjust to your data folder
MASTER_DEFAULT <- "Daten.xlsx"
SKIP_MAIN <- 1

SS_DEFAULT  <- "social_support_whr2025.xlsx"
QOG_DEFAULT <- "qog_std_ts_jan26.csv"

# V-Dem: set path below; prompts file.choose() if not found
USE_VDEM <- TRUE
VDEM_DEFAULT <- "V-Dem-CY-Core-v15.xlsx"
VDEM_YEAR_MAX <- 2024        # per country: most recent year <= VDEM_YEAR_MAX
VDEM_MODE <- "latest"        # "latest" | "align_to_whr" | "period_2010_2015"
VDEM_ALIGN_YEARS <- 2021:2023

# Sample size lock (EIU core sample):
STOP_ON_N_MISMATCH <- TRUE
EXPECTED_N_FILE <- file.path(BASE_DIR, "EXPECTED_N_ANALYSIS_SAMPLE.txt")
# To change N intentionally: set STOP_ON_N_MISMATCH <- FALSE or delete EXPECTED_N_FILE.

# Bootstrap (optional, slow):
DO_BOOTSTRAP <- TRUE
B_BOOT <- 1000

# Packages ----
req <- c(
  "readxl","readr","dplyr","tidyr","tibble","purrr","stringr","stringi","janitor",
  "countrycode","ggplot2","splines","MASS","broom",
  "lmtest","sandwich","car",
  "officer","flextable","modelsummary",
  "corrplot","reshape2",
  "glmnet","Matrix"
)

to_install <- setdiff(req, rownames(installed.packages()))
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)
invisible(lapply(req, library, character.only = TRUE))

# Optional (World map)
HAS_SF <- TRUE
if (!requireNamespace("sf", quietly = TRUE) ||
    !requireNamespace("rnaturalearth", quietly = TRUE) ||
    !requireNamespace("rnaturalearthdata", quietly = TRUE)) {
  HAS_SF <- FALSE
} else {
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
}

# WD + Output folders ----
if (!dir.exists(BASE_DIR)) {
  message("BASE_DIR not found. Using getwd(): ", getwd())
  BASE_DIR <- getwd()
}
setwd(BASE_DIR)

RUN_ID <- format(Sys.time(), "%Y%m%d_%H%M")
PATH_ROOT <- file.path(getwd(), paste0("out_analysis_", RUN_ID))
PATH_OUT  <- file.path(PATH_ROOT, "out")
FIG_DIR   <- file.path(PATH_OUT, "figs")
dir.create(PATH_OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR,  showWarnings = FALSE, recursive = TRUE)

save_plot <- function(p, filename, w=9.5, h=5.8, dpi=300){
  ggplot2::ggsave(filename = file.path(FIG_DIR, filename),
                  plot = p, width = w, height = h, dpi = dpi)
}

write_lines <- function(x, fn) writeLines(x, con = file.path(PATH_OUT, fn))

# Robust VCOV helpers
vc_hc3 <- function(m) sandwich::vcovHC(m, type = "HC3")

vc_cluster_from_data <- function(m, data, cluster_var){
  mf <- model.frame(m)
  idx <- suppressWarnings(as.integer(rownames(mf)))  # requires numeric rownames 1..n
  if (anyNA(idx)) stop("Cluster-VCOV: rownames(model.frame) are not numeric. Check rownames(data).")
  cl <- data[[cluster_var]][idx]
  sandwich::vcovCL(m, cluster = cl, type = "HC3")
}

mk_formula <- function(response, terms){
  terms <- terms[!is.na(terms) & terms != ""]
  if (length(terms) == 0) as.formula(paste(response, "~ 1")) else reformulate(terms, response = response)
}

find_col <- function(patterns, nms){
  hit <- nms[stringr::str_detect(nms, paste(patterns, collapse="|"))]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

# ISO3 mapping (robust) ----
norm_key <- function(x){
  x <- as.character(x)
  x <- stringi::stri_trans_general(x, "Latin-ASCII")
  x <- tolower(x)
  x <- gsub("[^a-z ]", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

apply_custom_iso <- function(cn, iso_vec){
  custom <- c(
    "united states"="USA","united states of america"="USA","usa"="USA",
    "united kingdom"="GBR","russia"="RUS","iran"="IRN","venezuela"="VEN","bolivia"="BOL",
    "cote d ivoire"="CIV","ivory coast"="CIV",
    "democratic republic of the congo"="COD","congo dem rep"="COD","dr congo"="COD",
    "republic of the congo"="COG","congo"="COG",
    "korea south"="KOR","south korea"="KOR","korea rep"="KOR","korea, rep"="KOR",
    "korea north"="PRK","north korea"="PRK",
    "tanzania"="TZA","czechia"="CZE","czech republic"="CZE",
    "egypt"="EGY","syrian arab republic"="SYR","syria"="SYR",
    "new zealand"="NZL"
  )
  k <- norm_key(cn)
  fix <- unname(custom[k])
  iso_vec[is.na(iso_vec) & !is.na(fix)] <- fix[is.na(iso_vec) & !is.na(fix)]
  iso_vec
}

best_iso3 <- function(cn){
  cn <- as.character(cn)

  try_map <- function(origin){
    out <- tryCatch(countrycode::countrycode(cn, origin = origin, destination = "iso3c", warn = FALSE),
                    error = function(e) rep(NA_character_, length(cn)))
    out
  }

  c1 <- try_map("country.name")
  c2 <- try_map("country.name.de")
  n1 <- sum(!is.na(c1))
  n2 <- sum(!is.na(c2))
  iso <- if (n1 >= n2) c1 else c2
  origin_used <- if (n1 >= n2) "country.name" else "country.name.de"

  iso <- apply_custom_iso(cn, iso)
  attr(iso, "origin_used") <- origin_used
  iso
}

# Load Master data ----
main_file <- MASTER_DEFAULT
if (!file.exists(main_file)) {
  message("Master file not found. Please select (xlsx/csv).")
  main_file <- file.choose()
}
is_csv  <- grepl("\\.csv$",  main_file, ignore.case = TRUE)
is_xlsx <- grepl("\\.xlsx?$", main_file, ignore.case = TRUE)
if (!is_csv && !is_xlsx) stop("Please select a .xlsx or .csv file as master table.")

raw <- if (is_xlsx) readxl::read_xlsx(main_file, skip = SKIP_MAIN) else readr::read_csv(main_file, show_col_types = FALSE)
raw <- raw |> janitor::clean_names() |> as.data.frame()
nms <- names(raw)

col_country    <- find_col(c("^country$","^land$"), nms)
col_eiu_overall<- find_col(c("eiu_overall_score","^eiu_overall$","eiu.*overall"), nms)
col_life       <- find_col(c("ladder_score_life_satisfaction","life_satisfaction","ladder"), nms)
must <- c(country=col_country, eiu_overall=col_eiu_overall, life_satisfaction=col_life)
if (any(is.na(must))) {
  print(must)
  stop("Required columns not found. Check column names or extend find_col patterns.")
}

col_dim_ep <- find_col(c("electoral_process_and_pluralism","electoral.*pluralism"), nms)
col_dim_fg <- find_col(c("functioning_of_government","functioning.*government"), nms)
col_dim_pp <- find_col(c("political_participation","political.*participation"), nms)
col_dim_pc <- find_col(c("political_culture","political.*culture"), nms)
col_dim_cl <- find_col(c("civil_liberties","civil.*libert"), nms)

col_gdp     <- find_col(c("gdp_per_capita_ppp","gdp.*ppp","gdp.*constant"), nms)
col_unemp   <- find_col(c("unemployment"), nms)
col_old     <- find_col(c("age_dependency_ratio_old","dependency.*old","old_age"), nms)
col_socprot <- find_col(c("total_expenditure_on_social_protection","social_protection","expenditure.*social"), nms)
col_region  <- find_col(c("^region$"), nms)
col_income  <- find_col(c("^income_group$","income.*group"), nms)

dat_mod <- tibble::tibble(
  country = as.character(raw[[col_country]]),
  life_satisfaction = suppressWarnings(as.numeric(raw[[col_life]])),
  eiu_overall = suppressWarnings(as.numeric(raw[[col_eiu_overall]]))
)

dims_map <- c(
  electoral_process_pluralism = col_dim_ep,
  functioning_government      = col_dim_fg,
  political_participation     = col_dim_pp,
  political_culture           = col_dim_pc,
  civil_liberties             = col_dim_cl
)
for (nm in names(dims_map)) {
  if (!is.na(dims_map[[nm]])) dat_mod[[nm]] <- suppressWarnings(as.numeric(raw[[dims_map[[nm]]]]))
}

if (!is.na(col_gdp))     dat_mod$gdp_pc_ppp <- suppressWarnings(as.numeric(raw[[col_gdp]]))
if (!is.na(col_unemp))   dat_mod$unemployment_rate <- suppressWarnings(as.numeric(raw[[col_unemp]]))
if (!is.na(col_old))     dat_mod$old_age_dep <- suppressWarnings(as.numeric(raw[[col_old]]))
if (!is.na(col_socprot)) dat_mod$social_protection_tot <- suppressWarnings(as.numeric(raw[[col_socprot]]))
if (!is.na(col_region))  dat_mod$region <- as.character(raw[[col_region]])
if (!is.na(col_income))  dat_mod$income_group <- as.character(raw[[col_income]])

dat_mod$iso3c <- best_iso3(dat_mod$country)
dat_mod$ln_gdp_pc_ppp <- if ("gdp_pc_ppp" %in% names(dat_mod)) {
  ifelse(is.finite(dat_mod$gdp_pc_ppp) & dat_mod$gdp_pc_ppp > 0, log(dat_mod$gdp_pc_ppp), NA_real_)
} else NA_real_

# unemployment scaling fix (robust)
if ("unemployment_rate" %in% names(dat_mod)) {
  dat_mod$unemployment_rate <- dplyr::case_when(
    is.na(dat_mod$unemployment_rate) ~ NA_real_,
    dat_mod$unemployment_rate == 0 ~ NA_real_,
    dat_mod$unemployment_rate > 100 ~ dat_mod$unemployment_rate / 1000,
    TRUE ~ dat_mod$unemployment_rate
  )
}

# WHR social support: load and merge ----
read_whr_social_support <- function(ss_file, iso3_ref = NULL){
  # robust: tries different skips AND enforces that the value column is numeric
  attempts <- list(list(skip=0, tag="skip0"),
                   list(skip=1, tag="skip1"),
                   list(skip=2, tag="skip2"))
  best <- NULL
  best_score <- -Inf
  best_meta <- list()

  for (a in attempts){
    tmp <- tryCatch(readxl::read_xlsx(ss_file, skip = a$skip, .name_repair = "unique"),
                    error = function(e) NULL)
    if (is.null(tmp) || ncol(tmp) < 2) next
    tmp <- tmp |> janitor::clean_names() |> as.data.frame()
    nms <- names(tmp)

    # choose country column (prefer 'land'/'country' names; fallback to first)
    col_country <- find_col(c("^land$", "^country$", "country_name", "country"), nms)
    if (is.na(col_country)) col_country <- nms[1]

    # choose value column among remaining columns: maximize numeric coverage
    value_candidates <- setdiff(nms, col_country)
    if (length(value_candidates) == 0) next

    # prioritize columns that look like social support
    pri <- value_candidates[grepl("support", value_candidates, ignore.case = TRUE)]
    ord <- unique(c(pri, value_candidates))

    best_val <- NULL
    best_num_n <- -Inf
    best_num_share_raw <- 0

    for (vcol in ord){
      x <- suppressWarnings(as.numeric(tmp[[vcol]]))
      num_n <- sum(!is.na(x))
      raw_n <- sum(!is.na(tmp[[vcol]]))
      num_share_raw <- ifelse(raw_n > 0, num_n / raw_n, 0)

      if (num_n > best_num_n){
        best_num_n <- num_n
        best_val <- vcol
        best_num_share_raw <- num_share_raw
      }
    }
    if (is.null(best_val)) next

    cand <- tmp |>
      dplyr::transmute(
        country_raw = as.character(.data[[col_country]]),
        social_support_01 = suppressWarnings(as.numeric(.data[[best_val]]))
      ) |>
      dplyr::filter(!is.na(.data$country_raw) & .data$country_raw != "") |>
      dplyr::filter(!tolower(.data$country_raw) %in% c("land","country"))

    cand$iso3c <- best_iso3(cand$country_raw)

    cand <- cand |>
      dplyr::filter(!is.na(.data$iso3c)) |>
      dplyr::group_by(.data$iso3c) |>
      dplyr::summarise(social_support_01 = dplyr::first(.data$social_support_01),
                       .groups = "drop")

    score_iso <- nrow(cand)
    score_hit <- if (!is.null(iso3_ref)) sum(cand$iso3c %in% iso3_ref) else 0L
    num_n <- sum(!is.na(cand$social_support_01))
    num_share <- ifelse(score_iso > 0, num_n / score_iso, 0)

    # Hard sanity: we need mostly numeric values; otherwise this skip is invalid
    if (num_share < 0.60){
      score <- -Inf
    } else {
      # prioritize numeric coverage heavily; then maximize overlap with analysis iso3_ref
      score <- 1000L * num_n + 10L * score_hit + score_iso
    }

    if (score > best_score){
      best <- cand
      best_score <- score
      best_meta <- list(skip=a$skip, tag=a$tag,
                        col_country=col_country, col_value=best_val,
                        score_iso=score_iso, score_hit=score_hit,
                        num_n=num_n, num_share=num_share,
                        num_share_raw=best_num_share_raw)
    }
  }

  if (is.null(best)) stop("WHR social_support: could not parse file (check header/columns).")

  attr(best, "meta") <- best_meta
  best
}

ss_file <- if (file.exists(SS_DEFAULT)) SS_DEFAULT else file.choose()
whr_ss <- read_whr_social_support(ss_file, iso3_ref = unique(na.omit(dat_mod$iso3c)))

meta <- attr(whr_ss, "meta")
write_lines(c(
  "=== WHR READ LOG ===",
  paste0("file: ", ss_file),
  paste0("skip_used: ", meta$skip, " (", meta$tag, ")"),
  paste0("col_country: ", meta$col_country),
  paste0("col_value: ", meta$col_value),
  paste0("score_iso: ", meta$score_iso, "  score_hit: ", meta$score_hit),
  paste0("num_n: ", meta$num_n, "  num_share: ", round(meta$num_share, 3), "  num_share_raw: ", round(meta$num_share_raw, 3)),
  paste0("rows (unique iso3): ", nrow(whr_ss))
), "whr_read_log.txt")

dat_mod <- dat_mod |> dplyr::left_join(whr_ss, by="iso3c")

# QoG/WGI6 mean 2022-2024 ----
qog_file <- if (file.exists(QOG_DEFAULT)) QOG_DEFAULT else file.choose()
needed <- c("ccodealp","year","wbgi_vae","wbgi_pve","wbgi_gee","wbgi_rqe","wbgi_rle","wbgi_cce")
qog <- readr::read_csv(qog_file, show_col_types = FALSE, col_select = dplyr::all_of(needed))

qog$ccodealp <- as.character(qog$ccodealp)
qog$year <- suppressWarnings(as.integer(qog$year))
wgi_cols <- c("wbgi_vae","wbgi_pve","wbgi_gee","wbgi_rqe","wbgi_rle","wbgi_cce")
for (cc in wgi_cols) qog[[cc]] <- suppressWarnings(as.numeric(qog[[cc]]))

qog_2224 <- qog |>
  dplyr::filter(year %in% 2022:2024) |>
  dplyr::group_by(ccodealp) |>
  dplyr::summarise(dplyr::across(dplyr::all_of(wgi_cols), ~ mean(.x, na.rm=TRUE)), .groups="drop") |>
  dplyr::mutate(qog_wgi6 = rowMeans(as.matrix(dplyr::pick(dplyr::all_of(wgi_cols))), na.rm=TRUE)) |>
  dplyr::rename(iso3c = ccodealp)

dat_mod <- dat_mod |> dplyr::left_join(qog_2224 |> dplyr::select(iso3c, qog_wgi6), by="iso3c")

# V-Dem: load and merge (optional) ----
maybe_rescale_vdem <- function(x){
  x <- suppressWarnings(as.numeric(readr::parse_number(as.character(x))))
  if (!any(is.finite(x))) return(x)
  mx <- max(x, na.rm = TRUE)
  # V-Dem values sometimes encoded as 0-1000 in Excel exports (e.g. 918 = 0.918)
  if (mx > 2) x <- x / 1000
  x
}

vdem_vars <- character(0)
vdem_raw_cache <- NULL
vdem_file_used <- NA_character_

read_vdem_raw <- function(vdem_file, iso3_ref = NULL){
  ext <- tolower(tools::file_ext(vdem_file))
  if (ext %in% c("csv","txt")) {
    tmp <- readr::read_csv(vdem_file, show_col_types = FALSE)
    tmp <- tmp |> janitor::clean_names() |> as.data.frame()
    attr(tmp, "meta") <- list(format="csv", skip_used=NA_integer_)
    return(tmp)
  }

  # xlsx: try skip 0..2 and pick best parse
  best <- NULL
  best_score <- -Inf
  best_meta <- list()

  for (sk in 0:2){
    tmp <- tryCatch(readxl::read_xlsx(vdem_file, skip = sk), error = function(e) NULL)
    if (is.null(tmp) || ncol(tmp) < 5) next
    tmp <- tmp |> janitor::clean_names() |> as.data.frame()
    nms <- names(tmp)

    col_iso  <- find_col(c("^country_text_id$","^country_textid$","^iso3c$","^iso3$"), nms)
    col_ctry <- find_col(c("^country_name$","^country$","country.*name"), nms)
    col_year <- find_col(c("^year$"), nms)

    # required V-Dem indices
    col_poly <- find_col(c("^v2x_polyarchy$"), nms)
    col_lib  <- find_col(c("^v2x_libdem$"), nms)
    col_part <- find_col(c("^v2x_partipdem$"), nms)
    col_del  <- find_col(c("^v2x_delibdem$"), nms)
    col_egal <- find_col(c("^v2x_egaldem$"), nms)

    # "score": counts available key columns + iso overlap
    score_cols <- sum(!is.na(c(col_year,col_poly,col_lib,col_part,col_del,col_egal)))
    # ISO3-coverage
    iso_vec <- NULL
    if (!is.na(col_iso)) iso_vec <- toupper(as.character(tmp[[col_iso]]))
    if (is.null(iso_vec) && !is.na(col_ctry)) {
      iso_vec <- best_iso3(as.character(tmp[[col_ctry]]))
    }
    hit <- if (!is.null(iso3_ref) && !is.null(iso_vec)) sum(iso_vec %in% iso3_ref, na.rm=TRUE) else 0L
    score <- 100L*score_cols + hit

    if (score > best_score){
      best <- tmp
      best_score <- score
      best_meta <- list(format="xlsx", skip_used=sk, score_cols=score_cols, score_hit=hit,
                        col_iso=col_iso, col_ctry=col_ctry, col_year=col_year,
                        col_poly=col_poly, col_lib=col_lib, col_part=col_part, col_del=col_del, col_egal=col_egal)
    }
  }

  if (is.null(best)) stop("V-Dem: could not read file (unexpected format).")

  attr(best, "meta") <- best_meta
  best
}

if (USE_VDEM) {

  vdem_file_used <- if (file.exists(VDEM_DEFAULT)) VDEM_DEFAULT else file.choose()
  vdem_raw_cache <- read_vdem_raw(vdem_file_used, iso3_ref = unique(na.omit(dat_mod$iso3c)))
  meta <- attr(vdem_raw_cache, "meta")

  # ISO3 + year robust
  nms <- names(vdem_raw_cache)
  getv <- function(df, nm) if (!is.na(nm) && nm %in% names(df)) df[[nm]] else rep(NA, nrow(df))

  iso3c <- NULL
  if (!is.null(meta$col_iso) && !is.na(meta$col_iso) && meta$col_iso %in% nms) {
    iso3c <- toupper(as.character(vdem_raw_cache[[meta$col_iso]]))
  } else if (!is.null(meta$col_ctry) && !is.na(meta$col_ctry) && meta$col_ctry %in% nms) {
    iso3c <- best_iso3(as.character(vdem_raw_cache[[meta$col_ctry]]))
  } else {
    # Fallback: suche irgendeine ISO-Spalte
    cand_iso <- find_col(c("^iso3c$","^iso3$","country_text_id","country_textid"), nms)
    if (!is.na(cand_iso)) iso3c <- toupper(as.character(vdem_raw_cache[[cand_iso]]))
  }

  year <- if (!is.null(meta$col_year) && !is.na(meta$col_year) && meta$col_year %in% nms) {
    suppressWarnings(as.integer(vdem_raw_cache[[meta$col_year]]))
  } else if ("year" %in% nms) {
    suppressWarnings(as.integer(vdem_raw_cache[["year"]]))
  } else {
    rep(NA_integer_, nrow(vdem_raw_cache))
  }

  # Core indices (create as NA if missing)
  v2x_polyarchy <- if (!is.null(meta$col_poly) && !is.na(meta$col_poly) && meta$col_poly %in% nms) vdem_raw_cache[[meta$col_poly]] else getv(vdem_raw_cache, "v2x_polyarchy")
  v2x_libdem    <- if (!is.null(meta$col_lib)  && !is.na(meta$col_lib)  && meta$col_lib  %in% nms) vdem_raw_cache[[meta$col_lib]]  else getv(vdem_raw_cache, "v2x_libdem")
  v2x_partipdem <- if (!is.null(meta$col_part) && !is.na(meta$col_part) && meta$col_part %in% nms) vdem_raw_cache[[meta$col_part]] else getv(vdem_raw_cache, "v2x_partipdem")
  v2x_delibdem  <- if (!is.null(meta$col_del)  && !is.na(meta$col_del)  && meta$col_del  %in% nms) vdem_raw_cache[[meta$col_del]]  else getv(vdem_raw_cache, "v2x_delibdem")
  v2x_egaldem   <- if (!is.null(meta$col_egal) && !is.na(meta$col_egal) && meta$col_egal %in% nms) vdem_raw_cache[[meta$col_egal]] else getv(vdem_raw_cache, "v2x_egaldem")

  vdem <- tibble::tibble(
    iso3c = iso3c,
    year  = year,
    v2x_polyarchy = maybe_rescale_vdem(suppressWarnings(as.numeric(v2x_polyarchy))),
    v2x_libdem    = maybe_rescale_vdem(suppressWarnings(as.numeric(v2x_libdem))),
    v2x_partipdem = maybe_rescale_vdem(suppressWarnings(as.numeric(v2x_partipdem))),
    v2x_delibdem  = maybe_rescale_vdem(suppressWarnings(as.numeric(v2x_delibdem))),
    v2x_egaldem   = maybe_rescale_vdem(suppressWarnings(as.numeric(v2x_egaldem)))
  ) |>
    dplyr::filter(!is.na(.data$iso3c), !is.na(.data$year))

  write_lines(c(
    "=== VDEM READ LOG ===",
    paste0("file: ", vdem_file_used),
    paste0("format: ", meta$format),
    paste0("skip_used: ", meta$skip_used),
    paste0("score_cols: ", meta$score_cols, "  score_hit: ", meta$score_hit),
    paste0("rows after iso3/year filter: ", nrow(vdem))
  ), "vdem_read_log.txt")

  vdem_mode_df <- function(mode = c("latest","period_2010_2015","align_to_whr"),
                           year_max = VDEM_YEAR_MAX,
                           align_years = VDEM_ALIGN_YEARS){
    mode <- match.arg(mode)
    base <- vdem |>
      dplyr::select(iso3c, year, v2x_polyarchy, v2x_libdem, v2x_partipdem, v2x_delibdem, v2x_egaldem)

    if (mode == "latest") {
      out <- base |>
        dplyr::filter(.data$year <= year_max) |>
        dplyr::group_by(.data$iso3c) |>
        dplyr::slice_max(.data$year, with_ties = FALSE) |>
        dplyr::ungroup() |>
        dplyr::select(-year)
    } else if (mode == "period_2010_2015") {   # 2010-2015 average, cf. Bromo et al. (2024)
      out <- base |>
        dplyr::filter(.data$year >= 2010, .data$year <= 2015) |>
        dplyr::group_by(.data$iso3c) |>
        dplyr::summarise(dplyr::across(dplyr::starts_with("v2x_"), ~ mean(.x, na.rm = TRUE)),
                         .groups = "drop")
    } else { # align_to_whr
      m <- base |>
        dplyr::filter(.data$year %in% align_years) |>
        dplyr::group_by(.data$iso3c) |>
        dplyr::summarise(dplyr::across(dplyr::starts_with("v2x_"), ~ mean(.x, na.rm = TRUE)),
                         .groups = "drop")
      fallback <- base |>
        dplyr::group_by(.data$iso3c) |>
        dplyr::slice_min(abs(.data$year - stats::median(align_years)), with_ties = FALSE) |>
        dplyr::ungroup() |>
        dplyr::select(-year)
      j <- dplyr::right_join(m, fallback, by="iso3c", suffix=c(".m",".fb"))
      for (nm in c("v2x_polyarchy","v2x_libdem","v2x_partipdem","v2x_delibdem","v2x_egaldem")) {
        j[[nm]] <- dplyr::coalesce(j[[paste0(nm, ".m")]], j[[paste0(nm, ".fb")]])
      }
      out <- j |> dplyr::select(iso3c, v2x_polyarchy, v2x_libdem, v2x_partipdem, v2x_delibdem, v2x_egaldem)
    }

    dplyr::rename(out,
      vdem_polyarchy   = v2x_polyarchy,
      vdem_liberal     = v2x_libdem,
      vdem_particip    = v2x_partipdem,
      vdem_deliber     = v2x_delibdem,
      vdem_egalitarian = v2x_egaldem
    )
  }

  # Merge V-Dem using selected mode
  vdem_df <- tryCatch(vdem_mode_df(mode = VDEM_MODE), error = function(e) NULL)
  if (!is.null(vdem_df) && nrow(vdem_df) > 0) {
    dat_mod <- dat_mod |> dplyr::left_join(vdem_df, by = "iso3c")
    vdem_vars <- intersect(c("vdem_polyarchy","vdem_liberal","vdem_particip","vdem_deliber","vdem_egalitarian"), names(dat_mod))
    # Fail-safe: if all vdem_* are NA after merge, skip V-Dem outputs
    if (all(dplyr::if_all(dat_mod[, vdem_vars, drop=FALSE], ~ is.na(.x)))) {
      warning("V-Dem merge: all vdem_* variables are NA. Check vdem_read_log.txt and ISO3 mapping.")
      vdem_vars <- character(0)
    }
  } else {
    warning("V-Dem could not be processed (vdem_df NULL/empty). Skipping V-Dem outputs.")
  }
}
# Analysis sample and controls ----
ctrl_candidates <- c("ln_gdp_pc_ppp","unemployment_rate","old_age_dep","social_protection_tot","qog_wgi6","social_support_01")
ctrl_vars <- intersect(ctrl_candidates, names(dat_mod))

eiu_dims <- c("electoral_process_pluralism","functioning_government",
              "political_participation","political_culture","civil_liberties")
eiu_dims <- intersect(eiu_dims, names(dat_mod))

# Enforce identical N across all EIU core models
need_common <- unique(c("country","iso3c","life_satisfaction","eiu_overall", eiu_dims, ctrl_vars))
if ("region" %in% names(dat_mod)) need_common <- unique(c(need_common, "region"))
if ("income_group" %in% names(dat_mod)) need_common <- unique(c(need_common, "income_group"))

dat_cc <- dat_mod |>
  dplyr::select(dplyr::all_of(need_common)) |>
  dplyr::filter(!is.na(.data$country) & !is.na(.data$iso3c)) |>
  dplyr::filter(dplyr::if_all(dplyr::all_of(setdiff(need_common, c("income_group","region"))), ~ !is.na(.x)))

# keep missing region/income as "Unknown" so FE/cluster does not drop cases
if ("region" %in% names(dat_cc)) dat_cc$region <- ifelse(is.na(dat_cc$region), "Unknown", dat_cc$region)
if ("income_group" %in% names(dat_cc)) dat_cc$income_group <- ifelse(is.na(dat_cc$income_group), "Unknown", dat_cc$income_group)

# Numeric rownames needed for cluster-VCOV
dat_cc <- as.data.frame(dat_cc)
dat_cc$.row <- seq_len(nrow(dat_cc))
rownames(dat_cc) <- dat_cc$.row

# Merge diagnostics
miss_iso3 <- dat_mod |> dplyr::filter(is.na(.data$iso3c) & !is.na(.data$country)) |> dplyr::distinct(country)
miss_whr  <- dat_mod |> dplyr::filter(is.na(.data$social_support_01) & !is.na(.data$iso3c)) |> dplyr::distinct(country, iso3c)
miss_qog  <- dat_mod |> dplyr::filter(is.na(.data$qog_wgi6) & !is.na(.data$iso3c)) |> dplyr::distinct(country, iso3c)
miss_vdem <- if (USE_VDEM && length(vdem_vars) > 0) {
  dat_mod |> dplyr::filter(!is.na(.data$iso3c)) |>
    dplyr::filter(dplyr::if_any(dplyr::all_of(vdem_vars), ~ is.na(.x))) |>
    dplyr::distinct(country, iso3c)
} else tibble::tibble()

readr::write_csv(miss_iso3, file.path(PATH_OUT, "check_missing_iso3.csv"))
readr::write_csv(miss_whr,  file.path(PATH_OUT, "check_missing_merge_whr.csv"))
readr::write_csv(miss_qog,  file.path(PATH_OUT, "check_missing_merge_qog.csv"))
if (USE_VDEM) readr::write_csv(miss_vdem, file.path(PATH_OUT, "check_missing_merge_vdem.csv"))

N_now <- nrow(dat_cc)

# Expected-N Check (EIU-Kernsample)
if (file.exists(EXPECTED_N_FILE)) {
  N_exp <- suppressWarnings(as.integer(readLines(EXPECTED_N_FILE)[1]))
  if (is.finite(N_exp) && !is.na(N_exp) && N_exp != N_now) {
    msg <- paste0("N mismatch: expected ", N_exp, " but got ", N_now,
                  ". Common cause: merge issue (WHR/QoG) or changed controls.")
    if (STOP_ON_N_MISMATCH) stop(msg) else warning(msg)
  }
} else {
  writeLines(as.character(N_now), EXPECTED_N_FILE)
}

qc <- tibble::tibble(
  metric = c("N_analysis_sample_EIU","controls_used",
             "missing_iso3","missing_whr","missing_qog","missing_vdem",
             "iso_origin_master","iso_origin_whr"),
  value = c(N_now,
            paste(ctrl_vars, collapse=", "),
            nrow(miss_iso3), nrow(miss_whr), nrow(miss_qog), nrow(miss_vdem),
            attr(dat_mod$iso3c, "origin_used"),
            attr(whr_ss$iso3c, "origin_used"))
)
readr::write_csv(qc, file.path(PATH_OUT, "00_qc_summary.csv"))

write_lines(c(
  paste0("N (analysis sample, EIU) = ", N_now),
  paste0("Controls used: ", paste(ctrl_vars, collapse=", ")),
  paste0("Missing ISO3: ", nrow(miss_iso3)),
  paste0("Missing WHR merges: ", nrow(miss_whr)),
  paste0("Missing QoG merges: ", nrow(miss_qog)),
  paste0("Missing V-Dem merges: ", nrow(miss_vdem))
), "checklist_report.txt")

# Table 2 (Summary stats) ----
desc_vars <- unique(c("life_satisfaction","eiu_overall", eiu_dims, ctrl_vars))
desc_vars <- desc_vars[desc_vars %in% names(dat_cc)]

desc <- dat_cc |>
  dplyr::summarise(dplyr::across(
    dplyr::all_of(desc_vars),
    list(N=~sum(!is.na(.x)), Mean=~mean(.x,na.rm=TRUE), SD=~sd(.x,na.rm=TRUE),
         Min=~min(.x,na.rm=TRUE), Max=~max(.x,na.rm=TRUE)),
    .names = "{.col}__{.fn}"
  )) |>
  tidyr::pivot_longer(everything(), names_to=c("var","stat"), names_sep="__", values_to="value") |>
  tidyr::pivot_wider(names_from=stat, values_from=value)

ft_desc <- flextable::flextable(desc) |>
  flextable::autofit() |>
  flextable::align(align="center", part="all") |>
  flextable::bold(part="header")

doc_desc <- officer::read_docx() |>
  officer::body_add_par("Table 2. Summary statistics (analysis sample)", style="heading 1") |>
  flextable::body_add_flextable(ft_desc)
print(doc_desc, target = file.path(PATH_OUT, "Table2_Summary_Stats.docx"))

# Table A1 (Country list) ----
app_countries <- dat_cc |> dplyr::distinct(country, iso3c) |> dplyr::arrange(country)
ftA1 <- flextable::flextable(app_countries) |>
  flextable::autofit() |>
  flextable::align(align="center", part="all") |>
  flextable::bold(part="header")

docA1 <- officer::read_docx() |>
  officer::body_add_par("Table A1. Country list (analysis sample)", style="heading 1") |>
  flextable::body_add_flextable(ftA1)
print(docA1, target = file.path(PATH_OUT, "TableA1_Country_List_EIU.docx"))

# Main models (Table 3) ----
m_overall_biv <- lm(life_satisfaction ~ eiu_overall, data = dat_cc)
m_overall_ctl <- lm(mk_formula("life_satisfaction", c("eiu_overall", ctrl_vars)), data = dat_cc)

dim_models <- list()
for (v in eiu_dims) dim_models[[v]] <- lm(mk_formula("life_satisfaction", c(v, ctrl_vars)), data = dat_cc)

models_main <- c(list("EIU overall (biv)"=m_overall_biv,
                      "EIU overall + ctrls"=m_overall_ctl),
                 setNames(dim_models, paste0("EIU dim: ", eiu_dims)))

# N overview
N_over <- tibble::tibble(
  model = names(models_main),
  N = purrr::map_int(models_main, nobs)
)
readr::write_csv(N_over, file.path(PATH_OUT, "00_model_N_overview.csv"))

ft_t3 <- modelsummary::modelsummary(
  models_main,
  vcov = function(m) vc_hc3(m),
  stars = TRUE,
  statistic = "({std.error})",
  output = "flextable"
)

doc_t3 <- officer::read_docx() |>
  officer::body_add_par("Table 3. EIU models (HC3 robust SE)", style="heading 1") |>
  flextable::body_add_flextable(ft_t3)
print(doc_t3, target = file.path(PATH_OUT, "Table3_EIU_dims_OLS.docx"))

# Simultaneous EIU main test: all 5 dimensions + controls ----
if (length(eiu_dims) == 5) {

  # Joint model: all five EIU dimensions + controls
  m_eiu_5 <- lm(mk_formula("life_satisfaction", c(eiu_dims, ctrl_vars)), data = dat_cc)

  ## Haupttabelle als DOCX ----
  ft_eiu5 <- modelsummary::modelsummary(
    list("EIU dimensions + controls" = m_eiu_5),
    vcov = function(m) vc_hc3(m),
    stars = TRUE,
    statistic = "({std.error})",
    output = "flextable"
  )

  doc_eiu5 <- officer::read_docx() |>
    officer::body_add_par("Table E1. EIU dimensions + controls (simultaneous model, HC3 robust SE)", style = "heading 1") |>
    flextable::body_add_flextable(ft_eiu5)

  print(doc_eiu5, target = file.path(PATH_OUT, "TableE1_EIU_5dims_OLS.docx"))

  ## Robuste Koeffizienten als CSV ----
  V_eiu5 <- vc_hc3(m_eiu_5)
  coefs_eiu5 <- broom::tidy(m_eiu_5) |>
    dplyr::mutate(
      std.error = sqrt(diag(V_eiu5)),
      statistic = estimate / std.error,
      p.value   = 2 * pt(abs(statistic), df = df.residual(m_eiu_5), lower.tail = FALSE),
      conf.low  = estimate - 1.96 * std.error,
      conf.high = estimate + 1.96 * std.error
    )

  readr::write_csv(coefs_eiu5, file.path(PATH_OUT, "TableE1_EIU_5dims_OLS_coefficients.csv"))

  ## Standardisierte Betas im simultanen Modell ----
  std_one_local <- function(x) as.numeric(scale(x))

  cols_joint <- unique(c("life_satisfaction", eiu_dims, ctrl_vars))
  d_joint <- dat_cc |>
    dplyr::select(dplyr::all_of(cols_joint)) |>
    tidyr::drop_na()

  d_joint_std <- d_joint |>
    dplyr::mutate(dplyr::across(where(is.numeric), std_one_local))

  m_eiu_5_std <- lm(mk_formula("life_satisfaction", c(eiu_dims, ctrl_vars)), data = d_joint_std)
  V_eiu5_std  <- vc_hc3(m_eiu_5_std)

  label_map_joint <- c(
    electoral_process_pluralism = "Electoral process & pluralism",
    functioning_government      = "Functioning of government",
    political_participation     = "Political participation",
    political_culture           = "Political culture",
    civil_liberties             = "Civil liberties"
  )

  std_joint_df <- tibble::tibble(
    term  = eiu_dims,
    beta  = as.numeric(coef(m_eiu_5_std)[eiu_dims]),
    se    = as.numeric(sqrt(diag(V_eiu5_std))[eiu_dims]),
    lo    = beta - 1.96 * se,
    hi    = beta + 1.96 * se,
    label = dplyr::recode(term, !!!label_map_joint),
    n     = nobs(m_eiu_5_std)
  )

  readr::write_csv(std_joint_df, file.path(PATH_OUT, "FigureE1_EIU_5dims_joint_std_coefs.csv"))

  ## Coefficient plot: simultaneous EIU main test ----
  p_eiu5 <- ggplot2::ggplot(std_joint_df, ggplot2::aes(x = beta, y = reorder(label, beta))) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = hi), height = 0.2) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Figure E1. Standardized effects of EIU dimensions (simultaneous model)",
      subtitle = "All five EIU dimensions entered jointly with the same controls (HC3 SE)",
      x = "Standardized coefficient (β)",
      y = ""
    )

  save_plot(p_eiu5, "figure_E1_EIU_5dims_joint_std_beta.png", w = 9.5, h = 5.8)

  ## Vergleich: separat vs simultan ----
  sep_joint_compare <- purrr::map_dfr(eiu_dims, function(term){
    cols <- unique(c("life_satisfaction", term, ctrl_vars))
    d <- dat_cc |>
      dplyr::select(dplyr::all_of(cols)) |>
      tidyr::drop_na() |>
      dplyr::mutate(dplyr::across(where(is.numeric), std_one_local))

    m_sep <- lm(mk_formula("life_satisfaction", c(term, ctrl_vars)), data = d)
    V_sep <- vc_hc3(m_sep)

    tibble::tibble(
      term = term,
      beta_sep = as.numeric(coef(m_sep)[term]),
      se_sep   = as.numeric(sqrt(diag(V_sep))[term])
    )
  })

  compare_eiu_df <- std_joint_df |>
    dplyr::select(term, label, beta_joint = beta, se_joint = se) |>
    dplyr::left_join(sep_joint_compare, by = "term") |>
    dplyr::mutate(delta_beta = beta_joint - beta_sep)

  readr::write_csv(compare_eiu_df, file.path(PATH_OUT, "TableE2_EIU_compare_separate_vs_joint.csv"))

  ## Additional diagnostics: joint EIU model ----
  vif_joint <- car::vif(m_eiu_5)

  if (is.matrix(vif_joint)) {
    vif_joint_df <- tibble::tibble(
      var  = rownames(vif_joint),
      GVIF = vif_joint[, 1],
      Df   = vif_joint[, 2],
      VIF  = vif_joint[, 3]
    )
  } else {
    vif_joint_df <- tibble::tibble(
      var = names(vif_joint),
      VIF = as.numeric(vif_joint)
    )
  }

  readr::write_csv(vif_joint_df, file.path(PATH_OUT, "TableE3_EIU_5dims_joint_VIF.csv"))

  kappa_joint <- kappa(model.matrix(m_eiu_5), exact = TRUE)
  write_lines(
    c(
      "Simultaneous EIU main test",
      paste0("N = ", nobs(m_eiu_5)),
      paste0("Condition number (kappa) = ", round(kappa_joint, 2))
    ),
    "TableE3_EIU_5dims_joint_kappa.txt"
  )

} else {
  message("Simultaneous EIU test skipped: not all 5 EIU dimensions available.")
}

# Table 3b (Region FE + clustered SE by region) ----
if ("region" %in% names(dat_cc)) {
  dat_cc$region <- as.factor(dat_cc$region)

  m_overall_regFE <- lm(mk_formula("life_satisfaction", c("eiu_overall", ctrl_vars, "region")), data = dat_cc)
  dim_regFE <- list()
  for (v in eiu_dims) dim_regFE[[v]] <- lm(mk_formula("life_satisfaction", c(v, ctrl_vars, "region")), data = dat_cc)

  models_regFE <- c(list("EIU overall + ctrls + region FE"=m_overall_regFE),
                    setNames(dim_regFE, paste0("EIU dim + ctrls + region FE: ", eiu_dims)))

  ft_t3b <- modelsummary::modelsummary(
    models_regFE,
    vcov = function(m) vc_cluster_from_data(m, dat_cc, "region"),
    stars = TRUE,
    statistic = "({std.error})",
    output = "flextable"
  )

  doc_t3b <- officer::read_docx() |>
    officer::body_add_par("Table 3b. Region fixed effects + region-clustered SE (HC3)", style="heading 1") |>
    flextable::body_add_flextable(ft_t3b)
  print(doc_t3b, target = file.path(PATH_OUT, "Table3b_EIU_dims_OLS_regionFE_cluster.docx"))

  
  b <- coef(m_overall_regFE)["eiu_overall"]
  se <- sqrt(diag(vc_cluster_from_data(m_overall_regFE, dat_cc, "region")))[ "eiu_overall" ]
  write_lines(c(
    "Region-FE, region-clustered (HC3):",
    paste0("beta(eiu_overall) = ", round(b, 4)),
    paste0("se_cluster = ", round(se, 4))
  ), "02c_m_ctrl_region_clustered.txt")
}

# Figure 1 (Standardized betas for EIU dimensions) ----
label_map <- c(
  electoral_process_pluralism = "Electoral process & pluralism",
  functioning_government      = "Functioning of government",
  political_participation     = "Political participation",
  political_culture           = "Political culture",
  civil_liberties             = "Civil liberties"
)

std_one <- function(x) as.numeric(scale(x))

get_std_beta <- function(term){
  cols <- unique(c("life_satisfaction", term, ctrl_vars))
  d <- dat_cc |> dplyr::select(dplyr::all_of(cols)) |> tidyr::drop_na()
  d_std <- d |> dplyr::mutate(dplyr::across(where(is.numeric), std_one))
  m <- lm(mk_formula("life_satisfaction", c(term, ctrl_vars)), data = d_std)
  V <- vc_hc3(m)
  b  <- coef(m)[term]
  se <- sqrt(diag(V))[term]
  tibble::tibble(term=term, beta=as.numeric(b), se=as.numeric(se),
                 lo=as.numeric(b-1.96*se), hi=as.numeric(b+1.96*se), n=nobs(m))
}

if (length(eiu_dims) > 0) {
  std_df <- purrr::map_dfr(eiu_dims, get_std_beta) |>
    dplyr::mutate(label = dplyr::recode(term, !!!label_map))
  readr::write_csv(std_df, file.path(PATH_OUT, "Figure1_EIU_std_coefs.csv"))

  p1 <- ggplot2::ggplot(std_df, ggplot2::aes(x = beta, y = reorder(label, beta))) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = hi), width = 0.2) +
    ggplot2::labs(
      title = "Figure 1. Standardized effects of EIU dimensions",
      subtitle = "Each coefficient from a separate OLS model with the same controls (HC3 SE)",
      x = "Standardized coefficient (β)", y = ""
    ) +
    ggplot2::theme_minimal()

  save_plot(p1, "figure_01_EIU_dims_std_beta.png")
}

# Figure 5 Spline (EIU overall) + spline lnGDP (optional) ----
spline_terms <- c("splines::ns(eiu_overall, df=3)", ctrl_vars)
m_spline <- lm(mk_formula("life_satisfaction", spline_terms), data = dat_cc)

grid <- data.frame(eiu_overall = seq(min(dat_cc$eiu_overall, na.rm=TRUE),
                                     max(dat_cc$eiu_overall, na.rm=TRUE),
                                     length.out=200))
for (cv in ctrl_vars) grid[[cv]] <- mean(dat_cc[[cv]], na.rm=TRUE)
grid$pred <- predict(m_spline, newdata = grid)

p_spline <- ggplot2::ggplot(dat_cc, ggplot2::aes(x = eiu_overall, y = life_satisfaction)) +
  ggplot2::geom_point(alpha=0.6) +
  ggplot2::geom_line(data=grid, ggplot2::aes(x=eiu_overall, y=pred), linewidth=1.1) +
  ggplot2::labs(title="Figure 5. Spline fit: EIU overall and life satisfaction",
                subtitle="Controls held at sample means",
                x="EIU overall", y="Life satisfaction (0–10)") +
  ggplot2::theme_minimal()
save_plot(p_spline, "fig05_spline_eiu.png")

if ("ln_gdp_pc_ppp" %in% names(dat_cc) && any(is.finite(dat_cc$ln_gdp_pc_ppp))) {
  m_spl_gdp <- lm(mk_formula("life_satisfaction",
                            c("splines::ns(ln_gdp_pc_ppp, df=3)",
                              setdiff(ctrl_vars, "ln_gdp_pc_ppp"),
                              "eiu_overall")), data = dat_cc)
  g2 <- data.frame(ln_gdp_pc_ppp = seq(min(dat_cc$ln_gdp_pc_ppp, na.rm=TRUE),
                                       max(dat_cc$ln_gdp_pc_ppp, na.rm=TRUE),
                                       length.out=200))
  g2$eiu_overall <- mean(dat_cc$eiu_overall, na.rm=TRUE)
  for (cv in setdiff(ctrl_vars, "ln_gdp_pc_ppp")) g2[[cv]] <- mean(dat_cc[[cv]], na.rm=TRUE)
  g2$pred <- predict(m_spl_gdp, newdata = g2)

  p_spl2 <- ggplot2::ggplot(dat_cc, ggplot2::aes(x = ln_gdp_pc_ppp, y = life_satisfaction)) +
    ggplot2::geom_point(alpha=0.6) +
    ggplot2::geom_line(data=g2, ggplot2::aes(x=ln_gdp_pc_ppp, y=pred), linewidth=1.1) +
    ggplot2::labs(title="Spline fit: ln(GDP pc PPP) and life satisfaction",
                  subtitle="EIU overall + other controls at means",
                  x="ln(GDP pc PPP)", y="Life satisfaction (0–10)") +
    ggplot2::theme_minimal()
  save_plot(p_spl2, "spline_lngdp.png")
}

# Multicollinearity diagnostics: VIF, condition number, sensitivity ----
rhs_all <- c(eiu_dims, ctrl_vars)
if (length(rhs_all) >= 2) {
  m_vif <- lm(mk_formula("life_satisfaction", rhs_all), data = dat_cc)
  vif_vals <- car::vif(m_vif)

  # car::vif kann Matrix (GVIF) liefern -> robust umwandeln
  if (is.matrix(vif_vals)) {
    # GVIF^(1/(2*Df)) as standardized measure
    vif_df <- tibble::tibble(
      var = rownames(vif_vals),
      GVIF = vif_vals[,1],
      Df   = vif_vals[,2],
      VIF  = vif_vals[,3]
    )
  } else {
    vif_df <- tibble::tibble(var = names(vif_vals), VIF = as.numeric(vif_vals))
  }

  readr::write_csv(vif_df, file.path(PATH_OUT, "_VIF_allDims.csv"))
  flag <- ifelse(any(vif_df$VIF >= 5, na.rm=TRUE), "FLAG: some VIF >= 5.", "OK: all VIF < 5.")
  write_lines(flag, "_VIF_flag.txt")

  # Condition number (Kappa) der Designmatrix
  X <- model.matrix(m_vif)
  kappa_val <- kappa(X, exact = TRUE)
  write_lines(paste0("Condition number (kappa) = ", round(kappa_val, 2)),
              "_condition_number_kappa.txt")

  # Sensitivity: drop high-VIF controls (>=5) and re-estimate EIU overall
  if ("eiu_overall" %in% rhs_all) {
    high_vif_ctrls <- vif_df |>
      dplyr::filter(.data$var %in% ctrl_vars, .data$VIF >= 5) |>
      dplyr::pull(.data$var)
    ctrl_sens <- setdiff(ctrl_vars, high_vif_ctrls)
    m_sens <- lm(mk_formula("life_satisfaction", c("eiu_overall", ctrl_sens)), data = dat_cc)

    V0 <- vc_hc3(m_overall_ctl)
    V1 <- vc_hc3(m_sens)
    out_sens <- tibble::tibble(
      model = c("baseline_controls","drop_highVIF_controls"),
      beta  = c(coef(m_overall_ctl)["eiu_overall"], coef(m_sens)["eiu_overall"]),
      se_hc3 = c(sqrt(diag(V0))["eiu_overall"], sqrt(diag(V1))["eiu_overall"]),
      n      = c(nobs(m_overall_ctl), nobs(m_sens)),
      dropped_controls = c(NA_character_, paste(high_vif_ctrls, collapse=", "))
    )
    readr::write_csv(out_sens, file.path(PATH_OUT, "_multicol_sensitivity_dropHighVIFControls.csv"))
  }
}

# VIF for V-Dem indices (optional)
if (USE_VDEM && length(vdem_vars) >= 2) {
  dat_vif_vdem <- dat_mod |>
    dplyr::select(dplyr::all_of(c("life_satisfaction", vdem_vars, ctrl_vars))) |>
    tidyr::drop_na()
  if (nrow(dat_vif_vdem) >= 10) {
    m_vif_vdem <- lm(mk_formula("life_satisfaction", c(vdem_vars, ctrl_vars)), data = dat_vif_vdem)
    vv <- car::vif(m_vif_vdem)
    if (is.matrix(vv)) {
      vv_df <- tibble::tibble(var = rownames(vv), GVIF = vv[,1], Df = vv[,2], VIF = vv[,3])
    } else {
      vv_df <- tibble::tibble(var = names(vv), VIF = as.numeric(vv))
    }
    readr::write_csv(vv_df, file.path(PATH_OUT, "_VIF_vdem.csv"))
  }
}

# Correlations (Figure A1): corrplot + heatmap + long csv ----
corr_vars <- unique(c("life_satisfaction","eiu_overall", eiu_dims, ctrl_vars))
corr_vars <- corr_vars[corr_vars %in% names(dat_cc)]
Xcorr <- dat_cc |> dplyr::select(dplyr::all_of(corr_vars)) |> as.data.frame()

cP <- suppressWarnings(cor(Xcorr, use="pairwise.complete.obs", method="pearson"))
cS <- suppressWarnings(cor(Xcorr, use="pairwise.complete.obs", method="spearman"))

cLong <- reshape2::melt(cP, varnames=c("var1","var2"), value.name="corr_pearson")
cLong2 <- reshape2::melt(cS, varnames=c("var1","var2"), value.name="corr_spearman")
cLong <- dplyr::left_join(cLong, cLong2, by=c("var1","var2"))
readr::write_csv(cLong, file.path(PATH_OUT, "FigureA1_Correlation_Matrix_long.csv"))

readr::write_csv(cLong, file.path(PATH_OUT, "01_correlations.csv"))

png(file.path(FIG_DIR, "corr_pearson.png"), width=1400, height=1100, res=180)
corrplot::corrplot(cP, method="color", type="upper", tl.cex=0.8, addCoef.col="black", number.cex=0.7)
dev.off()

png(file.path(FIG_DIR, "corr_spearman.png"), width=1400, height=1100, res=180)
corrplot::corrplot(cS, method="color", type="upper", tl.cex=0.8, addCoef.col="black", number.cex=0.7)
dev.off()

hm <- reshape2::melt(cP)
p_hm <- ggplot2::ggplot(hm, ggplot2::aes(x=Var1, y=Var2, fill=value)) +
  ggplot2::geom_tile() +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle=45, hjust=1)) +
  ggplot2::labs(title="Figure A1. Correlation heatmap (Pearson)", x="", y="", fill="corr")

save_plot(p_hm, "figure_A1_corr_heatmap.png", w=10, h=8)
save_plot(p_hm, "figure_A1_correlation_heatmap.png", w=10, h=8)

png(file.path(FIG_DIR, "figure_A1_corrplot.png"), width=1400, height=1100, res=180)
corrplot::corrplot(cP, method="color", type="upper", tl.cex=0.8)
dev.off()

# Alias im alten Namen:
if (file.exists(file.path(FIG_DIR, "figure_A1_correlation_heatmap.png"))) {
  file.copy(
    file.path(FIG_DIR, "figure_A1_correlation_heatmap.png"),
    file.path(FIG_DIR, "fig01_correlation_heatmap.png"),
    overwrite = TRUE
  )
}

# World map (Figure 02): EIU participation ----
if (HAS_SF && "political_participation" %in% names(dat_mod)) {
  world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
  map_df <- dat_mod |> dplyr::select(iso3c, political_participation) |> dplyr::distinct()
  world2 <- world |> dplyr::left_join(map_df, by=c("iso_a3"="iso3c"))

  p_map <- ggplot2::ggplot(world2) +
    ggplot2::geom_sf(ggplot2::aes(fill = political_participation), color="grey60", linewidth=0.1) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Figure 2. EIU political participation (country coverage)",
      subtitle = "EIU dimension values merged by ISO3",
      fill = "EIU participation"
    )
  save_plot(p_map, "figure_02_world_map_EIU_participation.png", w=11, h=6)
} else {
  message("World map skipped (sf/rnaturalearth not available or participation column missing).")
}

# Income subsamples + interaction (Figure A2 + TableR + fig08) ----
make_tertiles <- function(x){
  qs <- quantile(x, probs=c(1/3, 2/3), na.rm=TRUE)
  cut(x, breaks=c(-Inf, qs[1], qs[2], Inf), labels=c("Low","Middle","High"), include.lowest=TRUE)
}

if ("ln_gdp_pc_ppp" %in% names(dat_cc) && any(is.finite(dat_cc$ln_gdp_pc_ppp))) {
  dat_cc$income_tertile <- make_tertiles(dat_cc$ln_gdp_pc_ppp)
} else if ("income_group" %in% names(dat_cc)) {
  dat_cc$income_tertile <- as.factor(dat_cc$income_group)
} else {
  dat_cc$income_tertile <- NA
}

if (!all(is.na(dat_cc$income_tertile))) {

  get_std_beta_group <- function(g){
    d <- dat_cc |> dplyr::filter(income_tertile == g) |>
      dplyr::select(dplyr::all_of(unique(c("life_satisfaction","eiu_overall", ctrl_vars)))) |>
      tidyr::drop_na()
    d_std <- d |> dplyr::mutate(dplyr::across(where(is.numeric), std_one))
    m <- lm(mk_formula("life_satisfaction", c("eiu_overall", ctrl_vars)), data = d_std)
    V <- vc_hc3(m)
    b  <- coef(m)["eiu_overall"]; se <- sqrt(diag(V))["eiu_overall"]
    tibble::tibble(group=as.character(g), beta=as.numeric(b), se=as.numeric(se),
                   lo=as.numeric(b-1.96*se), hi=as.numeric(b+1.96*se), n=nobs(m))
  }

  grps <- levels(droplevels(as.factor(dat_cc$income_tertile)))
  inc_df <- purrr::map_dfr(grps, get_std_beta_group)
  readr::write_csv(inc_df, file.path(PATH_OUT, "FigureA2_IncomeGroups_std_coefs.csv"))

  p_inc <- ggplot2::ggplot(inc_df, ggplot2::aes(x=beta, y=reorder(group, beta))) +
    ggplot2::geom_vline(xintercept=0, linetype=2) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin=lo, xmax=hi), width=0.2) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title="Figure A2. Standardized EIU overall effect by income tertile",
                  subtitle="Separate OLS per tertile, same controls (HC3)",
                  x="Standardized coefficient (β)", y="")
  save_plot(p_inc, "figure_A2_income_tertiles_coeffs.png", w=9.5, h=5.8)
  save_plot(p_inc, "figure_A2_income_tertile_coefs.png", w=9.5, h=5.8)
  save_plot(p_inc, "fig_income_tertiles.png", w=9.5, h=5.8)

  # Interaction model table
  m_int <- lm(mk_formula("life_satisfaction", c("eiu_overall*income_tertile", ctrl_vars)), data = dat_cc)
  ft_int <- modelsummary::modelsummary(
    list("EIU overall × income tertile" = m_int),
    vcov = function(m) vc_hc3(m),
    stars = TRUE,
    statistic = "({std.error})",
    output = "flextable"
  )
  doc_int <- officer::read_docx() |>
    officer::body_add_par("Table R. Interaction: EIU overall by income group", style="heading 1") |>
    flextable::body_add_flextable(ft_int)
  print(doc_int, target = file.path(PATH_OUT, "TableR_Interaction_EIU_by_IncomeGroup.docx"))

  # FIXED interaction plot grid (kein Zeilen-Mismatch)
  xseq <- seq(min(dat_cc$eiu_overall, na.rm=TRUE),
              max(dat_cc$eiu_overall, na.rm=TRUE),
              length.out=200)

  gridI <- expand.grid(
    eiu_overall     = xseq,
    income_tertile  = grps,
    KEEP.OUT.ATTRS  = FALSE,
    stringsAsFactors = FALSE
  )
  for (cv in ctrl_vars) gridI[[cv]] <- mean(dat_cc[[cv]], na.rm = TRUE)
  gridI$pred <- predict(m_int, newdata = gridI)

  p_int <- ggplot2::ggplot(gridI, ggplot2::aes(x=eiu_overall, y=pred, linetype=income_tertile)) +
    ggplot2::geom_line(linewidth=1.0) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title="Interaction: EIU overall × income tertile",
                  x="EIU overall", y="Predicted life satisfaction",
                  linetype="Income tertile")
  save_plot(p_int, "fig08_interaction_income.png", w=9.5, h=5.8)
}

# Diagnostics: residual plots, Cook's D, leverage, DFBETAS ----
m_base <- m_overall_ctl

# Residuals vs fitted (+ loess)
p_res <- ggplot2::ggplot(data.frame(fitted=fitted(m_base), resid=residuals(m_base)),
                         ggplot2::aes(x=fitted, y=resid)) +
  ggplot2::geom_point(alpha=0.7) +
  ggplot2::geom_hline(yintercept=0, linetype=2) +
  ggplot2::geom_smooth(method="loess", se=FALSE) +
  ggplot2::theme_minimal() +
  ggplot2::labs(title="Residuals vs fitted (baseline)", x="Fitted", y="Residual")
save_plot(p_res, "fig03_residuals_fitted.png", w=8.5, h=5.5)

# QQ
qq <- data.frame(sample=residuals(m_base))
p_qq <- ggplot2::ggplot(qq, ggplot2::aes(sample=sample)) +
  ggplot2::stat_qq() +
  ggplot2::stat_qq_line() +
  ggplot2::theme_minimal() +
  ggplot2::labs(title="QQ plot of residuals (baseline)", x="Theoretical quantiles", y="Sample quantiles")
save_plot(p_qq, "fig04_residuals_qq.png", w=8.5, h=5.5)

# Neue Diagnostikplots
hat   <- hatvalues(m_base)
cook  <- cooks.distance(m_base)
rstd  <- rstandard(m_base)
rstud <- rstudent(m_base)
p_par <- length(coef(m_base))
n_obs <- nobs(m_base)

# D1: Scale-Location
df_sl <- data.frame(fitted = fitted(m_base), sqrt_abs_rstd = sqrt(abs(rstd)))
p_sl <- ggplot2::ggplot(df_sl, ggplot2::aes(x=fitted, y=sqrt_abs_rstd)) +
  ggplot2::geom_point(alpha=0.7) +
  ggplot2::geom_smooth(method="loess", se=FALSE) +
  ggplot2::theme_minimal() +
  ggplot2::labs(title="Scale-Location", x="Fitted", y="sqrt(|standardized residuals|)")
save_plot(p_sl, "figD1_scale_location.png", w=8.5, h=5.5)

# D2: Cook's distance by index
df_cd <- data.frame(idx = seq_along(cook), cook = cook)
thr_cd <- 4/n_obs
p_cd <- ggplot2::ggplot(df_cd, ggplot2::aes(x=idx, y=cook)) +
  ggplot2::geom_col() +
  ggplot2::geom_hline(yintercept = thr_cd, linetype=2) +
  ggplot2::theme_minimal() +
  ggplot2::labs(title="Cook's distance (baseline)", subtitle=paste0("Threshold 4/n = ", signif(thr_cd, 3)),
                x="Observation index", y="Cook's D")
save_plot(p_cd, "figD2_cooks_distance.png", w=9.0, h=5.5)

# D3: Leverage (hat) by index
df_hat <- data.frame(idx = seq_along(hat), hat = hat)
thr_hat <- 2*p_par/n_obs
p_hat <- ggplot2::ggplot(df_hat, ggplot2::aes(x=idx, y=hat)) +
  ggplot2::geom_col() +
  ggplot2::geom_hline(yintercept = thr_hat, linetype=2) +
  ggplot2::theme_minimal() +
  ggplot2::labs(title="Leverage (hat values) (baseline)",
                subtitle=paste0("Threshold 2p/n = ", signif(thr_hat, 3)),
                x="Observation index", y="Leverage (hat)")
save_plot(p_hat, "figD3_leverage_hat.png", w=9.0, h=5.5)

# D4: Residuals vs Leverage with Cook contours (approx)
df_rl <- data.frame(hat = hat, rstud = rstud)
hseq <- seq(min(hat, na.rm=TRUE) + 1e-6, max(hat, na.rm=TRUE) - 1e-6, length.out=200)
cook_levels <- c(0.5, 1)
contours <- purrr::map_dfr(cook_levels, function(D){
  r <- sqrt(D * p_par * (1 - hseq) / hseq)
  tibble::tibble(hat=hseq, r= r, D=D, sign="+") |>
    dplyr::bind_rows(tibble::tibble(hat=hseq, r= -r, D=D, sign="-"))
})

p_rl <- ggplot2::ggplot(df_rl, ggplot2::aes(x=hat, y=rstud)) +
  ggplot2::geom_point(alpha=0.7) +
  ggplot2::geom_hline(yintercept=0, linetype=2) +
  ggplot2::geom_vline(xintercept=thr_hat, linetype=2) +
  ggplot2::geom_line(data=contours, ggplot2::aes(x=hat, y=r, group=interaction(D,sign)), linewidth=0.6) +
  ggplot2::theme_minimal() +
  ggplot2::labs(title="Residuals vs Leverage (baseline)",
                subtitle="Cook contours (approx) for D = 0.5 and 1; vertical line = 2p/n",
                x="Leverage (hat)", y="Studentized residuals")
save_plot(p_rl, "figD4_resid_vs_leverage_cook.png", w=9.0, h=5.8)

# DFBETAS for EIU overall
if ("eiu_overall" %in% names(coef(m_base))) {
  dfb <- dfbetas(m_base)[, "eiu_overall"]
  df_dfb <- data.frame(idx = seq_along(dfb), dfbetas = as.numeric(dfb))
  thr_dfb <- 2/sqrt(n_obs)
  p_dfb <- ggplot2::ggplot(df_dfb, ggplot2::aes(x=idx, y=dfbetas)) +
    ggplot2::geom_hline(yintercept=c(-thr_dfb, thr_dfb), linetype=2) +
    ggplot2::geom_point(alpha=0.7) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title="DFBETAS for EIU overall (baseline)",
                  subtitle=paste0("Rule-of-thumb ±2/sqrt(n) = ±", signif(thr_dfb, 3)),
                  x="Observation index", y="DFBETAS (EIU overall)")
  save_plot(p_dfb, "figD5_dfbetas_eiu_overall.png", w=9.0, h=5.5)
}

# D6: Influence plot (car)
png(file.path(FIG_DIR, "figD6_influencePlot.png"), width=1400, height=1100, res=180)
suppressWarnings(try(car::influencePlot(m_base,
                                       main="Influence Plot (baseline)",
                                       sub="Circle size ~ Cook's distance"),
                     silent = TRUE))
dev.off()

# Robustness Table 4 ----
# ---- Robust coefficient extraction (handles backticks, wrappers, renamed terms) ----
escape_regex <- function(x){
  x <- gsub("\\", "\\\\", x)
  gsub("([\.^$|()\[\]{}*+?])", "\\\1", x, perl = TRUE)
}

resolve_term_in_names <- function(term, nms){
  if (is.null(nms) || length(nms) == 0) return(NA_character_)
  if (term %in% nms) return(term)
  bt <- paste0("`", term, "`")
  if (bt %in% nms) return(bt)
  for (w in c(paste0("scale(", term, ")"), paste0("I(", term, ")"))) {
    if (w %in% nms) return(w)
  }
  # token regex (avoid partial matches)
  tesc <- escape_regex(term)
  pat  <- paste0("(^|[^A-Za-z0-9_])`?", tesc, "`?($|[^A-Za-z0-9_])")
  hits <- nms[grepl(pat, nms)]
  if (length(hits) >= 1) {
    strip <- gsub("`", "", hits, fixed = TRUE)
    if (any(strip == term)) return(hits[which(strip == term)[1]])
    return(hits[1])
  }
  NA_character_
}

resolve_term_in_model <- function(m, term){
  resolve_term_in_names(term, names(stats::coef(m)))
}

get_b_se_p_hc3 <- function(m, term){
  cn <- resolve_term_in_model(m, term)
  if (is.na(cn)) return(c(b = NA_real_, se = NA_real_, p = NA_real_))
  V  <- vc_hc3(m)
  b  <- unname(stats::coef(m)[cn])
  se <- unname(sqrt(diag(V))[cn])
  if (is.na(b) || is.na(se) || se == 0) return(c(b = b, se = se, p = NA_real_))
  t  <- b / se
  p  <- 2 * stats::pt(abs(t), df = stats::df.residual(m), lower.tail = FALSE)
  c(b = b, se = se, p = p)
}

get_b_se_p_rreg <- function(mr, term){
  if (is.null(mr)) return(c(b = NA_real_, se = NA_real_, p = NA_real_))
  smr <- tryCatch(summary(mr), error = function(e) NULL)
  if (is.null(smr) || is.null(smr$coefficients)) return(c(b = NA_real_, se = NA_real_, p = NA_real_))
  rn <- rownames(smr$coefficients)
  cn <- resolve_term_in_names(term, rn)
  if (is.na(cn)) return(c(b = NA_real_, se = NA_real_, p = NA_real_))
  b  <- unname(stats::coef(mr)[cn])
  se <- unname(smr$coefficients[cn, "Std. Error"])
  if (is.na(b) || is.na(se) || se == 0) return(c(b = b, se = se, p = NA_real_))
  p  <- 2 * stats::pnorm(abs(b / se), lower.tail = FALSE)
  c(b = b, se = se, p = p)
}

drop_influential <- function(m){
  cd <- cooks.distance(m)
  thr <- 4/length(cd)
  keep <- which(cd <= thr | is.na(cd))
  list(keep=keep, thr=thr, cd=cd)
}

jackknife_term <- function(term, f_terms, data_use){
  n <- nrow(data_use)
  betas <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    d_i <- data_use[-i, , drop=FALSE]
    m_i <- lm(mk_formula("life_satisfaction", f_terms), data = d_i)
    cn_i <- resolve_term_in_model(m_i, term)
    betas[i] <- if (!is.na(cn_i) && (cn_i %in% names(coef(m_i)))) unname(coef(m_i)[cn_i]) else NA_real_
  }
  tibble::tibble(
    term = term,
    jack_mean = mean(betas, na.rm=TRUE),
    jack_sd   = sd(betas, na.rm=TRUE),
    jack_min  = min(betas, na.rm=TRUE),
    jack_max  = max(betas, na.rm=TRUE),
    jack_sign_share = mean(sign(betas) == sign(mean(betas, na.rm=TRUE)), na.rm=TRUE)
  )
}

exclude_top10 <- function(d){
  ord <- order(d$life_satisfaction, decreasing=TRUE)
  drop <- ord[seq_len(min(10, length(ord)))]
  d[-drop, , drop=FALSE]
}

roll_term <- function(term){
  f_terms <- c(term, ctrl_vars)

  # baseline
  m0 <- lm(mk_formula("life_satisfaction", f_terms), data = dat_cc)
  b0 <- get_b_se_p_hc3(m0, term)
  # rlm
  mr <- tryCatch(MASS::rlm(mk_formula("life_satisfaction", f_terms), data = dat_cc, maxit=200),
                 error=function(e) NULL)
  rr <- get_b_se_p_rreg(mr, term)

  # excluding influential (Cook's D based on m0)
  infl <- drop_influential(m0)
  d_ex <- dat_cc[infl$keep, , drop=FALSE]
  m_ex <- lm(mk_formula("life_satisfaction", f_terms), data = d_ex)
  b_ex <- get_b_se_p_hc3(m_ex, term)

  # excluding top10 DV
  d_t10 <- exclude_top10(dat_cc)
  m_t10 <- lm(mk_formula("life_satisfaction", f_terms), data = d_t10)
  b_t10 <- get_b_se_p_hc3(m_t10, term)

  # jackknife summary
  j <- jackknife_term(term, f_terms, dat_cc)

  tibble::tibble(
    term = term,
    beta_ols = b0["b"], se_ols = b0["se"], p_ols = b0["p"],
    beta_rreg = rr["b"], se_rreg = rr["se"], p_rreg = rr["p"],
    beta_excl_infl = b_ex["b"], se_excl_infl = b_ex["se"], p_excl_infl = b_ex["p"],
    beta_excl_top10 = b_t10["b"], se_excl_top10 = b_t10["se"], p_excl_top10 = b_t10["p"],
    jack_sign_share = j$jack_sign_share,
    N = nobs(m0),
    N_excl_infl = nobs(m_ex),
    N_excl_top10 = nobs(m_t10)
  )
}

terms_for_roll <- c("eiu_overall", eiu_dims)
roll_df <- purrr::map_dfr(terms_for_roll, roll_term)

# FDR correction across dimensions
if (length(eiu_dims) > 0) {
  pdims <- roll_df$p_ols[roll_df$term %in% eiu_dims]
  roll_df$p_fdr_dims <- NA_real_
  roll_df$p_fdr_dims[roll_df$term %in% eiu_dims] <- p.adjust(pdims, method="fdr")
}

readr::write_csv(roll_df, file.path(PATH_OUT, "Table4_robustness_rollup.csv"))

ft_roll <- flextable::flextable(roll_df) |>
  flextable::autofit() |>
  flextable::align(align="center", part="all") |>
  flextable::bold(part="header")

doc_roll <- officer::read_docx() |>
  officer::body_add_par("Table 4. Robustness checks (EIU overall + dimensions)", style="heading 1") |>
  flextable::body_add_flextable(ft_roll)
print(doc_roll, target = file.path(PATH_OUT, "Table4_Robustness.docx"))

# Table A2 (OLS vs rreg) for dims ----
ols_list <- setNames(dim_models, paste0("OLS: ", eiu_dims))
rreg_list <- list()
for (v in eiu_dims) {
  rreg_list[[paste0("RREG: ", v)]] <- tryCatch(
    MASS::rlm(mk_formula("life_satisfaction", c(v, ctrl_vars)), data = dat_cc, maxit=200),
    error=function(e) NULL
  )
}
rreg_list <- rreg_list[!purrr::map_lgl(rreg_list, is.null)]

ft_a2_ols <- modelsummary::modelsummary(
  ols_list,
  vcov = function(m) vc_hc3(m),
  stars = TRUE,
  statistic = "({std.error})",
  output = "flextable"
)

doc_a2 <- officer::read_docx() |>
  officer::body_add_par("Table A2. OLS (HC3) – EIU dimensions", style="heading 1") |>
  flextable::body_add_flextable(ft_a2_ols)

if (length(rreg_list) > 0) {
  ft_a2_r <- modelsummary::modelsummary(
    rreg_list,
    stars = TRUE,
    statistic = "({std.error})",
    output = "flextable"
  )
  doc_a2 <- doc_a2 |>
    officer::body_add_par("Table A2 (cont.). Robust regression (rlm) – EIU dimensions", style="heading 1") |>
    flextable::body_add_flextable(ft_a2_r)
}
print(doc_a2, target = file.path(PATH_OUT, "TableA2_EIU_OLS_rreg.docx"))

# Table A4 Jackknife (summary) + leaves csvs ----
jack_all <- list()

for (term in terms_for_roll) {
  f_terms <- c(term, ctrl_vars)
  n <- nrow(dat_cc)
  betas <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    d_i <- dat_cc[-i, , drop=FALSE]
    m_i <- lm(mk_formula("life_satisfaction", f_terms), data = d_i)
    betas[i] <- coef(m_i)[term]
  }
  leaves <- tibble::tibble(
    left_out_country = dat_cc$country,
    left_out_iso3 = dat_cc$iso3c,
    beta_loo = betas
  )
  fn <- paste0("02e_jackknife_leaves_", term, ".csv")
  readr::write_csv(leaves, file.path(PATH_OUT, fn))

  jack_all[[term]] <- tibble::tibble(
    term=term,
    mean=mean(betas, na.rm=TRUE),
    sd=sd(betas, na.rm=TRUE),
    min=min(betas, na.rm=TRUE),
    max=max(betas, na.rm=TRUE),
    sign_share = mean(sign(betas) == sign(mean(betas, na.rm=TRUE)), na.rm=TRUE)
  )
}

jack_sum <- dplyr::bind_rows(jack_all)
ft_j <- flextable::flextable(jack_sum) |>
  flextable::autofit() |>
  flextable::align(align="center", part="all") |>
  flextable::bold(part="header")
doc_j <- officer::read_docx() |>
  officer::body_add_par("Table A4. Jackknife (leave-one-out) summary", style="heading 1") |>
  flextable::body_add_flextable(ft_j)
print(doc_j, target = file.path(PATH_OUT, "TableA4_EIU_Jackknife.docx"))

file.copy(file.path(PATH_OUT, "02e_jackknife_leaves_eiu_overall.csv"),
          file.path(PATH_OUT, "02e_jackknife_beta_eiu_ctrl_leaves.csv"),
          overwrite = TRUE)

# Jackknife Plot
leaves_path <- file.path(PATH_OUT, "02e_jackknife_leaves_eiu_overall.csv")
if (file.exists(leaves_path)) {
  leaves <- readr::read_csv(leaves_path, show_col_types = FALSE)
  beta_full <- unname(coef(m_overall_ctl)["eiu_overall"])
  ci <- stats::quantile(leaves$beta_loo, probs = c(0.025, 0.975), na.rm = TRUE)

  p_jack <- ggplot2::ggplot(leaves, ggplot2::aes(x = seq_along(beta_loo), y = beta_loo)) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::geom_hline(yintercept = beta_full, linetype = 2) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Jackknife (leave-one-out): beta(EIU overall + controls)",
      subtitle = paste0("Full-sample β=", round(beta_full, 3),
                        " | 95% LOO interval [", round(ci[1], 3), ", ", round(ci[2], 3), "]"),
      x = "Left-out observation index",
      y = "β (EIU overall)"
    )

  save_plot(p_jack, "fig02e_jackknife_beta_eiu_ctrl.png", w=9.5, h=5.8)
}

# Table A5 Excluding Top10 DV ----
d_t10 <- exclude_top10(dat_cc)
m_overall_t10 <- lm(mk_formula("life_satisfaction", c("eiu_overall", ctrl_vars)), data = d_t10)
dim_t10 <- list()
for (v in eiu_dims) dim_t10[[v]] <- lm(mk_formula("life_satisfaction", c(v, ctrl_vars)), data = d_t10)

models_t10 <- c(list("EIU overall + ctrls (excl top10)"=m_overall_t10),
                setNames(dim_t10, paste0("EIU dim + ctrls (excl top10): ", eiu_dims)))

ft_t10 <- modelsummary::modelsummary(
  models_t10,
  vcov = function(m) vc_hc3(m),
  stars = TRUE,
  statistic = "({std.error})",
  output = "flextable"
)

doc_t10 <- officer::read_docx() |>
  officer::body_add_par("Table A5. Excluding top-10 life satisfaction observations", style="heading 1") |>
  flextable::body_add_flextable(ft_t10)
print(doc_t10, target = file.path(PATH_OUT, "TableA5_EIU_ExcludingTop10.docx"))

# Ridge + Lasso (Table R1 + CV figs) ----
xvars <- unique(c("eiu_overall", eiu_dims, ctrl_vars))
X <- dat_cc |> dplyr::select(dplyr::all_of(xvars)) |> as.data.frame()
mm <- model.matrix(~ . , data = X)
y <- dat_cc$life_satisfaction

cv_ridge <- glmnet::cv.glmnet(mm, y, alpha = 0, nfolds = 10, standardize = TRUE)
ridge_fit <- glmnet::glmnet(mm, y, alpha = 0, lambda = cv_ridge$lambda.min, standardize = TRUE)
ridge_coef <- as.matrix(coef(ridge_fit))
ridge_df <- tibble::tibble(term = rownames(ridge_coef), coef = as.numeric(ridge_coef[,1])) |>
  dplyr::filter(term != "(Intercept)")
readr::write_csv(ridge_df, file.path(PATH_OUT, "TableR1_Ridge_Coefficients.csv"))

ft_ridge <- flextable::flextable(ridge_df) |>
  flextable::autofit() |>
  flextable::align(align="center", part="all") |>
  flextable::bold(part="header")
doc_ridge_coef <- officer::read_docx() |>
  officer::body_add_par("Table R1. Ridge coefficients (lambda.min)", style="heading 1") |>
  flextable::body_add_flextable(ft_ridge)
print(doc_ridge_coef, target = file.path(PATH_OUT, "TableR1_Ridge_Coefficients.docx"))

ridge_info <- tibble::tibble(
  item = c("lambda.min","lambda.1se","n"),
  value = c(cv_ridge$lambda.min, cv_ridge$lambda.1se, length(y))
)
ft_ridge2 <- flextable::flextable(ridge_info) |> flextable::autofit()
doc_ridge <- officer::read_docx() |>
  officer::body_add_par("Table R1. Ridge (CV summary)", style="heading 1") |>
  flextable::body_add_flextable(ft_ridge2)
print(doc_ridge, target = file.path(PATH_OUT, "TableR1_Ridge.docx"))

png(file.path(FIG_DIR, "fig06_ridge_cv.png"), width=1400, height=900, res=180)
plot(cv_ridge)
dev.off()

cv_lasso <- glmnet::cv.glmnet(mm, y, alpha = 1, nfolds = 10, standardize = TRUE)
png(file.path(FIG_DIR, "fig07_lasso_cv.png"), width=1400, height=900, res=180)
plot(cv_lasso)
dev.off()

p_ridge <- ggplot2::ggplot(ridge_df, ggplot2::aes(x=coef, y=reorder(term, coef))) +
  ggplot2::geom_vline(xintercept=0, linetype=2) +
  ggplot2::geom_point() +
  ggplot2::theme_minimal() +
  ggplot2::labs(title="Figure R1. Ridge coefficients (lambda.min)", x="Coefficient", y="")
save_plot(p_ridge, "figure_R1_ridge_coefs.png", w=9.5, h=7)

# PCA (Tables R2/R3 + loadings plot) ----
if (length(eiu_dims) >= 3) {
  Z <- dat_cc |> dplyr::select(dplyr::all_of(eiu_dims)) |> as.data.frame()
  pca <- prcomp(Z, center=TRUE, scale.=TRUE)

  load <- tibble::tibble(var = rownames(pca$rotation), PC1 = pca$rotation[,1], PC2 = pca$rotation[,2])
  scores <- tibble::tibble(PC1 = pca$x[,1], PC2 = pca$x[,2])
  dat_cc$EIU_PC1 <- scores$PC1

  ft_load <- flextable::flextable(load) |> flextable::autofit() |> flextable::bold(part="header")
  doc_pcaL <- officer::read_docx() |>
    officer::body_add_par("Table R2. PCA loadings (EIU dimensions)", style="heading 1") |>
    flextable::body_add_flextable(ft_load)
  print(doc_pcaL, target = file.path(PATH_OUT, "TableR2_PCA_Loadings.docx"))

  p_load <- ggplot2::ggplot(load, ggplot2::aes(x=PC1, y=reorder(var, PC1))) +
    ggplot2::geom_vline(xintercept=0, linetype=2) +
    ggplot2::geom_point() +
    ggplot2::theme_minimal() +
    ggplot2::labs(title="Figure R2. PCA loadings (PC1)", x="Loading on PC1", y="")
  save_plot(p_load, "figure_R2_pca_loadings.png", w=9.5, h=5.8)
  save_plot(p_load, "pc1_loadings.png", w=9.5, h=5.8)

  pca_sum <- tibble::tibble(
    component = paste0("PC", 1:length(pca$sdev)),
    stdev = pca$sdev,
    var_explained = (pca$sdev^2) / sum(pca$sdev^2)
  )
  ft_pca <- flextable::flextable(head(pca_sum, 5)) |> flextable::autofit() |> flextable::bold(part="header")
  doc_pca <- officer::read_docx() |>
    officer::body_add_par("Table R2. PCA summary (EIU dimensions)", style="heading 1") |>
    flextable::body_add_flextable(ft_pca)
  print(doc_pca, target = file.path(PATH_OUT, "TableR2_PCA.docx"))

  m_pc1 <- lm(mk_formula("life_satisfaction", c("EIU_PC1", ctrl_vars)), data = dat_cc)
  ft_pc1m <- modelsummary::modelsummary(
    list("OLS: PC1 index + controls" = m_pc1),
    vcov = function(m) vc_hc3(m),
    stars = TRUE,
    statistic = "({std.error})",
    output = "flextable"
  )
  doc_pc1m <- officer::read_docx() |>
    officer::body_add_par("Table R3. OLS using PCA-based democracy index (PC1)", style="heading 1") |>
    flextable::body_add_flextable(ft_pc1m)
  print(doc_pc1m, target = file.path(PATH_OUT, "TableR3_PCA_Index_OLS.docx"))
}

# Optional Bootstrap beta (overall + controls) ----
if (DO_BOOTSTRAP) {
  f0 <- mk_formula("life_satisfaction", c("eiu_overall", ctrl_vars))
  bet <- rep(NA_real_, B_BOOT)
  for (b in seq_len(B_BOOT)) {
    idx <- sample(seq_len(nrow(dat_cc)), replace=TRUE)
    d_b <- dat_cc[idx, , drop=FALSE]
    m_b <- lm(f0, data = d_b)
    bet[b] <- coef(m_b)["eiu_overall"]
  }
  boot_df <- tibble::tibble(beta = bet)
  readr::write_csv(boot_df, file.path(PATH_OUT, "02d_bootstrap_beta_eiu_ctrl.csv"))

  ci <- quantile(boot_df$beta, probs=c(0.025,0.975), na.rm=TRUE)
  write_lines(paste0("Bootstrap 95% CI: [", round(ci[1],4), ", ", round(ci[2],4), "]"),
              "02d_bootstrap_beta_eiu_ctrl_ci.txt")

  p_boot <- ggplot2::ggplot(boot_df, ggplot2::aes(x=beta)) +
    ggplot2::geom_histogram(bins=40) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title="Bootstrap distribution of beta (EIU overall + controls)",
                  x="beta", y="count")
  save_plot(p_boot, "fig02d_bootstrap_beta_eiu_ctrl.png", w=8.5, h=5.5)
}

# Figure 2: coefficient plot EIU-only vs EIU+controls ----
tidy_hc3 <- function(m, model_name){
  V <- vc_hc3(m)
  est <- coef(m)
  se <- sqrt(diag(V))
  # align names robustly
  se <- se[names(est)]
  stat <- est / se
  pval <- 2*pt(abs(stat), df = df.residual(m), lower.tail = FALSE)
  tibble::tibble(
    term = names(est),
    estimate = as.numeric(est),
    std.error = as.numeric(se),
    statistic = as.numeric(stat),
    p.value = as.numeric(pval),
    model = model_name
  ) |>
    dplyr::mutate(
      lo = .data$estimate - 1.96*.data$std.error,
      hi = .data$estimate + 1.96*.data$std.error
    )
}

coef_df <- dplyr::bind_rows(
  tidy_hc3(m_overall_biv, "EIU only"),
  tidy_hc3(m_overall_ctl, "EIU + controls")
) |>
  dplyr::filter(.data$term != "(Intercept)")

readr::write_csv(coef_df, file.path(PATH_OUT, "02_coefplot_eiu_data.csv"))

p_coef <- ggplot2::ggplot(coef_df, ggplot2::aes(x = .data$estimate, y = reorder(.data$term, .data$estimate))) +
  ggplot2::geom_vline(xintercept = 0, linetype = 2) +
  ggplot2::geom_point() +
  ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data$lo, xmax = .data$hi), width = 0.2) +
  ggplot2::facet_wrap(~model, scales = "free_y") +
  ggplot2::theme_minimal() +
  ggplot2::labs(
    title = "Koeffizientenplot (EIU)",
    x = "Koeffizient (HC3 robuste SE)",
    y = ""
  )
save_plot(p_coef, "fig02_coef_eiu.png", w = 11, h = 6)

# V-Dem Robustheit (optional) ----
if (USE_VDEM && length(vdem_vars) >= 3) {

  # V-Dem sample: not forced to match EIU N, but documented
  need_vdem <- unique(c("country","iso3c","life_satisfaction", vdem_vars, ctrl_vars))
  dat_cc_vdem <- dat_mod |>
    dplyr::select(dplyr::all_of(need_vdem)) |>
    dplyr::filter(!is.na(.data$country) & !is.na(.data$iso3c)) |>
    tidyr::drop_na()

  readr::write_csv(tibble::tibble(N_vdem_sample = nrow(dat_cc_vdem)),
                   file.path(PATH_OUT, "00_vdem_sample_N.csv"))

  # Modelle: einzelne V-Dem-Dimensionen (jeweils separat, analog Fig.1)
  get_std_beta_vdem <- function(term){
    cols <- unique(c("life_satisfaction", term, ctrl_vars))
    d <- dat_cc_vdem |> dplyr::select(dplyr::all_of(cols)) |> tidyr::drop_na()
    d_std <- d |> dplyr::mutate(dplyr::across(where(is.numeric), std_one))
    m <- lm(mk_formula("life_satisfaction", c(term, ctrl_vars)), data = d_std)
    V <- vc_hc3(m)
    b  <- coef(m)[term]; se <- sqrt(diag(V))[term]
    tibble::tibble(term=term, beta=as.numeric(b), se=as.numeric(se),
                   lo=as.numeric(b-1.96*se), hi=as.numeric(b+1.96*se), n=nobs(m))
  }

  std_vdem <- purrr::map_dfr(vdem_vars, get_std_beta_vdem)
  readr::write_csv(std_vdem, file.path(PATH_OUT, "FigureV1_VDem_std_coefs.csv"))

  p_vdem_std <- ggplot2::ggplot(std_vdem, ggplot2::aes(x = beta, y = reorder(term, beta))) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2) +
    ggplot2::geom_point() +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = hi), width = 0.2) +
    ggplot2::labs(
      title = "Figure V1. Standardized effects of V-Dem indices",
      subtitle = "Each coefficient from a separate OLS model with the same controls (HC3 SE)",
      x = "Standardized coefficient (β)", y = ""
    ) +
    ggplot2::theme_minimal()
  save_plot(p_vdem_std, "figure_V1_vdem_std_betas.png", w=9.5, h=5.8)

  # 5-Dimensionsmodell (alle V-Dem-Indizes zusammen)
  m_vdem_5 <- lm(mk_formula("life_satisfaction", c(vdem_vars, ctrl_vars)), data = dat_cc_vdem)

  ft_vdem <- modelsummary::modelsummary(
    list("V-Dem indices + controls" = m_vdem_5),
    vcov = function(m) vc_hc3(m),
    stars = TRUE,
    statistic = "({std.error})",
    output = "flextable"
  )
  doc_vdem <- officer::read_docx() |>
    officer::body_add_par("Table V1. V-Dem indices + controls (HC3 robust SE)", style="heading 1") |>
    flextable::body_add_flextable(ft_vdem)
  print(doc_vdem, target = file.path(PATH_OUT, "TableV1_VDem_5indices_OLS.docx"))

  # Vergleich: EIU-Dims ohne K. vs V-Dem ohne K. vs V-Dem mit K.
  m_eiu_dims_only <- if (length(eiu_dims) >= 3) {
    lm(mk_formula("life_satisfaction", eiu_dims), data = dat_cc)
  } else NULL
  m_vdem_only <- lm(mk_formula("life_satisfaction", vdem_vars), data = dat_cc_vdem)

  models_cmp <- list()
  if (!is.null(m_eiu_dims_only)) models_cmp[["EIU dims (no ctrls)"]] <- m_eiu_dims_only
  models_cmp[["V-Dem (no ctrls)"]] <- m_vdem_only
  models_cmp[["V-Dem + ctrls"]] <- m_vdem_5

  ft_cmp <- modelsummary::modelsummary(
    models_cmp,
    vcov = function(m) vc_hc3(m),
    stars = TRUE,
    statistic = "({std.error})",
    output = "flextable"
  )
  doc_cmp <- officer::read_docx() |>
    officer::body_add_par("Table V2. Compare EIU vs V-Dem models", style="heading 1") |>
    flextable::body_add_flextable(ft_cmp)
  print(doc_cmp, target = file.path(PATH_OUT, "TableV2_Compare_EIU_VDem.docx"))

  # Coef plot (analog fig09 in altem Skript)
  Vtab <- tidy_hc3(m_vdem_5, "V-Dem indices + ctrls") |>
    dplyr::filter(stringr::str_detect(.data$term, "^vdem_"), .data$term != "(Intercept)")
  p_vdem <- ggplot2::ggplot(Vtab, ggplot2::aes(x=reorder(.data$term, .data$estimate), y=.data$estimate))+
    ggplot2::geom_point()+
    ggplot2::geom_errorbar(ggplot2::aes(ymin=.data$estimate-1.96*.data$std.error,
                                        ymax=.data$estimate+1.96*.data$std.error), width=.15)+
    ggplot2::coord_flip()+
    ggplot2::theme_minimal() +
    ggplot2::labs(x="", y="Coefficient (HC3)", title="V-Dem indices (with controls)")
  save_plot(p_vdem, "fig09_coef_vdem.png", w=8.5, h=5.5)

  # Robustness rollup (V-Dem): analog zu Table4
  roll_term_vdem <- function(term){
    f_terms <- c(term, ctrl_vars)

    m0 <- lm(mk_formula("life_satisfaction", f_terms), data = dat_cc_vdem)
    b0 <- get_b_se_p_hc3(m0, term)
    mr <- tryCatch(MASS::rlm(mk_formula("life_satisfaction", f_terms), data = dat_cc_vdem, maxit=200),
                   error=function(e) NULL)
    rr <- get_b_se_p_rreg(mr, term)

    infl <- drop_influential(m0)
    d_ex <- dat_cc_vdem[infl$keep, , drop=FALSE]
    m_ex <- lm(mk_formula("life_satisfaction", f_terms), data = d_ex)
    b_ex <- get_b_se_p_hc3(m_ex, term)

    d_t10 <- exclude_top10(dat_cc_vdem)
    m_t10 <- lm(mk_formula("life_satisfaction", f_terms), data = d_t10)
    b_t10 <- get_b_se_p_hc3(m_t10, term)

    j <- jackknife_term(term, f_terms, dat_cc_vdem)

    tibble::tibble(
      term = term,
      beta_ols = b0["b"], se_ols = b0["se"], p_ols = b0["p"],
      beta_rreg = rr["b"], se_rreg = rr["se"], p_rreg = rr["p"],
      beta_excl_infl = b_ex["b"], se_excl_infl = b_ex["se"], p_excl_infl = b_ex["p"],
      beta_excl_top10 = b_t10["b"], se_excl_top10 = b_t10["se"], p_excl_top10 = b_t10["p"],
      jack_sign_share = j$jack_sign_share,
      N = nobs(m0),
      N_excl_infl = nobs(m_ex),
      N_excl_top10 = nobs(m_t10)
    )
  }

  roll_vdem <- purrr::map_dfr(vdem_vars, roll_term_vdem)
  readr::write_csv(roll_vdem, file.path(PATH_OUT, "TableV3_VDem_robustness_rollup.csv"))

  ft_rv <- flextable::flextable(roll_vdem) |>
    flextable::autofit() |>
    flextable::align(align="center", part="all") |>
    flextable::bold(part="header")
  doc_rv <- officer::read_docx() |>
    officer::body_add_par("Table V3. Robustness checks (V-Dem indices)", style="heading 1") |>
    flextable::body_add_flextable(ft_rv)
  print(doc_rv, target = file.path(PATH_OUT, "TableV3_VDem_Robustness.docx"))

} else if (USE_VDEM) {
  message("V-Dem outputs skipped (too few merged V-Dem variables or insufficient N).")
}

# Manifest + Final notes ----
# Manifest of all output files
all_files <- list.files(PATH_OUT, recursive = TRUE, full.names = TRUE)
man <- tibble::tibble(
  file = gsub(paste0("^", PATH_OUT, "/?"), "", all_files),
  bytes = file.info(all_files)$size
) |>
  dplyr::arrange(.data$file)
readr::write_csv(man, file.path(PATH_OUT, "00_manifest_out_files.csv"))

write_lines(c(
  "=== FINAL NOTES ===",
  paste0("Output root: ", PATH_ROOT),
  paste0("Output folder: ", PATH_OUT),
  paste0("Figures: ", FIG_DIR),
  paste0("N (common analysis sample, EIU) = ", N_now),
  paste0("Controls used: ", paste(ctrl_vars, collapse=", ")),
  paste0("V-Dem enabled: ", USE_VDEM),
  "",
  "If N mismatch occurs unexpectedly: check check_missing_merge_whr/qog and EXPECTED_N_ANALYSIS_SAMPLE.txt."
), "analysis_notes.txt")

message("Done. Outputs: ", PATH_OUT, " | Figures: ", FIG_DIR)
