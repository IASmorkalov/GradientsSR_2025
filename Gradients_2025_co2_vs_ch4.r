# Soil CO2 efflux and CH4 uptake along two pollution gradients.

# 0. Configuration

INPUT_FILE <- "SRdata_Revda_Karabash_2025.xlsx"

OUTPUT_DIR <- "results_CO2_CH4_2025"

SEED <- 20251012

T_STEP <- 1.0
H_STEP <- 0.05
T_RADIUS <- 2.0
H_RADIUS <- 0.08
MIN_LOCAL_SITES <- 2L

RF_TREES <- 700L
RF_MTRY <- 2L
RF_MIN_NODE <- 5L
XGB_NROUNDS <- 350L
XGB_ETA <- 0.03
XGB_MAX_DEPTH <- 3L
XGB_MIN_CHILD <- 5
XGB_SUBSAMPLE <- 0.8
XGB_COLSAMPLE <- 0.9

BOOT_B <- 500L

required_packages <- c("tidyverse", "ranger", "xgboost", "readxl", "lubridate")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                               logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Install missing packages first: install.packages(c(",
    paste(sprintf('"%s"', missing_packages), collapse = ", "), "))"
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(ranger)
  library(xgboost)
})

set.seed(SEED)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

area_labels <- c(`1` = "MUCS", `2` = "KCS")
zone_labels <- c(`1` = "background", `2` = "buffer", `3` = "impact")

# 1. Helper functions

mean_or_na <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

# Weighted population SD; all zones contribute to the scaling.
weighted_center_scale <- function(y, w) {
  ok <- is.finite(y) & is.finite(w) & w > 0
  mu <- weighted.mean(y[ok], w[ok])
  s <- sqrt(sum(w[ok] * (y[ok] - mu)^2) / sum(w[ok]))
  if (!is.finite(s) || s <= 0) stop("Response has zero or invalid weighted SD.")
  list(center = mu, scale = s)
}

# Give each site the same total weight.
add_cluster_weights <- function(d, cluster_col = "Site") {
  d %>%
    group_by(.data[[cluster_col]]) %>%
    mutate(.w = 1 / n()) %>%
    ungroup() %>%
    mutate(.w = .w / mean(.w))
}

x_matrix <- function(d) {
  dd <- d %>%
    transmute(
      Zona = factor(as.character(Zona), levels = c("1", "2", "3")),
      Tst = as.numeric(Tst),
      Hums = as.numeric(Hums)
    )
  model.matrix(~ Zona + Tst + Hums - 1, data = dd)
}

fit_ml <- function(d, response, algorithm, cluster_col = "Site", seed = SEED) {
  d2 <- d %>%
    filter(
      is.finite(.data[[response]]), is.finite(Tst), is.finite(Hums),
      Zona %in% c(1, 2, 3)
    ) %>%
    mutate(Zona = factor(as.character(Zona), levels = c("1", "2", "3"))) %>%
    add_cluster_weights(cluster_col)

  sc <- weighted_center_scale(d2[[response]], d2$.w)
  d2 <- d2 %>% mutate(.z = (.data[[response]] - sc$center) / sc$scale)

  if (algorithm == "RF") {

    model <- ranger(
      .z ~ Zona + Tst + Hums,
      data = d2,
      num.trees = RF_TREES,
      mtry = min(RF_MTRY, 3L),
      min.node.size = RF_MIN_NODE,
      case.weights = d2$.w,
      importance = "permutation",
      respect.unordered.factors = "partition",
      seed = seed,
      num.threads = max(1L, parallel::detectCores(logical = FALSE) - 1L)
    )
  } else if (algorithm == "XGB") {

    dm <- xgb.DMatrix(
      data = x_matrix(d2),
      label = d2$.z,
      weight = d2$.w
    )
    model <- xgb.train(
      data = dm,
      nrounds = XGB_NROUNDS,
      params = list(
        objective = "reg:squarederror",
        eval_metric = "rmse",
        eta = XGB_ETA,
        max_depth = XGB_MAX_DEPTH,
        min_child_weight = XGB_MIN_CHILD,
        subsample = XGB_SUBSAMPLE,
        colsample_bytree = XGB_COLSAMPLE,
        nthread = max(1L, parallel::detectCores(logical = FALSE) - 1L),
        seed = seed
      ),
      verbose = 0
    )
  } else {
    stop("Unknown algorithm: ", algorithm)
  }

  structure(
    list(model = model, algorithm = algorithm, center = sc$center, scale = sc$scale),
    class = "gas_ml_fit"
  )
}

predict_ml <- function(object, newdata, scale = c("z", "original")) {
  scale <- match.arg(scale)
  nd <- newdata %>%
    mutate(Zona = factor(as.character(Zona), levels = c("1", "2", "3")))

  pred_z <- if (object$algorithm == "RF") {
    predict(object$model, data = nd)$predictions
  } else {
    as.numeric(predict(object$model, x_matrix(nd)))
  }

  if (scale == "z") pred_z else pred_z * object$scale + object$center
}

metric_rmse <- function(obs, pred) sqrt(mean((obs - pred)^2, na.rm = TRUE))
metric_mae <- function(obs, pred) mean(abs(obs - pred), na.rm = TRUE)

q025 <- function(x) quantile(x, 0.025, na.rm = TRUE, names = FALSE)
q975 <- function(x) quantile(x, 0.975, na.rm = TRUE, names = FALSE)

# 2. Import and aggregate data

raw <- readxl::read_excel(INPUT_FILE, sheet = 1)

required_columns <- c(
  "Area", "Zona", "Site", "Plot", "Date",
  "Rs", "fch4", "Tst", "Tsa", "Hums"
)
missing_columns <- setdiff(required_columns, names(raw))
if (length(missing_columns) > 0) {
  stop("Missing columns in input data: ", paste(missing_columns, collapse = ", "))
}

dat <- raw %>%
  mutate(
    Area = as.integer(Area),
    Zona = as.integer(Zona),
    Plot = as.character(Plot),
    Site = as.character(Site),
    Date = as.character(Date)
  ) %>%
  group_by(Area, Zona, Site, Plot, Date) %>%
  summarise(
    Rs = mean_or_na(Rs),
    fch4 = mean_or_na(fch4),
    Tst = mean_or_na(Tst),
    Tsa = mean_or_na(Tsa),
    Hums = mean_or_na(Hums),
    n_source_rows = n(),
    .groups = "drop"
  ) %>%
  mutate(

    # Preserve the original absolute-flux transformation.
    CH4_uptake = abs(fch4),
    Area_name = recode(as.character(Area), !!!area_labels),
    Zona_name = recode(as.character(Zona), !!!zone_labels)
  ) %>%
  filter(if_all(c(Rs, CH4_uptake, Tst, Hums), is.finite))

if (nrow(dat) == 0) stop("No complete observations after aggregation.")
if (any(dat$fch4 > 0, na.rm = TRUE)) {
  warning("Positive fch4 values were found. CH4_uptake = -fch4 will be negative for them.")
}
if (!all(c(1, 2) %in% unique(dat$Area))) stop("Both Area 1 and Area 2 are required.")
if (!all(c(1, 3) %in% unique(dat$Zona))) stop("Background and impact zones are required.")

readr::write_csv(dat, file.path(OUTPUT_DIR, "data_2025_aggregated.csv"))

coverage <- dat %>%
  count(Area_name, Zona_name, Site, name = "n_plot_dates")
readr::write_csv(coverage, file.path(OUTPUT_DIR, "sampling_coverage.csv"))

# 3. Leave-one-site-out cross-validation

# Keep all plots and dates from a site in the same validation fold.
run_loso <- function(d_area, response, algorithm) {
  sites <- sort(unique(d_area$Site))

  map_dfr(seq_along(sites), function(k) {
    test_site <- sites[k]
    train <- d_area %>% filter(Site != test_site)
    test <- d_area %>% filter(Site == test_site)

    if (!all(c(1, 3) %in% as.integer(as.character(unique(train$Zona))))) {
      return(tibble())
    }

    fit <- fit_ml(
      train, response = response, algorithm = algorithm,
      cluster_col = "Site", seed = SEED + k
    )

    tibble(
      Area = unique(d_area$Area),
      Area_name = unique(d_area$Area_name),
      response = response,
      Gas = if_else(response == "Rs", "CO2", "CH4 uptake"),
      algorithm = algorithm,
      test_site = test_site,
      observed_z = (test[[response]] - fit$center) / fit$scale,
      predicted_z = predict_ml(fit, test, scale = "z")
    )
  })
}

cv_predictions <- crossing(
  Area = sort(unique(dat$Area)),
  response = c("Rs", "CH4_uptake"),
  algorithm = c("RF", "XGB")
) %>%
  pmap_dfr(function(Area, response, algorithm) {
    area_i <- Area
    message("LOSO-CV: Area=", area_i, ", response=", response, ", model=", algorithm)
    run_loso(dat %>% filter(.data$Area == area_i), response, algorithm)
  })

if (nrow(cv_predictions) == 0) stop("LOSO-CV produced no predictions.")
readr::write_csv(cv_predictions, file.path(OUTPUT_DIR, "cv_predictions.csv"))

cv_metrics_by_task <- cv_predictions %>%
  group_by(Area, Area_name, response, Gas, algorithm) %>%
  summarise(
    n_sites = n_distinct(test_site),
    n = n(),
    RMSE_z = metric_rmse(observed_z, predicted_z),
    MAE_z = metric_mae(observed_z, predicted_z),
    bias_z = mean(predicted_z - observed_z),
    .groups = "drop"
  )

model_selection <- cv_metrics_by_task %>%
  group_by(algorithm) %>%
  summarise(
    mean_RMSE_z = mean(RMSE_z),
    mean_MAE_z = mean(MAE_z),
    max_RMSE_z = max(RMSE_z),
    .groups = "drop"
  ) %>%
  arrange(mean_RMSE_z, mean_MAE_z)

chosen_algorithm <- model_selection$algorithm[1]
message("Chosen common algorithm: ", chosen_algorithm)

readr::write_csv(cv_metrics_by_task, file.path(OUTPUT_DIR, "cv_metrics_by_task.csv"))
readr::write_csv(model_selection, file.path(OUTPUT_DIR, "cv_model_selection.csv"))
writeLines(chosen_algorithm, file.path(OUTPUT_DIR, "chosen_algorithm.txt"))

# 4. Supported temperature-moisture grid

# Retain nodes locally supported by at least two sites in each comparison zone.
make_supported_grid <- function(d_area) {
  bg <- d_area %>% filter(Zona == 1)
  imp <- d_area %>% filter(Zona == 3)

  t_low <- max(min(bg$Tst), min(imp$Tst))
  t_high <- min(max(bg$Tst), max(imp$Tst))
  h_low <- max(min(bg$Hums), min(imp$Hums))
  h_high <- min(max(bg$Hums), max(imp$Hums))

  t_seq <- seq(ceiling(t_low / T_STEP) * T_STEP,
               floor(t_high / T_STEP) * T_STEP, by = T_STEP)
  h_seq <- seq(ceiling(h_low / H_STEP) * H_STEP,
               floor(h_high / H_STEP) * H_STEP, by = H_STEP)

  if (length(t_seq) == 0 || length(h_seq) == 0) {
    stop("No overlapping rectangular T-H range for Area ", unique(d_area$Area))
  }

  grid <- tidyr::crossing(Tst = t_seq, Hums = h_seq)

  count_local_sites <- function(zone_code) {
    dz <- d_area %>% filter(Zona == zone_code)
    map2_int(grid$Tst, grid$Hums, function(tt, hh) {
      dz %>%
        filter(abs(Tst - tt) <= T_RADIUS, abs(Hums - hh) <= H_RADIUS) %>%
        summarise(n = n_distinct(Site)) %>%
        pull(n)
    })
  }

  grid %>%
    mutate(
      n_sites_background = count_local_sites(1),
      n_sites_impact = count_local_sites(3),

      supported = n_sites_background >= MIN_LOCAL_SITES &
        n_sites_impact >= MIN_LOCAL_SITES
    ) %>%
    filter(supported) %>%
    mutate(
      Area = unique(d_area$Area),
      Area_name = unique(d_area$Area_name)
    )
}

supported_grids <- map_dfr(sort(unique(dat$Area)), function(a) {
  g <- make_supported_grid(dat %>% filter(Area == a))
  if (nrow(g) == 0) {
    stop(
      "No supported T-H grid points for Area ", a,
      ". Reduce MIN_LOCAL_SITES or enlarge T_RADIUS/H_RADIUS after inspection."
    )
  }
  g
})
readr::write_csv(supported_grids, file.path(OUTPUT_DIR, "supported_T_H_grid.csv"))

# 5. Final models and pollution contrasts
response_key <- tribble(
  ~response,       ~Gas,
  "Rs",            "CO2",
  "CH4_uptake",    "CH4 uptake"
)

final_models <- list()
contrast_surfaces <- list()

for (a in sort(unique(dat$Area))) {
  da <- dat %>% filter(Area == a)
  ga <- supported_grids %>% filter(Area == a)

  for (r in response_key$response) {
    gas_name <- response_key$Gas[match(r, response_key$response)]
    fit <- fit_ml(da, r, chosen_algorithm, cluster_col = "Site", seed = SEED + a)
    final_models[[paste(a, r, sep = "_")]] <- fit

    pred_bg <- ga %>% mutate(Zona = 1L)
    pred_imp <- ga %>% mutate(Zona = 3L)

    contrast_surfaces[[paste(a, r, sep = "_")]] <- ga %>%
      transmute(
        Area, Area_name, Tst, Hums,
        n_sites_background, n_sites_impact,
        response = r,
        Gas = gas_name,
        pred_background_z = predict_ml(fit, pred_bg, "z"),
        pred_impact_z = predict_ml(fit, pred_imp, "z"),

        contrast_z = pred_impact_z - pred_background_z
      )
  }
}

saveRDS(final_models, file.path(OUTPUT_DIR, "final_models.rds"))
contrasts <- bind_rows(contrast_surfaces)
readr::write_csv(contrasts, file.path(OUTPUT_DIR, "contrast_surfaces.csv"))

# Positive differences indicate larger absolute standardized CH4 contrasts.
delta_surface <- contrasts %>%
  select(Area, Area_name, Tst, Hums, Gas, contrast_z) %>%
  pivot_wider(names_from = Gas, values_from = contrast_z) %>%
  mutate(
    delta_CH4_minus_CO2 = abs(`CH4 uptake`) - abs(CO2),
    stronger_indicator = case_when(
      delta_CH4_minus_CO2 > 0 ~ "CH4 uptake",
      delta_CH4_minus_CO2 < 0 ~ "CO2",
      TRUE ~ "equal"
    )
  )
readr::write_csv(delta_surface, file.path(OUTPUT_DIR, "delta_CH4_minus_CO2_surface.csv"))

point_gas_summary <- contrasts %>%
  group_by(Area, Area_name, Gas) %>%
  summarise(
    median_abs_contrast = median(abs(contrast_z)),
    median_signed_contrast = median(contrast_z),
    sign_consistency = max(mean(contrast_z > 0), mean(contrast_z < 0)),
    n_grid = n(),
    .groups = "drop"
  )

point_winner_summary <- delta_surface %>%
  group_by(Area, Area_name) %>%
  summarise(
    median_delta_CH4_minus_CO2 = median(delta_CH4_minus_CO2),
    prop_grid_CH4_stronger = mean(delta_CH4_minus_CO2 > 0),
    prop_grid_CO2_stronger = mean(delta_CH4_minus_CO2 < 0),
    n_grid = n(),
    .groups = "drop"
  )

# 6. Site-cluster bootstrap

# Resample whole sites without zone stratification; retain all nested records.
make_boot_sample <- function(d_area) {
  sites <- unique(d_area$Site)
  sampled <- sample(sites, length(sites), replace = TRUE)

  map2_dfr(sampled, seq_along(sampled), function(s, copy_id) {
    d_area %>%
      filter(Site == s) %>%
      mutate(.boot_cluster = paste0("copy_", copy_id, "_", s))
  })
}

# Refit both gases on paired bootstrap samples; keep the original grid fixed.
bootstrap_gas_metrics <- vector("list", BOOT_B * length(unique(dat$Area)))
bootstrap_delta_metrics <- vector("list", BOOT_B * length(unique(dat$Area)))
bootstrap_delta_surface <- vector("list", BOOT_B * length(unique(dat$Area)))
idx <- 0L

for (b in seq_len(BOOT_B)) {
  if (b == 1L || b %% 25L == 0L) message("Bootstrap ", b, "/", BOOT_B)

  for (a in sort(unique(dat$Area))) {
    idx <- idx + 1L
    da <- dat %>% filter(Area == a)
    ga <- supported_grids %>% filter(Area == a)
    db <- make_boot_sample(da)

    if (!all(c(1, 3) %in% unique(db$Zona))) next

    one_rep <- tryCatch({

      fit_co2 <- fit_ml(
        db, "Rs", chosen_algorithm,
        cluster_col = ".boot_cluster", seed = SEED + 10000L + b + a
      )
      fit_ch4 <- fit_ml(
        db, "CH4_uptake", chosen_algorithm,
        cluster_col = ".boot_cluster", seed = SEED + 20000L + b + a
      )

      bg <- ga %>% mutate(Zona = 1L)
      imp <- ga %>% mutate(Zona = 3L)
      c_co2 <- predict_ml(fit_co2, imp, "z") - predict_ml(fit_co2, bg, "z")
      c_ch4 <- predict_ml(fit_ch4, imp, "z") - predict_ml(fit_ch4, bg, "z")

      delta <- abs(c_ch4) - abs(c_co2)

      gas_metrics <- tibble(
        iteration = b,
        Area = a,
        Area_name = unique(da$Area_name),
        Gas = c("CO2", "CH4 uptake"),
        median_abs_contrast = c(median(abs(c_co2)), median(abs(c_ch4))),
        median_signed_contrast = c(median(c_co2), median(c_ch4)),
        sign_consistency = c(
          max(mean(c_co2 > 0), mean(c_co2 < 0)),
          max(mean(c_ch4 > 0), mean(c_ch4 < 0))
        )
      )

      delta_metrics <- tibble(
        iteration = b,
        Area = a,
        Area_name = unique(da$Area_name),
        median_delta_CH4_minus_CO2 = median(delta),
        prop_grid_CH4_stronger = mean(delta > 0),
        prop_grid_CO2_stronger = mean(delta < 0)
      )

      delta_grid <- ga %>%
        transmute(
          iteration = b, Area, Area_name, Tst, Hums,
          contrast_CO2_z = c_co2,
          contrast_CH4_z = c_ch4,
          delta_CH4_minus_CO2 = delta
        )

      list(gas = gas_metrics, delta = delta_metrics, grid = delta_grid)
    }, error = function(e) {
      warning("Bootstrap failed: iteration=", b, ", Area=", a, ": ", e$message)
      NULL
    })

    if (!is.null(one_rep)) {
      bootstrap_gas_metrics[[idx]] <- one_rep$gas
      bootstrap_delta_metrics[[idx]] <- one_rep$delta
      bootstrap_delta_surface[[idx]] <- one_rep$grid
    }
  }
}

boot_gas <- bind_rows(bootstrap_gas_metrics)
boot_delta <- bind_rows(bootstrap_delta_metrics)
boot_grid <- bind_rows(bootstrap_delta_surface)

if (nrow(boot_gas) == 0 || nrow(boot_delta) == 0) {
  stop("All bootstrap replicates failed. Inspect zone representation and model settings.")
}

readr::write_csv(boot_gas, file.path(OUTPUT_DIR, "bootstrap_gas_metrics.csv"))
readr::write_csv(boot_delta, file.path(OUTPUT_DIR, "bootstrap_delta_metrics.csv"))

boot_gas_ci <- boot_gas %>%
  group_by(Area, Area_name, Gas) %>%
  summarise(
    valid_boot = n(),
    median_abs_boot = median(median_abs_contrast),
    median_abs_low = q025(median_abs_contrast),
    median_abs_high = q975(median_abs_contrast),
    median_signed_boot = median(median_signed_contrast),
    median_signed_low = q025(median_signed_contrast),
    median_signed_high = q975(median_signed_contrast),
    sign_consistency_boot = median(sign_consistency),
    sign_consistency_low = q025(sign_consistency),
    sign_consistency_high = q975(sign_consistency),
    .groups = "drop"
  )

indicator_summary <- point_gas_summary %>%
  left_join(boot_gas_ci, by = c("Area", "Area_name", "Gas"))
readr::write_csv(indicator_summary, file.path(OUTPUT_DIR, "indicator_summary.csv"))

boot_delta_ci <- boot_delta %>%
  group_by(Area, Area_name) %>%
  summarise(
    valid_boot = n(),
    median_delta_boot = median(median_delta_CH4_minus_CO2),
    median_delta_low = q025(median_delta_CH4_minus_CO2),
    median_delta_high = q975(median_delta_CH4_minus_CO2),
    prop_CH4_stronger_boot = median(prop_grid_CH4_stronger),
    prop_CH4_stronger_low = q025(prop_grid_CH4_stronger),
    prop_CH4_stronger_high = q975(prop_grid_CH4_stronger),
    probability_CH4_median_effect_larger = mean(median_delta_CH4_minus_CO2 > 0),
    .groups = "drop"
  )

winner_summary <- point_winner_summary %>%
  left_join(boot_delta_ci, by = c("Area", "Area_name")) %>%
  mutate(

    result = case_when(
      median_delta_low > 0 ~ "CH4 uptake is the stronger indicator",
      median_delta_high < 0 ~ "CO2 is the stronger indicator",
      TRUE ~ "No clear overall winner"
    )
  )
readr::write_csv(winner_summary, file.path(OUTPUT_DIR, "winner_summary.csv"))

boot_grid_ci <- boot_grid %>%
  group_by(Area, Area_name, Tst, Hums) %>%
  summarise(
    valid_boot = n(),
    delta_median = median(delta_CH4_minus_CO2),
    delta_low = q025(delta_CH4_minus_CO2),
    delta_high = q975(delta_CH4_minus_CO2),
    p_CH4_stronger = mean(delta_CH4_minus_CO2 > 0),
    certain_winner = case_when(
      delta_low > 0 ~ "CH4 uptake",
      delta_high < 0 ~ "CO2",
      TRUE ~ "uncertain"
    ),
    .groups = "drop"
  )
readr::write_csv(boot_grid_ci, file.path(OUTPUT_DIR, "delta_surface_bootstrap_CI.csv"))

# 7. Model-based figures

contrast_limits <- max(abs(contrasts$contrast_z), na.rm = TRUE)

p_contrasts <- ggplot(contrasts, aes(Tst, Hums, fill = contrast_z)) +
  geom_tile() +
  facet_grid(Gas ~ Area_name) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    limits = c(-contrast_limits, contrast_limits),
    name = "Impact - background\n(SD units)"
  ) +
  labs(
    x = "Soil temperature, °C", y = expression("Soil moisture, " * m^3 / m^3),
    title = "Pollution contrasts under supported hydrothermal conditions"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill = "grey92"))

ggsave(
  file.path(OUTPUT_DIR, "Fig1_pollution_contrasts.png"),
  p_contrasts, width = 10, height = 7, dpi = 320
)

delta_limits <- max(abs(boot_grid_ci$delta_median), na.rm = TRUE)
p_delta <- ggplot(boot_grid_ci, aes(Tst, Hums, fill = delta_median)) +
  geom_tile() +
  geom_point(
    data = boot_grid_ci %>% filter(certain_winner == "uncertain"),
    shape = 4, size = 0.9, alpha = 0.7
  ) +
  facet_wrap(~ Area_name) +
  scale_fill_gradient2(
    low = "#D95F02", mid = "white", high = "#1B9E77", midpoint = 0,
    limits = c(-delta_limits, delta_limits),
    name = "|CH4 effect| - |CO2 effect|\n(SD units)"
  ) +
  labs(
    x = "Soil temperature, °C", y = expression("Soil moisture, " * m^3 / m^3),
    title = "Which gas is the stronger pollution indicator?",
    subtitle = "Green: CH4 uptake; orange: CO2; crosses: 95% bootstrap interval includes zero"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(), strip.background = element_rect(fill = "grey92"))

ggsave(
  file.path(OUTPUT_DIR, "Fig2_indicator_advantage.png"),
  p_delta, width = 10, height = 4.8, dpi = 320
)

p_summary <- indicator_summary %>%
  ggplot(
    aes(
      x = median_abs_contrast,
      y = interaction(Gas, Area_name),
      color = Gas,
      shape=Gas
    )
  ) +
  geom_errorbar(
    aes(
      xmin = median_abs_low,
      xmax = median_abs_high
    ),
    orientation = "y",
    width = 0.16,
    linewidth = 0.7
  ) +
  geom_point(size = 2.8) +
  scale_color_manual(
    values = c(
      "CO2" = "#454545",
      "CH4 uptake"   = "#0072B2"
    )
  ) +
  labs(
    x = "Median absolute pollution contrast (SD units)",
    y = NULL,
    title = "Overall indicator strength",
    subtitle = "Points: final models; intervals: 95% site-cluster bootstrap"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

p_summary

ggsave(
  file.path(OUTPUT_DIR, "Fig3_indicator_strength.png"),
  p_summary, width = 8.5, height = 4.5, dpi = 320
)

# 8. Reproducibility

writeLines(capture.output(sessionInfo()), file.path(OUTPUT_DIR, "sessionInfo.txt"))

message("Analysis completed. Results are in: ", normalizePath(OUTPUT_DIR))
message("Main answer: results_CO2_CH4_2025/winner_summary.csv")
message("Main map: results_CO2_CH4_2025/Fig2_indicator_advantage.png")

library(lubridate)
# Accept Excel-imported ISO dates and day-month-year text dates.
dat$Date <- as.Date(lubridate::parse_date_time(
  dat$Date, orders = c("ymd", "dmy", "ymd HMS", "dmy HMS"), quiet = TRUE
))
if (anyNA(dat$Date)) stop("Some dates could not be parsed. Check the Date column.")

# 9. Descriptive seasonal dynamics

library(dplyr)
library(tidyr)
library(ggplot2)

mean_safe <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

if (!inherits(dat$Date, "Date")) {
  stop(
    "Convert dat$Date to Date first, ",
    "using the date format in the input workbook."
  )
}

dat_long <- dat %>%
  mutate(
    Zona_name = factor(
      tolower(trimws(as.character(Zona_name))),
      levels = c("background", "buffer", "impact")
    ),
    Area_name = factor(
      Area_name,
      levels = c("MUCS", "KCS")
    )
  ) %>%
  pivot_longer(
    cols = c(Rs, CH4_uptake),
    names_to = "Gas_name",
    values_to = "GAS"
  ) %>%
  mutate(
    Gas_name = factor(
      Gas_name,
      levels = c("Rs", "CH4_uptake")
    )
  ) %>%
  filter(
    !is.na(Date),
    !is.na(Area_name),
    !is.na(Zona_name),
    is.finite(GAS)
  )

# Average plots within sites before calculating zone means and among-site SD.
gas_site_date <- dat_long %>%
  group_by(Area_name, Zona_name, Gas_name, Date, Site) %>%
  summarise(
    GAS = mean_safe(GAS),
    n_plots = n_distinct(Plot),
    .groups = "drop"
  )

gas_date_summary <- gas_site_date %>%
  group_by(Area_name, Zona_name, Gas_name, Date) %>%
  summarise(
    n_sites = sum(is.finite(GAS)),
    GAS_sd = sd(GAS, na.rm = TRUE),
    GAS = mean_safe(GAS),
    .groups = "drop"
  ) %>%
  mutate(
    GAS_low = GAS - GAS_sd,
    GAS_high = GAS + GAS_sd
  )

zone_colors <- c(
  "background" = "#228B63",
  "buffer"     = "#C49A00",
  "impact"     = "#C43C39"
)

gas_labels <- c(
  "Rs" = "CO₂ emission (µmol m⁻² s⁻¹)",
  "CH4_uptake" = "CH₄ uptake (nmol m⁻² s⁻¹)"
)

smooth_data <- gas_date_summary %>%
  group_by(Area_name, Zona_name, Gas_name) %>%
  filter(n_distinct(Date) >= 5L) %>%
  ungroup()

p_dynamics <- ggplot(
  gas_date_summary,
  aes(
    x = Date,
    y = GAS,
    color = Zona_name,
    group = Zona_name
  )
) +
  geom_hline(
    yintercept = 0,
    color = "grey75",
    linewidth = 0.35,
    linetype = "dashed"
  ) +
  geom_errorbar(
    aes(ymin = GAS_low, ymax = GAS_high),
    width = 2,
    linewidth = 0.45,
    alpha = 0.65,
    na.rm = TRUE
  ) +
  geom_smooth(
    data = smooth_data,
    method = "loess",

    se = FALSE,
    linewidth = 0.9
  ) +
  geom_point(
    size = 2.2,
    alpha = 0.9
  ) +
  facet_wrap(
    vars(Gas_name, Area_name),
    ncol = 2,
    scales = "free_y",
    labeller = labeller(Gas_name = gas_labels)
  ) +
  scale_color_manual(
    values = zone_colors,
    breaks = c("background", "buffer", "impact"),
    labels = c("Background", "Buffer", "Impact"),
    drop = FALSE
  ) +
  scale_x_date(
    date_breaks = "1 month",
    date_labels = "%m",
    expand = expansion(mult = c(0.025, 0.025))
  ) +
  labs(
    x = "Month",
    y = NULL,
    color = "Pollution zone",
    caption = paste(
      "Points: means across sites; error bars: ±1 SD among site means.",
      "Lines: descriptive cubic regression splines."
    )
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(
      fill = "grey95",
      color = "grey70"
    ),
    strip.text = element_text(size = 10),
    axis.text = element_text(color = "black"),
    panel.spacing = grid::unit(1, "lines"),
    plot.caption = element_text(
      hjust = 0,
      size = 9,
      color = "grey35"
    )
  )

p_dynamics

dir.create(
  "results_CO2_CH4_2025",
  showWarnings = FALSE,
  recursive = TRUE
)

ggsave(
  filename = "results_CO2_CH4_2025/Fig4_field_dynamics.png",
  plot = p_dynamics,
  width = 11,
  height = 7,
  units = "in",
  dpi = 600,
  bg = "white"
)

write.csv(
  gas_date_summary,
  file = "results_CO2_CH4_2025/field_dynamics_summary.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
