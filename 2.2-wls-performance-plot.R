suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(broom)
  library(purrr); library(ggplot2)
})

# ======================
# 0) Load & helpers
# ======================
data_path <- "/Users/yuchang/PhD/PhD_research/swiss_heat_exposure/data/revision_r1/muni_all_vars_final.csv"
dat <- read_csv(data_path, show_col_types = FALSE)

z    <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
mm01 <- function(x){r <- range(x, na.rm = TRUE); (x - r[1])/(r[2]-r[1])}
AICc <- function(mod){
  k <- length(coef(mod))                # number of parameters incl. intercept
  a <- AIC(mod)
  n <- length(residuals(mod))
  a + (2 * k * (k + 1)) / pmax(n - k - 1, 1)  # guard small denom
}

# ======================
# 1) Prepare variables
# ======================
raw_needed <- c(
  # outcomes
  "lst_rbn","avgt25_","htwv_p_",
  # physical drivers
  "dem_men","imprv_m","evi_men","bhght_m","albd_mn",
  # socio-demographics
  "income","age65_","age80_","forgn19","for_eu","socisss","av_lvng",
  # weights
  "num_rsd"
)

dat <- dat %>%
  mutate(across(where(is.character), ~na_if(.,"X"))) %>%
  mutate(across(intersect(names(.), raw_needed), ~suppressWarnings(as.numeric(.)))) %>%
  mutate(
    for_eu_prop  = ifelse(for_eu>1, for_eu/100, for_eu),
    age65_79     = age65_ - age80_,
    percen_eu    = forgn19 * for_eu_prop,
    percen_noneu = forgn19 * (1 - for_eu_prop),
    w            = log(pmax(num_rsd,0)+1),
    income_ln    = log(pmax(income,0)+1),
    LST_z = z(lst_rbn),
    HWD_z = z(avgt25_),
    HWP_z = z(htwv_p_)
  ) %>%
  mutate(CHEI_equal = mm01(as.numeric(as.matrix(select(., LST_z,HWD_z,HWP_z)) %*% rep(1/3,3))))

# Predictor groups (names fixed)
physical_vars <- c("dem_men","evi_men","bhght_m","albd_mn","imprv_m")
socio_vars    <- c("income_ln","age65_79","age80_","percen_eu","percen_noneu","socisss","av_lvng")
socio_vars_no <- setdiff(socio_vars, "sociss")   # "no_socio" version

# Drop impossible values and keep weights
dat <- dat %>%
  mutate(across(c("CHEI_equal",physical_vars,socio_vars,"num_rsd"), ~ifelse(is.infinite(.),NA,.))) %>%
  mutate(w = log(num_rsd + 1))

# Create a z-scored copy for modeling
dat_step <- dat %>% mutate(across(c(all_of(physical_vars), all_of(socio_vars)), z))

# ======================
# 2) Runner for models
# ======================
fit_three_models <- function(df, dv = "CHEI_equal", socio_set, run_tag){
  # m1: socio only
  f1 <- reformulate(socio_set, response = dv)
  # m2: + dem_men (geo)
  f2 <- reformulate(c(socio_set, "dem_men"), response = dv)
  # m3: + landscape (evi/bldg/albd/impervious)
  f3 <- reformulate(c(socio_set, "dem_men","evi_men","bhght_m","albd_mn","imprv_m"), response = dv)
  
  # ensure complete cases for each formula
  d1 <- df %>% drop_na(all_of(all.vars(f1)))
  d2 <- df %>% drop_na(all_of(all.vars(f2)))
  d3 <- df %>% drop_na(all_of(all.vars(f3)))
  
  m1 <- lm(f1, data = d1, weights = w)
  m2 <- lm(f2, data = d2, weights = w)
  m3 <- lm(f3, data = d3, weights = w)
  
  diag_one <- function(mod, data_used, model_tag){
    gl <- glance(mod)
    pred <- predict(mod)
    obs  <- mod$model[[dv]]
    tibble(
      run      = run_tag,
      model    = model_tag,
      dep_var  = dv,
      n        = nrow(data_used),
      R2       = gl$r.squared,
      Adj_R2   = gl$adj.r.squared,
      RMSE     = sqrt(mean((obs - pred)^2, na.rm = TRUE)),
      MAE      = mean(abs(obs - pred), na.rm = TRUE),
      F_stat   = gl$statistic,
      F_p      = gl$p.value,
      AIC      = AIC(mod),
      AICc     = AICc(mod),
      BIC      = BIC(mod)
    )
  }
  
  bind_rows(
    diag_one(m1, d1, "m1_socio"),
    diag_one(m2, d2, "m2_+geo"),
    diag_one(m3, d3, "m3_+geo+urban")
  )
}

# ======================
# 3) Run WITH and NO socio
# ======================
diag_with    <- fit_three_models(dat_step, dv = "CHEI_equal", socio_set = socio_vars,    run_tag = "with_socio")

ols_diag_all <- bind_rows(diag_with) %>%
  arrange(run, model)

print(ols_diag_all)

# Save if you want
write_csv(ols_diag_all, "/Users/yuchang/PhD/PhD_research/swiss_heat_exposure/data/revision_r1/ols_wls_diagnostics_CHEI.csv")

# ================================
# Faceted stepwise predictions (uses m1/m2/m3)
# ================================
socio_labels <- c(
  income_ln   = "Income (log)",
  age65_79    = "% Age 65–79",
  age80_      = "% Age 80+",
  percen_eu   = "% EU foreigners",
  percen_noneu= "% Non–EU foreigners",
  socioassis  = "Social assistance",
  ave_living  = "Avg. living space"
)

means_row <- dat_step %>%
  summarise(across(c(all_of(socio_vars), dem_men,evi_men,bhght_m,albd_mn,imprv_m), ~mean(.x, na.rm = TRUE)))

make_grid <- function(var, n = 100) {
  rng <- range(dat_step[[var]], na.rm = TRUE)
  base <- means_row[rep(1, n), ]
  base[[var]] <- seq(rng[1], rng[2], length.out = n)
  base$var <- var; base$x <- base[[var]]; base
}
pred_grid <- purrr::map_dfr(socio_vars, make_grid)

pred_long <- pred_grid %>%
  group_by(var) %>% tidyr::nest() %>%
  mutate(
    pred_m1 = purrr::map(data, ~ as_tibble(predict(m1, newdata = .x, se.fit = TRUE)) %>% mutate(x = .x$x)),
    pred_m2 = purrr::map(data, ~ as_tibble(predict(m2, newdata = .x, se.fit = TRUE)) %>% mutate(x = .x$x)),
    pred_m3 = purrr::map(data, ~ as_tibble(predict(m3, newdata = .x, se.fit = TRUE)) %>% mutate(x = .x$x))
  ) %>%
  select(-data) %>%
  pivot_longer(starts_with("pred_"), names_to = "Model", values_to = "pred") %>%
  mutate(Model = dplyr::recode(Model,
                               pred_m1 = "Model 1 (socio)",
                               pred_m2 = "Model 2 (+geo)",
                               pred_m3 = "Model 3 (+geo+urban)")) %>%
  unnest(pred) %>%
  rename(Predicted = fit, se = se.fit) %>%
  mutate(
    lwr = Predicted - 1.96 * se,
    upr = Predicted + 1.96 * se,
    var_lab = factor(var, levels = socio_vars, labels = socio_labels[socio_vars])
  ) %>%
  ungroup()

means_for_vline <- dat_step %>%
  summarise(across(all_of(socio_vars), ~mean(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "var", values_to = "xmean") %>%
  mutate(var_lab = factor(var, levels = socio_vars, labels = socio_labels[socio_vars]))

# 1) Build labels (β and p) – reuse your code up to `anno_df <- ...`
get_beta_p <- function(model, model_label) {
  broom::tidy(model) %>%
    dplyr::filter(term %in% socio_vars) %>%
    dplyr::transmute(
      var      = term,
      var_lab  = factor(term, levels = socio_vars, labels = socio_labels[socio_vars]),
      Model    = model_label,
      beta     = estimate,
      p.value,
      label    = paste0("coefficient = ", round(beta, 3), ", ",
                        ifelse(is.na(p.value), "p = —",
                               ifelse(p.value < 0.001, "p < 0.001",
                                      paste0("p = ", round(p.value, 3)))))
    )
}

anno_df <- dplyr::bind_rows(
  get_beta_p(m1, "Model 1 (socio)"),
  get_beta_p(m2, "Model 2 (+geo)"),
  get_beta_p(m3, "Model 3 (+geo+landscape)")
)

# 2) Stack the three model labels inside each facet AT A FIXED CORNER.
#    Choose corner by setting CORNER <- "tl" (top-left) or "tr" (top-right)
CORNER <- "tl"  # change to "tr" for top-right

order_lvls <- c("Model 1 (socio)", "Model 2 (+geo)", "Model 3 (+geo+landscape)")
anno_df <- anno_df %>%
  dplyr::mutate(
    row   = as.integer(factor(Model, levels = order_lvls)),
    x_inf = ifelse(CORNER == "tr", Inf, -Inf),
    y_inf = Inf,
    # pad from the corner (hjust/vjust are scale-free)
    hjust = ifelse(CORNER == "tr", 1.02, -0.02),      # small inward offset
    vjust = c(1.50, 3.00, 4.50)[row]   # now 1.40 apart                # vertical stacking
  )

p <- ggplot(pred_long, aes(x = x, y = Predicted, color = Model, fill = Model)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, color = NA, show.legend = TRUE) +
  geom_line(linewidth = 0.8) +
  geom_vline(data = means_for_vline, aes(xintercept = xmean),
             linetype = "dashed", size = 0.4, color = "grey50") +
  geom_text(data = anno_df,
            aes(x = x_inf, y = y_inf, label = label, color = Model, hjust = hjust, vjust = vjust),
            inherit.aes = FALSE, size = 3) +
  facet_wrap(~ var_lab, scales = "free_x", ncol = 4) +
  labs(title = NULL,
       y = "Composite heat exposure index", x = "Standardized socio-demographic variables") +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    
    # >>> place legend in the white gap to the right between rows
    legend.position = c(0.76, 0.20),     # (x, y) in [0,1] of the whole plot
    legend.justification = c(0, 0.5),    # anchor legend's left-middle to that point
    legend.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(t = 6, r = 18, b = 6, l = 6) # a little right margin
  ) +
  guides(color = guide_legend(override.aes = list(linetype = 1, linewidth = 1)))
p




suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(broom)
  library(purrr); library(ggplot2); library(tibble)
})

# ======================
# 0) Paths & helpers
# ======================
data_path <- "/Users/yuchang/PhD/PhD_research/swiss_heat_exposure/data/revision_r1/muni_all_vars_final.csv"

z     <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
mm01  <- function(x){r <- range(x, na.rm = TRUE); (x - r[1])/(r[2]-r[1])}
AICc_ <- function(mod){
  k <- length(coef(mod)); a <- AIC(mod); n <- length(residuals(mod))
  a + (2 * k * (k + 1)) / pmax(n - k - 1, 1)
}

# ======================
# 1) Load & preprocess
# ======================
raw_needed <- c(
  # outcomes / inputs to CHEI
  "lst_rbn","avgt25_","htwv_p_",
  # physical drivers
  "dem_men","imprv_m","evi_men","bhght_m","albd_mn",
  # socio-demographics
  "income","age65_","age80_","forgn19","for_eu","socisss","av_lvng",
  # weights
  "num_rsd"
)

dat <- read_csv(data_path, show_col_types = FALSE) %>%
  mutate(across(where(is.character), ~na_if(.,"X"))) %>%
  mutate(across(all_of(raw_needed), ~suppressWarnings(as.numeric(.)))) %>%
  mutate(
    # splits & transforms
    for_eu_prop  = ifelse(for_eu > 1, for_eu/100, for_eu),
    age65_79     = age65_ - age80_,
    percen_eu    = forgn19 * for_eu_prop,
    percen_noneu = forgn19 * (1 - for_eu_prop),
    income_ln    = log(pmax(income,0)+1),
    # weights: log(population)
    w            = log(pmax(num_rsd, 1)),
    # z of heat indicators
    LST_z = z(lst_rbn),
    HWD_z = z(avgt25_),
    HWP_z = z(htwv_p_)
  ) %>%
  # Equal-weight CHEI (main DV here)
  mutate(CHEI_equal = mm01(as.numeric(as.matrix(select(., LST_z,HWD_z,HWP_z)) %*% rep(1/3,3)))) %>%
  # clean infs
  mutate(across(c("CHEI_equal", raw_needed), ~ifelse(is.infinite(.), NA, .)))

# ======================
# 2) Define variable sets & z-standardize predictors
# ======================
physical_vars <- c("dem_men","evi_men","bhght_m","albd_mn","imprv_m")
socio_vars    <- c("income_ln","age65_79","age80_","percen_eu","percen_noneu","socisss","av_lvng")

# z-scored predictors for modeling; keep weight
dat_step <- dat %>%
  mutate(across(c(all_of(physical_vars), all_of(socio_vars)), z)) %>%
  # ensure weight present & finite
  mutate(w = ifelse(is.infinite(w) | is.na(w), NA, w))

# ======================
# 3) Fit three WLS models (M1 socio → M2 +geo → M3 +geo+urban)
# ======================
DV <- "CHEI_equal"

f1 <- reformulate(socio_vars, response = DV)
f2 <- reformulate(c(socio_vars, "dem_men"), response = DV)
f3 <- reformulate(c(socio_vars, "dem_men","evi_men","bhght_m","albd_mn","imprv_m"), response = DV)

d1 <- dat_step %>% drop_na(all_of(all.vars(f1)), w)
d2 <- dat_step %>% drop_na(all_of(all.vars(f2)), w)
d3 <- dat_step %>% drop_na(all_of(all.vars(f3)), w)

m1 <- lm(f1, data = d1, weights = w)
m2 <- lm(f2, data = d2, weights = w)
m3 <- lm(f3, data = d3, weights = w)

diag_tbl <- function(m, d, tag){
  gl <- glance(m); pred <- predict(m); obs <- m$model[[DV]]
  tibble(
    Model = tag, DV = DV, n = nrow(d),
    R2 = gl$r.squared, Adj_R2 = gl$adj.r.squared,
    RMSE = sqrt(mean((obs - pred)^2, na.rm = TRUE)),
    MAE  = mean(abs(obs - pred), na.rm = TRUE),
    AIC  = AIC(m), AICc = AICc_(m), BIC = BIC(m)
  )
}

diagnostics <- bind_rows(
  diag_tbl(m1, d1, "M1: socio"),
  diag_tbl(m2, d2, "M2: +geo"),
  diag_tbl(m3, d3, "M3: +geo+urban")
)
print(diagnostics)

# ======================
# 4) Prediction curves with 95% CIs
#     - vary one socio var over its observed range
#     - hold others (and physicals for M2/M3) at sample means
# ======================
socio_labels <- c(
  income_ln    = "Income (log)",
  age65_79     = "Age 65–79",
  age80_       = "Age 80+",
  percen_eu    = "EU foreigners",
  percen_noneu = "Non-EU foreigners",
  socisss      = "Social assistance",
  av_lvng      = "Living space per capita"
)

means_row <- dat_step %>%
  summarise(across(c(all_of(socio_vars), all_of(physical_vars)), ~mean(.x, na.rm = TRUE)))

make_grid <- function(var, n = 120) {
  rng  <- range(dat_step[[var]], na.rm = TRUE)
  base <- means_row[rep(1, n), , drop = FALSE]
  base[[var]] <- seq(rng[1], rng[2], length.out = n)
  base$var <- var; base$x <- base[[var]]
  base
}
pred_grid <- purrr::map_dfr(socio_vars, make_grid)

pred_long <- pred_grid %>%
  group_by(var) %>% tidyr::nest() %>%
  mutate(
    pred_m1 = purrr::map(data, ~ as_tibble(predict(m1, newdata = .x, se.fit = TRUE)) %>% mutate(x = .x$x)),
    pred_m2 = purrr::map(data, ~ as_tibble(predict(m2, newdata = .x, se.fit = TRUE)) %>% mutate(x = .x$x)),
    pred_m3 = purrr::map(data, ~ as_tibble(predict(m3, newdata = .x, se.fit = TRUE)) %>% mutate(x = .x$x))
  ) %>%
  select(-data) %>%
  pivot_longer(starts_with("pred_"), names_to = "Model", values_to = "pred") %>%
  mutate(Model = dplyr::recode(Model,
                               pred_m1 = "Model 1 (socio)",
                               pred_m2 = "Model 2 (+geo)",
                               pred_m3 = "Model 3 (+geo+urban)")) %>%
  unnest(pred) %>%
  rename(Predicted = fit, se = se.fit) %>%
  mutate(
    lwr = Predicted - 1.96 * se,
    upr = Predicted + 1.96 * se,
    var_lab = factor(var, levels = socio_vars, labels = socio_labels[socio_vars])
  ) %>%
  ungroup()

# Optional: vertical line at mean of each socio var
means_for_vline <- dat_step %>%
  summarise(across(all_of(socio_vars), ~mean(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "var", values_to = "xmean") %>%
  mutate(var_lab = factor(var, levels = socio_vars, labels = socio_labels[socio_vars]))

# ======================
# 5) Annotate each facet with β and p for M1/M2/M3
# ======================
get_beta_p <- function(model, model_label) {
  broom::tidy(model) %>%
    filter(term %in% socio_vars) %>%
    transmute(
      var      = term,
      var_lab  = factor(term, levels = socio_vars, labels = socio_labels[socio_vars]),
      Model    = model_label,
      beta     = estimate,
      p.value,
      label    = paste0("coef = ", round(beta, 3), ", ",
                        ifelse(is.na(p.value), "p = —",
                               ifelse(p.value < 0.001, "p < 0.001",
                                      paste0("p = ", round(p.value, 3)))))
    )
}
anno_df <- bind_rows(
  get_beta_p(m1, "Model 1 (socio)"),
  get_beta_p(m2, "Model 2 (+geo)"),
  get_beta_p(m3, "Model 3 (+geo+urban)")
)

# Position stacked labels in top-left (change CORNER to "tr" for top-right)
CORNER <- "tl"
order_lvls <- c("Model 1 (socio)", "Model 2 (+geo)", "Model 3 (+geo+urban)")
anno_df <- anno_df %>%
  mutate(
    row   = as.integer(factor(Model, levels = order_lvls)),
    x_inf = ifelse(CORNER == "tr", Inf, -Inf),
    y_inf = Inf,
    hjust = ifelse(CORNER == "tr", 1.02, -0.02),
    vjust = c(1.50, 3.00, 4.50)[row]
  )

# ======================
# 6) Plot
# ======================
p <- ggplot(pred_long, aes(x = x, y = Predicted, color = Model, fill = Model)) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.15, color = NA, show.legend = TRUE) +
  geom_line(linewidth = 0.8) +
  geom_vline(data = means_for_vline, aes(xintercept = xmean),
             linetype = "dashed", linewidth = 0.4, color = "grey50") +
  geom_text(data = anno_df,
            aes(x = x_inf, y = y_inf, label = label, color = Model, hjust = hjust, vjust = vjust),
            inherit.aes = FALSE, size = 3) +
  facet_wrap(~ var_lab, scales = "free_x", ncol = 4) +
  labs(y = "Composite heat exposure (CHEI, equal weights)", x = "Standardized socio-demographic variable") +
  theme_minimal() +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "bottom"
  ) +
  guides(color = guide_legend(override.aes = list(linetype = 1, linewidth = 1)))

print(p)

