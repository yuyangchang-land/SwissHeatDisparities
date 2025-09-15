#######################compare the 3 models by dividing the physical environmental vars and socio-demographic vars
# ================================
# Libraries
# ================================
library(dplyr)
library(readr)
library(car)        # vif()
library(lmtest)     # coeftest()
library(sandwich)   # vcovHC()
library(broom)      # tidy()
library(ggplot2)    # optional plot

# ================================
# 1) Read & derive variables
# ================================
df <- read_csv("/Users/yuchang/PhD/PhD_research/swiss_heat_exposure/data/revision_r1/muni_all_vars_panel_raster.csv")

# Derive: age65_79, EU/non‑EU shares (as discussed)
dat <- df %>%
  mutate(
    for_eu_prop = ifelse(for_eu > 1, for_eu/100, for_eu),  # ensure proportion
    age65_79    = age65_ - age80_,
    percen_eu    = foreign19 * for_eu_prop,
    percen_noneu = foreign19 * (1 - for_eu_prop)
  )

# Outcomes
outcomes <- c(
  LST = "lst_urban",
  HWD = "hetwv_d",
  HWP = "hetwv_p"
)

# Predictors (final sets you’re using)
physical_vars <- c("dem_mean","imp_mean","evi_mean","bldg_mean","albd_mean")
socio_vars    <- c("income","age65_79","age80_","percen_noneu","percen_eu",
                   "socioassis","ave_living")

# Keep needed cols; drop NAs
need <- c(unname(outcomes), physical_vars, socio_vars, "num_reside")
dat  <- dat %>% select(any_of(need)) %>% na.omit()

# ================================
# 2) Weights & scaling (match your old rules)
# ================================
# WEIGHTS: raw log(pop + 1), no normalization (your old manuscript style)
#dat <- dat %>% mutate(w = log(num_reside + 1))
dat <- dat %>% mutate(w = sqrt(num_reside))

# LOG INCOME before scaling; then z‑score predictors & outcomes (NOT weights)
std <- function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
dat <- dat %>%
  mutate(income = log(income + 1)) %>%  # <-- log first (matching your old approach)
  mutate(across(c(all_of(physical_vars), all_of(socio_vars), all_of(unname(outcomes))), std))

# ================================
# 3) Helper to run WLS + robust SE
# ================================
run_wls <- function(yvar, xvars, label, data = dat) {
  f   <- as.formula(paste(yvar, "~", paste(xvars, collapse = "+")))
  mod <- lm(f, data = data, weights = w)
  
  # robust (HC3) SE table
  vc  <- vcovHC(mod, type = "HC3")
  tidy_coef <- broom::tidy(coeftest(mod, vcov. = vc)) %>%
    rename(estimate = estimate, std.error = std.error,
           statistic = statistic, p.value = p.value) %>%
    mutate(model = label, outcome = yvar)
  
  # diagnostics
  res  <- residuals(mod)
  diag <- tibble(
    model   = label,
    outcome = yvar,
    n       = nobs(mod),
    R2      = summary(mod)$r.squared,
    Adj_R2  = summary(mod)$adj.r.squared,
    RMSE    = sqrt(mean(res^2)),
    MAE     = mean(abs(res)),
    AIC     = AIC(mod)
  )
  list(diag = diag, coefs = tidy_coef, fit = mod)
}

# ================================
# 4) Multicollinearity (VIF)
# ================================
vif_m0 <- vif(lm(as.formula(paste(outcomes[1], "~", paste(physical_vars, collapse="+"))), data=dat))
vif_m1 <- vif(lm(as.formula(paste(outcomes[1], "~", paste(socio_vars,    collapse="+"))), data=dat))
vif_m2 <- vif(lm(as.formula(paste(outcomes[1], "~", paste(c(physical_vars, socio_vars), collapse="+"))), data=dat))
print(list(VIF_M0 = vif_m0, VIF_M1 = vif_m1, VIF_M2 = vif_m2))

# ================================
# 5) Fit models M0, M1, M2, M3
# ================================
diag_list <- list(); coef_list <- list()

for (y in unname(outcomes)) {
  # M0 physical-only
  r0 <- run_wls(y, physical_vars, "M0_physical")
  diag_list[[length(diag_list)+1]] <- r0$diag
  coef_list[[length(coef_list)+1]] <- r0$coefs
  
  # M1 socio-only
  r1 <- run_wls(y, socio_vars, "M1_socio")
  diag_list[[length(diag_list)+1]] <- r1$diag
  coef_list[[length(coef_list)+1]] <- r1$coefs
  
  # M2 combined
  r2 <- run_wls(y, c(physical_vars, socio_vars), "M2_combined")
  diag_list[[length(diag_list)+1]] <- r2$diag
  coef_list[[length(coef_list)+1]] <- r2$coefs
  
  # M3 residualization (label by outcome so rows are clear)
  mod_phys <- lm(as.formula(paste(y, "~", paste(physical_vars, collapse="+"))), data=dat, weights=w)
  dat$resid_phys <- residuals(mod_phys)
  r3 <- run_wls("resid_phys", socio_vars, paste0("M3_residualization_", y))
  diag_list[[length(diag_list)+1]] <- r3$diag
  coef_list[[length(coef_list)+1]] <- r3$coefs
}

diagnostics  <- dplyr::bind_rows(diag_list) %>%
  dplyr::group_by(outcome) %>%
  dplyr::mutate(delta_AIC = AIC - min(AIC)) %>%
  dplyr::arrange(outcome, delta_AIC) %>%
  dplyr::ungroup()

coefficients <- dplyr::bind_rows(coef_list)

# ================================
# 6) Save outputs
# ================================
write_csv(diagnostics,  "wls_diagnostics.csv")
write_csv(coefficients, "wls_coefficients.csv")

print(diagnostics)

#######################mapping coefficients for model 2

# ================================
# 7) OPTIONAL: Plot M2 coefficients with 95% CI
# ================================
coefs <- coefficients

m2 <- coefs %>%
  filter(model == "M2_combined", term != "(Intercept)") %>%
  mutate(
    outcome = dplyr::recode(outcome,
                            "lst_urban" = "Land surface temperature",
                            "hetwv_d"   = "Heat warning days",
                            "hetwv_p"   = "Heatwave probability"),
    # set factor order here:
    outcome = factor(outcome,
                     levels = c("Land surface temperature",
                                "Heat warning days",
                                "Heatwave probability")),
    var_label = dplyr::recode(term,
                              "dem_mean"     = "Elevation",
                              "imp_mean"     = "Impervious %",
                              "evi_mean"     = "Vegetation (EVI)",
                              "bldg_mean"    = "Urban form / height",
                              "albd_mean"    = "Surface albedo",
                              "income"       = "Income (log)",
                              "age65_79"     = "Age 65–79",
                              "age80_"       = "Age 80+",
                              "percen_eu"    = "EU foreigners",
                              "percen_noneu" = "Non‑EU foreigners",
                              "socioassis"   = "Social assistance",
                              "ave_living"   = "Living space per capita",
                              .default       = term
    ),
    ci_low  = estimate - 1.96 * std.error,
    ci_high = estimate + 1.96 * std.error,
    sig_cat = dplyr::case_when(
      p.value < 0.05 & estimate > 0 ~ "Positive (p<0.05)",
      p.value < 0.05 & estimate < 0 ~ "Negative (p<0.05)",
      TRUE                           ~ "Not significant"
    )
  )

var_order <- c("Impervious %","Urban form / height","Surface albedo","Vegetation (EVI)","Elevation",
               "Income (log)","EU foreigners","Non‑EU foreigners","Social assistance","Living space per capita",
               "Age 65–79","Age 80+")
m2 <- m2 %>% mutate(var_label = factor(var_label, levels = var_order))

cols <- c("Positive (p<0.05)"="firebrick3",
          "Negative (p<0.05)"="#f4a300",
          "Not significant"="grey60")

p <- ggplot(m2, aes(x = var_label, y = estimate, fill = sig_cat)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2, color = "black") +
  facet_wrap(~ outcome, nrow = 1, scales = "free_x") +
  coord_flip() +
  scale_fill_manual(values = cols, name = NULL) +
  labs(x = NULL, y = "Standardized estimate (95% CI)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

# Save figure (optional)
ggsave("M2_coefficients.png", p, width = 10, height = 3.6, dpi = 300)
ggsave("M2_coefficients.pdf", p, width = 10, height = 3.6)

print(p)

# ================================
# Plot M3 (residualization) coefficients
# ================================
library(dplyr)
library(ggplot2)

coefs <- coefficients  # from your earlier code

m3 <- coefs %>%
  dplyr::filter(grepl("^M3_residualization_", model), term != "(Intercept)") %>%
  dplyr::mutate(
    # map back to original outcome names using the model string
    outcome = dplyr::case_when(
      grepl("lst_urban", model) ~ "Land surface temperature",
      grepl("hetwv_d",   model) ~ "Heat warning days",
      grepl("hetwv_p",   model) ~ "Heatwave probability",
      TRUE ~ "Residuals"
    ),
    # enforce facet order
    outcome = factor(outcome,
                     levels = c("Land surface temperature",
                                "Heat warning days",
                                "Heatwave probability")),
    # nice labels for socio variables (force dplyr::recode to avoid car::recode clash)
    var_label = dplyr::recode(
      term,
      "income"       = "Income (log)",
      "age65_79"     = "Age 65–79",
      "age80_"       = "Age 80+",
      "percen_eu"    = "EU foreigners",
      "percen_noneu" = "Non‑EU foreigners",
      "socioassis"   = "Social assistance",
      "ave_living"   = "Living space per capita",
      .default       = term
    ),
    ci_low  = estimate - 1.96 * std.error,
    ci_high = estimate + 1.96 * std.error,
    sig_cat = dplyr::case_when(
      p.value < 0.05 & estimate > 0 ~ "Positive (p<0.05)",
      p.value < 0.05 & estimate < 0 ~ "Negative (p<0.05)",
      TRUE                           ~ "Not significant"
    )
  )

# Order bars (socio vars only) for readability
var_order_m3 <- c("Income (log)",
                  "EU foreigners","Non‑EU foreigners",
                  "Social assistance","Living space per capita",
                  "Age 65–79","Age 80+")
m3 <- m3 %>% dplyr::mutate(var_label = factor(var_label, levels = var_order_m3))

cols <- c("Positive (p<0.05)"="firebrick3",
          "Negative (p<0.05)"="#f4a300",
          "Not significant"="grey60")

p_m3 <- ggplot(m3, aes(x = var_label, y = estimate, fill = sig_cat)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high), width = 0.2, color = "black") +
  facet_wrap(~ outcome, nrow = 1, scales = "free_x") +
  coord_flip() +
  scale_fill_manual(values = cols, name = NULL) +
  labs(x = NULL,
       y = "Standardized estimate (95% CI)",
       title = "Stage 2 (M3): Socio‑demographic associations after controlling for physical variables") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

print(p_m3)
