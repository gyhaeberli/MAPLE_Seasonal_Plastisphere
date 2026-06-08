################################################################################
# GLLVM PRESENCE-ABSENCE CO-OCCURRENCE PIPELINE
# FOULING (biofilm) + WATER communities in ONE script
#
# This script runs, for BOTH communities:
#   Section 3-A  Pollock decomposition (design vs residual correlation)
#   Pollock scatterplots
#   Quadrant summary
#   Top "CoOc pairs + surprise scatterplot
#   Circos chord diagrams (by association type, and by Pollock quadrant)
#
# HOW THE TWO COMMUNITIES DIFFER (everything else is identical):
#   - Data + models are stored in different files, loaded in different ways.
#   - Environmental formula:
#         fouling = site + season + substrate + all 2-way interactions
#         water   = site + season   (no substrate)
#   - Output directories.
# Because the analysis itself is the same, it is written ONCE as the function
# `run_pollock_circos()` and called once per community in a loop at the bottom.
################################################################################

# ── Libraries (only those needed up to the circos plots) ──────────────────────
library(gllvm)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(scales)
library(RColorBrewer)
library(ggrepel)
library(circlize)
library(flextable)  
library(officer)     # supporting package for flextable Word output

cat("=== COMBINED GLLVM CO-OCCURRENCE SCRIPT (fouling + water) ===\n\n")

setwd("~/Github/MAPLE_Seasonal_Plastisphere/Scripts/03.StatisticalAnalysis")


# ==============================================================================
# GLOBAL SETTINGS (shared by both communities)
# ==============================================================================

# Where the saved data + fitted GLLVM models live
SAVE_DIR_BASE <- "~/Github/MAPLE_Seasonal_Plastisphere/Processed_data/gllvm_models"

# |residual correlation| above this counts as a "strong" association.
# Used both to flag strong pairs and to build the circos diagrams.
R_STRONG <- 0.5

# Circos drawing constants (identical in both original scripts)
CANVAS_SIZE  <- 12      # inches
PLOT_DPI     <- 300
GAP_DEGREES  <- 7       # gap between sectors, in degrees
LABEL_CEX    <- 1.5     # sector label font size
LABEL_OFFSET <- 0.1    # how far labels sit outside the track


# ── patch_model() ─────────────────────────────────────────────────────────────
# The fitted gllvm objects were saved without a fully populated $call. Several
# downstream gllvm functions read fields from $call, so we copy the relevant
# pieces back in. include_X = FALSE for the null model (it has no predictors).
patch_model <- function(fit, num_lv, include_X = TRUE) {
  fit$call$y           <- fit$y
  fit$call$num.lv      <- num_lv
  fit$call$studyDesign <- fit$studyDesign
  if (!is.null(fit$offset))          fit$call$offset <- fit$offset
  if (include_X && !is.null(fit$X))  fit$call$X      <- fit$X
  fit
}


# ── draw_label_track() ────────────────────────────────────────────────────────
# Draws the outer class-name labels on a circos plot.
draw_label_track <- function() {
  circos.trackPlotRegion(
    track.index  = 1,
    track.height = 0.08,
    panel.fun = function(x, y) {
      xlim   <- get.cell.meta.data("xlim")
      ylim   <- get.cell.meta.data("ylim")
      sector <- get.cell.meta.data("sector.index")
      circos.text(
        mean(xlim), ylim[2] + LABEL_OFFSET, sector,
        facing = "clockwise", niceFacing = TRUE,
        adj = c(0, 0.5), cex = LABEL_CEX, font = 2
      )
    },
    bg.border = NA
  )
}


# ── track_dirs() / ensure_dirs() ──────────────────────────────────────────────
# Resolve and (optionally) create the output folders for a community.
track_dirs <- function(cfg) {
  out <- cfg$out_dir
  list(out   = out,
       table = file.path(out, "Tables"),
       fig   = file.path(out, "Figures"))
}
ensure_dirs <- function(cfg) {
  d <- track_dirs(cfg)
  for (p in c(d$out, d$table, d$fig))
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
  invisible(d)
}


# ==============================================================================
# DATA + MODEL LOADERS (the part that genuinely differs between communities)
# Each returns: list(model_data = ..., best_model = ..., null_model = ...)
# ==============================================================================

# --- Fouling / biofilm: single saved list + separate model files; num.lv = 2 --
load_fouling <- function() {
  model_data <- readRDS(file.path(SAVE_DIR_BASE, "data_for_gllvm_FINAL.rds"))
  
  best_model <- patch_model(
    readRDS(file.path(SAVE_DIR_BASE, "GLLVM_final_incM2.rds")),
    num_lv = 2, include_X = TRUE
  )
  null_model <- patch_model(
    readRDS(file.path(SAVE_DIR_BASE, "GLLVM_null_incM7.rds")),
    num_lv = 2, include_X = FALSE
  )
  list(model_data = model_data, best_model = best_model, null_model = null_model)
}

# --- Water: one bundled environment file; num.lv from params$NUM_LV -----------
# NOTE: assumes env_water_WF_for_HPC.rds has top-level metadata, tax_mat,
# otu_id_to_label, top_otus, otu_filt, and a nested list `params` with NUM_LV.
load_water <- function() {
  env_list <- readRDS(file.path(SAVE_DIR_BASE, "env_water_WF_for_HPC.rds"))
  
  num_lv <- env_list$params$NUM_LV
  
  model_data <- list(
    metadata        = env_list$metadata,
    tax_mat         = env_list$tax_mat,
    otu_id_to_label = env_list$otu_id_to_label,
    top_otus        = env_list$top_otus,
    otu_filt        = env_list$otu_filt
  )
  
  best_model <- patch_model(
    readRDS(file.path(SAVE_DIR_BASE, "GLLVM_final_water_WF_add.rds")),
    num_lv = num_lv, include_X = TRUE
  )
  null_model <- patch_model(
    readRDS(file.path(SAVE_DIR_BASE, "GLLVM_null_water_WF.rds")),
    num_lv = num_lv, include_X = FALSE
  )
  list(model_data = model_data, best_model = best_model, null_model = null_model)
}


# ==============================================================================
# PER-COMMUNITY CONFIGURATION
# ==============================================================================

tracks <- list(
  
  fouling = list(
    label            = "Fouling (biofilm)",
    loader           = load_fouling,
    env_formula      = ~ site + season + substrate +
      site:season + site:substrate + season:substrate,
    trace_predictors = "site + season + substrate",
    out_dir          = "~/Github/MAPLE_Seasonal_Plastisphere/Results/gllvm_results",
    file_tag         = "biof",
    n_surprise       = 100,   # pairs highlighted/labelled on the surprise plot + CSV
    n_circos         = 50     # pairs aggregated to class level for the circos plots
  ),
  
  water = list(
    label            = "Water",
    loader           = load_water,
    env_formula      = ~ site + season,
    trace_predictors = "site + season",
    out_dir          = "~/Github/MAPLE_Seasonal_Plastisphere/Results/gllvm_results/water",
    file_tag         = "water",
    n_surprise       = 100,
    n_circos         = 50
  )
)




################################################################################
#   build_pollock_pairs()    -> Pollock pairs table + CSV   (the shared input)
#   plot_pollock_quadrant()  -> Pollock scatterplot, strong pairs by quadrant
#   plot_pollock_sameclass() -> Pollock scatterplot, same-class pairs in blue
#   summarise_quadrants()    -> Quadrant summary table (CSV + Word)
#   top_surprising_pairs()   -> Top "surprising" pairs table + CSV
#   plot_pollock_surprise()  -> Pollock surprise scatterplot
#   circos_by_association()  -> Circos chord diagram, positive vs negative
#   circos_by_quadrant()     -> Circos chord diagram, coloured by quadrant
#
# DATA FLOW:
#   build_pollock_pairs() produces `pair_df`. EVERY other step takes
#   (cfg, pair_df) and is self-contained, so you can run any single step on its
#   own. The two circos plots also need a class-level aggregation, built once by
#   the helper prepare_circos_links() and shared between them.
#
################################################################################

# ── SETUP (run once) ──────────────────────────────────────────────────────────
cfg_f <- tracks$fouling;  ensure_dirs(cfg_f);  loaded_f <- cfg_f$loader()
cfg_w <- tracks$water;    ensure_dirs(cfg_w);  loaded_w <- cfg_w$loader()

# ==============================================================================
# STEP 1: build_pollock_pairs()
# Computes design + residual correlations, the trace reduction, and the
# OTU-pair table; writes the Pollock pairs CSV. Returns `pair_df`, the shared
# input for every other step.
# ==============================================================================

build_pollock_pairs <- function(cfg, loaded) {
  
  cat("--- [1] Pollock pairs (", cfg$label, ") ---\n\n", sep = "")
  d <- track_dirs(cfg)
  
  model_data <- loaded$model_data
  best_model <- loaded$best_model
  null_model <- loaded$null_model
  
  # ── Environmental correlation ───────────────────────────────────────────────
  B     <- coef(best_model)$Xcoef
  X_raw <- best_model$X
  
  # Same dummy-coded design matrix the model used internally ([, -1] drops intercept)
  X_numeric <- model.matrix(cfg$env_formula, data = X_raw)[, -1]
  
  if (!identical(colnames(X_numeric), colnames(B))) {
    cat("Mismatch — only in X_numeric:", setdiff(colnames(X_numeric), colnames(B)), "\n")
    cat("Only in B:",                    setdiff(colnames(B), colnames(X_numeric)), "\n")
    stop("Design matrix columns do not match coefficient matrix. Check factor levels.")
  }
  
  fitted_env <- X_numeric %*% t(B)
  env_cor    <- cor(fitted_env)
  diag(env_cor) <- NA
  
  # ── Residual correlation ─────────────────────────────────────────────────────
  res_cor <- getResidualCor(best_model)
  diag(res_cor) <- NA
  
  # ── Trace reduction (how much co-occurrence the predictors explain) ──────────
  rcov_best <- getResidualCov(best_model, adjust = 0)
  rcov_null <- getResidualCov(null_model, adjust = 0)
  trace_reduction_pct <- round((1 - rcov_best$trace / rcov_null$trace) * 100, 1)
  
  cat("Residual covariance trace — null model:", round(rcov_null$trace, 3), "\n")
  cat("Residual covariance trace — best model:", round(rcov_best$trace, 3), "\n")
  cat("Reduction by", cfg$trace_predictors, ":", trace_reduction_pct, "%\n\n")
  
  # ── Taxonomy lookup ──────────────────────────────────────────────────────────
  otu_labels  <- colnames(res_cor)
  label_to_id <- setNames(names(model_data$otu_id_to_label),
                          model_data$otu_id_to_label)
  otu_ids <- label_to_id[otu_labels]
  
  tax_cols_use <- intersect(
    c("Division", "Class", "Order", "Family", "Genus", "Species"),
    colnames(model_data$tax_mat)
  )
  
  tax_lookup <- model_data$tax_mat[otu_ids, tax_cols_use, drop = FALSE] %>%
    as.data.frame() %>%
    rownames_to_column("OTU_ID") %>%
    mutate(OTU_Label = model_data$otu_id_to_label[OTU_ID])
  
  # ── OTU-pair data frame (upper triangle only) ────────────────────────────────
  pairs <- which(upper.tri(res_cor), arr.ind = TRUE)
  
  pair_df <- data.frame(
    OTU_j   = otu_labels[pairs[, 1]],
    OTU_k   = otu_labels[pairs[, 2]],
    res_cor = res_cor[pairs],
    env_cor = env_cor[pairs],
    stringsAsFactors = FALSE
  ) %>%
    left_join(tax_lookup %>% select(OTU_Label, Class_j = Class, Division_j = Division),
              by = c("OTU_j" = "OTU_Label")) %>%
    left_join(tax_lookup %>% select(OTU_Label, Class_k = Class, Division_k = Division),
              by = c("OTU_k" = "OTU_Label")) %>%
    mutate(
      same_class    = !is.na(Class_j) & !is.na(Class_k) & Class_j == Class_k,
      same_division = !is.na(Division_j) & !is.na(Division_k) & Division_j == Division_k,
      # Quadrant: Q1(++) shared niche, Q2(+-) competition, Q3(-+) facilitation,
      #           Q4(--) niche differentiation.
      quadrant = case_when(
        env_cor >= 0 & res_cor >= 0 ~ "Q1",
        env_cor >= 0 & res_cor <  0 ~ "Q2",
        env_cor <  0 & res_cor >= 0 ~ "Q3",
        env_cor <  0 & res_cor <  0 ~ "Q4"
      ),
      surprise   = abs(res_cor) - abs(env_cor),   # co-occurrence score
      strong_res = abs(res_cor) > R_STRONG,
      concordant = sign(res_cor) == sign(env_cor)
    ) %>%
    filter(!is.na(res_cor), !is.na(env_cor))
  
  cat("Total OTU pairs:", nrow(pair_df), "\n")
  cat("Strong residual pairs (|r| >", R_STRONG, "):",
      sum(pair_df$strong_res),
      sprintf("(%.1f%%)\n", sum(pair_df$strong_res) / nrow(pair_df) * 100))
  cat("Quadrant breakdown:\n")
  print(table(pair_df$quadrant))
  cat("\n")
  
  write.csv(pair_df,
            file.path(d$table, paste0("S3A_pollock_pairs_", cfg$file_tag, ".csv")),
            row.names = FALSE)
  
  # Carry the trace reduction along as an attribute (handy for titles/reporting)
  attr(pair_df, "trace_reduction_pct") <- trace_reduction_pct
  pair_df
}




pair_df_f <- build_pollock_pairs(cfg_f, loaded_f)
pair_df_w <- build_pollock_pairs(cfg_w, loaded_w)

# ==================================loaded = # =================================
# STEP 2a: plot_pollock_quadrant()   — strong pairs coloured by quadrant
# ==============================================================================

plot_pollock_quadrant <- function(cfg, pair_df) {
  
  cat("--- [2a] Pollock scatterplot: quadrant (", cfg$label, ") ---\n", sep = "")
  d <- track_dirs(cfg)
  
  quad_cols_pt <- c("Q1" = "#2166AC", "Q2" = "#F4A582",
                    "Q3" = "#2AB5A0", "Q4" = "#B2182B")
  
  p <- ggplot(pair_df, aes(x = env_cor, y = res_cor)) +
    geom_hline(yintercept = 0, linewidth = 0.5, colour = "grey40") +
    geom_vline(xintercept = 0, linewidth = 0.5, colour = "grey40") +
    geom_point(data = filter(pair_df, !strong_res),
               colour = "grey85", alpha = 0.3, size = 1.2) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "grey60", linewidth = 0.4) +
    geom_point(data = filter(pair_df, strong_res),
               aes(colour = quadrant), alpha = 0.75, size = 1.2) +
    annotate("text", x =  0.65, y =  0.9, label = "Q1", size = 2.7, colour = "grey25") +
    annotate("text", x =  0.65, y = -0.9, label = "Q2", size = 2.7, colour = "grey25") +
    annotate("text", x = -0.65, y =  0.9, label = "Q3", size = 2.7, colour = "grey25") +
    annotate("text", x = -0.65, y = -0.9, label = "Q4", size = 2.7, colour = "grey25") +
    scale_colour_manual(name = "Quadrant", values = quad_cols_pt) +
    coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1)) +
    labs(title = paste0("Pollock decomposition — ", cfg$label),
         x = "Design-based correlation\n",
         y = "Residual correlation\n") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), legend.position = "none")
  
  print(p)
  ggsave(file.path(d$fig, paste0("pollock_all_highlighted_", cfg$file_tag, ".png")),
         p, width = 8, height = 8, dpi = 200)
  invisible(p)
}



plot_pollock_quadrant(cfg_f, pair_df_f)
plot_pollock_quadrant(cfg_w, pair_df_w)




# ==============================================================================
# STEP 2b: plot_pollock_sameclass()  — same-class pairs in blue
# ==============================================================================

plot_pollock_sameclass <- function(cfg, pair_df) {
  
  cat("--- [2b] Pollock scatterplot: same class (", cfg$label, ") ---\n", sep = "")
  d <- track_dirs(cfg)
  
  p <- ggplot(pair_df, aes(x = env_cor, y = res_cor)) +
    geom_hline(yintercept = 0, linewidth = 0.5, colour = "grey40") +
    geom_vline(xintercept = 0, linewidth = 0.5, colour = "grey40") +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "grey60", linewidth = 0.4) +
    geom_point(data = filter(pair_df, !same_class),
               colour = "grey70", alpha = 0.25, size = 0.7) +
    geom_point(data = filter(pair_df, same_class),
               colour = "blue", alpha = 0.55, size = 0.8) +
    annotate("text", x =  0.65, y =  0.9, label = "Q1", size = 2.7, colour = "grey25") +
    annotate("text", x =  0.65, y = -0.9, label = "Q2", size = 2.7, colour = "grey25") +
    annotate("text", x = -0.65, y =  0.9, label = "Q3", size = 2.7, colour = "grey25") +
    annotate("text", x = -0.65, y = -0.9, label = "Q4", size = 2.7, colour = "grey25") +
    coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1)) +
    labs(
      # title    = paste0("Pollock decomposition — ", cfg$label),
      # subtitle = paste0("Each point = one OTU pair (n = ", nrow(pair_df), ")\n",
      #                   "Blue = same taxonomic class | Dashed diagonal = env_cor == res_cor"),
      x = paste0("Design-based correlation\n(shared response to ", cfg$trace_predictors, ")"),
      y = "Residual correlation\n(co-occurrence beyond fixed effects)"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))
  
  print(p)
  ggsave(file.path(d$fig, paste0("pollock_sameclass_", cfg$file_tag, ".png")),
         p, width = 8, height = 8, dpi = 200)
  invisible(p)
}




plot_pollock_sameclass(cfg_f, pair_df_f)
plot_pollock_sameclass(cfg_w, pair_df_w)


# ==============================================================================
# STEP 3: summarise_quadrants()  — quadrant summary table (CSV + Word)
# ==============================================================================

summarise_quadrants <- function(cfg, pair_df) {
  
  cat("--- [3] Quadrant summary (", cfg$label, ") ---\n", sep = "")
  d <- track_dirs(cfg)
  
  quadrant_summary <- pair_df %>%
    filter(!is.na(quadrant)) %>%
    group_by(quadrant) %>%
    summarise(
      n_pairs        = n(),
      pct_pairs      = round(n() / nrow(pair_df) * 100, 1),
      n_strong       = sum(strong_res),
      pct_strong     = round(sum(strong_res) / sum(pair_df$strong_res) * 100, 1),
      mean_res_cor   = round(mean(res_cor),  3),
      mean_env_cor   = round(mean(env_cor),  3),
      mean_surprise  = round(mean(surprise), 3),
      n_same_class   = sum(same_class),
      pct_same_class = round(sum(same_class) / sum(pair_df$same_class) * 100, 1),
      .groups        = "drop"
    )
  
  print(quadrant_summary, n = Inf)
  cat("\n")
  
  write.csv(quadrant_summary,
            file.path(d$table, paste0("S3A_quadrant_summary_", cfg$file_tag, ".csv")),
            row.names = FALSE)
  
  quad_sum <- flextable(quadrant_summary) %>%
    set_header_labels(
      quadrant = "Quadrant", n_pairs = "N pairs", pct_pairs = "% pairs",
      n_strong = "N strong", pct_strong = "% strong",
      mean_res_cor = "Mean res. cor.", mean_env_cor = "Mean env. cor.",
      mean_surprise = "Mean CoOc Score",
      n_same_class = "N same class", pct_same_class = "% same class"
    ) %>%
    theme_vanilla() %>%
    autofit()
  save_as_docx(quad_sum,
               path = file.path(d$table, paste0("quadrant_summary_", cfg$file_tag, ".docx")))
  
  invisible(quadrant_summary)
}




quad_summary_f <- summarise_quadrants(cfg_f, pair_df_f)
quad_summary_w <- summarise_quadrants(cfg_w, pair_df_w)



# ==============================================================================
# STEP 4: top_surprising_pairs()  — top "surprising" pairs table + CSV
# ==============================================================================

top_surprising_pairs <- function(cfg, pair_df, n = cfg$n_surprise) {
  
  cat("--- [4] Top surprising pairs (", cfg$label, ") ---\n", sep = "")
  d <- track_dirs(cfg)
  
  top_surprise <- pair_df %>%
    filter(strong_res) %>%
    arrange(desc(surprise)) %>%
    slice_head(n = n) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))
  
  cat("Top", n, "most surprising pairs (high residual, low environmental):\n")
  print(top_surprise %>% select(OTU_j, OTU_k, Class_j, Class_k,
                                env_cor, res_cor, surprise, quadrant))
  cat("\n")
  
  write.csv(top_surprise,
            file.path(d$table, paste0("S3A_top_surprising_pairs_", cfg$file_tag, ".csv")),
            row.names = FALSE)
  
  invisible(top_surprise)
}



top_surprise_f <- top_surprising_pairs(cfg_f, pair_df_f)
top_surprise_w <- top_surprising_pairs(cfg_w, pair_df_w)



# ==============================================================================
# STEP 5: plot_pollock_surprise()  — Pollock surprise scatterplot
# Recomputes the top-n subset internally so it can run standalone. Lower `n`
# (e.g. 30) if the labels are too crowded.
# ==============================================================================

plot_pollock_surprise <- function(cfg, pair_df, n = cfg$n_surprise) {
  
  cat("--- [5] Pollock surprise scatterplot (", cfg$label, ") ---\n", sep = "")
  d <- track_dirs(cfg)
  
  top_surprise <- pair_df %>%
    filter(strong_res) %>%
    arrange(desc(surprise)) %>%
    slice_head(n = n) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))
  
  p <- ggplot() +
    geom_point(data = pair_df, aes(x = env_cor, y = res_cor),
               colour = "grey80", alpha = 0.25, size = 0.6) +
    geom_point(data = filter(pair_df, strong_res), aes(x = env_cor, y = res_cor),
               colour = "grey50", alpha = 0.4, size = 0.9) +
    geom_point(data = top_surprise, aes(x = env_cor, y = res_cor),
               colour = "orange", size = 3, alpha = 0.9) +
    ggrepel::geom_text_repel(
      data = top_surprise,
      aes(x = env_cor, y = res_cor,
          label = paste0(sub("_ASV_.*", "", OTU_j), " × ", sub("_ASV_.*", "", OTU_k))),
      colour = "darkred", size = 4, max.overlaps = 50,
      segment.size = 0.3, segment.alpha = 0.5, box.padding = 0.4
    ) +
    geom_hline(yintercept = 0, linewidth = 0.5, colour = "grey30") +
    geom_vline(xintercept = 0, linewidth = 0.5, colour = "grey30") +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "grey50", linewidth = 0.4) +
    geom_hline(yintercept =  R_STRONG, linetype = "dotted", colour = "grey40", linewidth = 0.4) +
    geom_hline(yintercept = -R_STRONG, linetype = "dotted", colour = "grey40", linewidth = 0.4) +
    annotate("text", x =  0.65, y =  0.95, label = "Q1", size = 2.8, colour = "grey30") +
    annotate("text", x =  0.65, y = -0.95, label = "Q2", size = 2.8, colour = "grey30") +
    annotate("text", x = -0.65, y =  0.95, label = "Q3", size = 2.8, colour = "grey30") +
    annotate("text", x = -0.65, y = -0.95, label = "Q4", size = 2.8, colour = "grey30") +
    coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1)) +
    labs(
      title = paste0("Pollock decomposition — ", cfg$label, " (top ", n, " surprising pairs)"),
      x = paste0("Design-based correlation\n(shared response to ", cfg$trace_predictors, ")"),
      y = "Residual correlation\n(co-occurrence beyond fixed effects)"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 8))
  
  print(p)
  ggsave(file.path(d$fig, paste0("pollock_surprise_", cfg$file_tag, ".png")),
         p, width = 10, height = 10, dpi = 200)
  invisible(p)
}


plot_pollock_surprise(cfg_f, pair_df_f)
plot_pollock_surprise(cfg_w, pair_df_w)


# ==============================================================================
# HELPER: prepare_circos_links()
# Aggregates the top-n_circos surprising (classified) pairs to the class level
# and builds the colour palettes. Used by BOTH circos functions. Returns a list
# with everything the chord diagrams need.
# ==============================================================================

prepare_circos_links <- function(cfg, pair_df) {
  
  circos_pairs <- pair_df %>%
    filter(strong_res,
           !is.na(Class_j), Class_j != "",
           !is.na(Class_k), Class_k != "") %>%
    arrange(desc(surprise)) %>%
    slice_head(n = cfg$n_circos) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)),
           association = if_else(res_cor > 0, "positive", "negative"))
  
  class_links <- circos_pairs %>%
    group_by(Class_j, Class_k, association) %>%
    summarise(
      mean_res_cor      = round(mean(abs(res_cor)), 3),
      mean_surprise     = round(mean(surprise), 3),
      n_pairs           = n(),
      dominant_quadrant = names(sort(table(quadrant), decreasing = TRUE))[1],
      .groups = "drop"
    ) %>%
    mutate(link_weight = mean_res_cor * n_pairs)   # strength x frequency
  
  all_classes <- sort(unique(c(class_links$Class_j, class_links$Class_k)))
  n_classes   <- length(all_classes)
  
  if (n_classes <= 8) {
    sector_cols <- setNames(brewer.pal(max(n_classes, 3), "Set2")[seq_len(n_classes)],
                            all_classes)
  } else if (n_classes <= 12) {
    sector_cols <- setNames(brewer.pal(n_classes, "Paired"), all_classes)
  } else {
    sector_cols <- setNames(colorRampPalette(brewer.pal(12, "Paired"))(n_classes),
                            all_classes)
  }
  
  list(
    class_links       = class_links,
    chord_df          = class_links %>%
      select(from = Class_j, to = Class_k, value = link_weight) %>%
      as.data.frame(),
    all_classes       = all_classes,
    sector_cols       = sector_cols,
    link_col_positive = adjustcolor("darkblue", alpha.f = 0.6),
    link_col_negative = adjustcolor("red2",     alpha.f = 0.7),
    quadrant_cols     = c(
      "Q1" = adjustcolor("#2166AC", alpha.f = 0.65),
      "Q2" = adjustcolor("#F4A582", alpha.f = 0.65),
      "Q3" = adjustcolor("#2AB5A0", alpha.f = 0.65),
      "Q4" = adjustcolor("#B2182B", alpha.f = 0.65)
    )
  )
}


links_f <- prepare_circos_links(cfg_f, pair_df_f)
links_w <- prepare_circos_links(cfg_w, pair_df_w)


# ==============================================================================
# STEP 6: circos_by_association()  — links coloured positive vs negative
# `links` defaults to a fresh prepare_circos_links() so the function is
# standalone; the orchestrator passes a pre-built one to avoid recomputing.
# ==============================================================================

circos_by_association <- function(cfg, pair_df, links = prepare_circos_links(cfg, pair_df)) {
  
  cat("--- [6] Circos: by association (", cfg$label, ") ---\n", sep = "")
  d <- track_dirs(cfg)
  
  # Give more canvas room for water (longer/more class names)
  canvas_lim <- if (cfg$file_tag == "water") 1.9 else 1.6
  
  png(file.path(d$fig, paste0("circos_by_association_", cfg$file_tag, ".png")),
      width = CANVAS_SIZE, height = CANVAS_SIZE, units = "in", res = PLOT_DPI)
  par(mar = c(10, 2, 4, 2))
  circos.clear()
  circos.par(gap.after = GAP_DEGREES, start.degree = 90, clock.wise = TRUE,
             track.margin = c(0.01, 0.05),
             canvas.xlim = c(-canvas_lim, canvas_lim), 
             canvas.ylim = c(-canvas_lim, canvas_lim),
             points.overflow.warning = FALSE)
  chordDiagram(
    x                     = links$chord_df,
    grid.col              = links$sector_cols,
    col                   = ifelse(links$class_links$association == "positive",
                                   links$link_col_positive, links$link_col_negative),
    transparency          = 0, directional = 0,
    annotationTrack       = "grid", annotationTrackHeight = 0.04,
    link.border           = NA, order = links$all_classes,
    link.sort             = TRUE, link.decreasing = TRUE,
    preAllocateTracks     = list(track.height = 0.08)
  )
  draw_label_track()
  par(xpd = TRUE)
  circos.clear()
  dev.off()
  cat("  saved.\n")
  invisible(links$class_links)
}


circos_by_association(cfg_f, pair_df_f, links_f)
circos_by_association(cfg_w, pair_df_w, links_w)



# ==============================================================================
# STEP 6.5: circos_by_assocaition_split()  — balanced positive and negative links
# ==============================================================================

prepare_circos_split <- function(cfg, pair_df, n_per_sign = 25) {
  
  top_positive <- pair_df %>%
    filter(strong_res,
           !is.na(Class_j), Class_j != "",
           !is.na(Class_k), Class_k != "",
           res_cor > 0) %>%
    arrange(desc(res_cor)) %>%
    slice_head(n = n_per_sign) %>%
    mutate(association = "positive")
  
  top_negative <- pair_df %>%
    filter(strong_res,
           !is.na(Class_j), Class_j != "",
           !is.na(Class_k), Class_k != "",
           res_cor < 0) %>%
    arrange(res_cor) %>%
    slice_head(n = n_per_sign) %>%
    mutate(association = "negative")
  
  circos_pairs <- bind_rows(top_positive, top_negative) %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))
  
  cat("  Circos subset:", nrow(top_positive), "positive +",
      nrow(top_negative), "negative pairs\n")
  
  class_links <- circos_pairs %>%
    group_by(Class_j, Class_k, association) %>%
    summarise(
      mean_res_cor      = round(mean(abs(res_cor)), 3),
      mean_surprise     = round(mean(surprise), 3),
      n_pairs           = n(),
      dominant_quadrant = names(sort(table(quadrant), decreasing = TRUE))[1],
      .groups = "drop"
    ) %>%
    mutate(link_weight = mean_res_cor * n_pairs)
  
  all_classes <- sort(unique(c(class_links$Class_j, class_links$Class_k)))
  n_classes   <- length(all_classes)
  
  if (n_classes <= 8) {
    sector_cols <- setNames(brewer.pal(max(n_classes, 3), "Set2")[seq_len(n_classes)],
                            all_classes)
  } else if (n_classes <= 12) {
    sector_cols <- setNames(brewer.pal(n_classes, "Paired"), all_classes)
  } else {
    sector_cols <- setNames(colorRampPalette(brewer.pal(12, "Paired"))(n_classes),
                            all_classes)
  }
  
  list(
    class_links       = class_links,
    chord_df          = class_links %>%
      select(from = Class_j, to = Class_k, value = link_weight) %>%
      as.data.frame(),
    all_classes       = all_classes,
    sector_cols       = sector_cols,
    link_col_positive = adjustcolor("darkblue", alpha.f = 0.6),
    link_col_negative = adjustcolor("red2",     alpha.f = 0.7),
    quadrant_cols     = c(
      "Q1" = adjustcolor("#2166AC", alpha.f = 0.65),
      "Q2" = adjustcolor("#F4A582", alpha.f = 0.65),
      "Q3" = adjustcolor("#2AB5A0", alpha.f = 0.65),
      "Q4" = adjustcolor("#B2182B", alpha.f = 0.65)
    )
  )
}

# Rebuild both link objects, then rerun the split circos
links_f <- prepare_circos_split(cfg_f, pair_df_f, n_per_sign = 50)
links_w <- prepare_circos_split(cfg_w, pair_df_w, n_per_sign = 50)





circos_by_association_split <- function(cfg, pair_df, links = prepare_circos_links(cfg, pair_df)) {
  
  cat("--- [6.5] Circos: positive vs negative split (", cfg$label, ") ---\n", sep = "")
  d <- track_dirs(cfg)
  
  canvas_lim <- if (cfg$file_tag == "water") 1.9 else 1.6
  
  # Separate link tables for each sign
  pos_links <- links$class_links %>% filter(association == "positive")
  neg_links <- links$class_links %>% filter(association == "negative")
  
  # Class sets per sign (sectors differ between panels)
  pos_classes <- sort(unique(c(pos_links$Class_j, pos_links$Class_k)))
  neg_classes <- sort(unique(c(neg_links$Class_j, neg_links$Class_k)))
  
  make_chord_df <- function(cl) {
    cl %>% select(from = Class_j, to = Class_k, value = link_weight) %>% as.data.frame()
  }
  
  png(file.path(d$fig, paste0("circos_split_posneg_", cfg$file_tag, ".png")),
      width  = CANVAS_SIZE * 2,   # two panels side by side
      height = CANVAS_SIZE,
      units  = "in", res = PLOT_DPI)
  
  layout(matrix(1:2, nrow = 1))
  
  # ── LEFT: positive associations ──────────────────────────────────────────────
  par(mar = c(6, 2, 6, 2))
  circos.clear()
  circos.par(gap.after = GAP_DEGREES, start.degree = 90, clock.wise = TRUE,
             track.margin = c(0.01, 0.05),
             canvas.xlim = c(-canvas_lim, canvas_lim),
             canvas.ylim = c(-canvas_lim, canvas_lim),
             points.overflow.warning = FALSE)
  
  chordDiagram(
    x                   = make_chord_df(pos_links),
    grid.col            = links$sector_cols[pos_classes],
    col                 = links$link_col_positive,
    transparency        = 0, directional = 0,
    annotationTrack     = "grid", annotationTrackHeight = 0.04,
    link.border         = NA, order = pos_classes,
    link.sort           = TRUE, link.decreasing = TRUE,
    preAllocateTracks   = list(track.height = 0.08)
  )
  draw_label_track()
  title("Positive associations", cex.main = 1.6, font.main = 2, line = 3)
  circos.clear()
  
  # ── RIGHT: negative associations ─────────────────────────────────────────────
  par(mar = c(6, 2, 6, 2))
  circos.clear()
  circos.par(gap.after = GAP_DEGREES, start.degree = 90, clock.wise = TRUE,
             track.margin = c(0.01, 0.05),
             canvas.xlim = c(-canvas_lim, canvas_lim),
             canvas.ylim = c(-canvas_lim, canvas_lim),
             points.overflow.warning = FALSE)
  
  chordDiagram(
    x                   = make_chord_df(neg_links),
    grid.col            = links$sector_cols[neg_classes],
    col                 = links$link_col_negative,
    transparency        = 0, directional = 0,
    annotationTrack     = "grid", annotationTrackHeight = 0.04,
    link.border         = NA, order = neg_classes,
    link.sort           = TRUE, link.decreasing = TRUE,
    preAllocateTracks   = list(track.height = 0.08)
  )
  draw_label_track()
  title("Negative associations", cex.main = 1.6, font.main = 2, line = 3)
  circos.clear()
  
  dev.off()
  cat("  saved.\n")
  invisible(links$class_links)
}

# Calls (after links_f / links_w are built):
circos_by_association_split(cfg_f, pair_df_f, links_f)
circos_by_association_split(cfg_w, pair_df_w, links_w)

# ==============================================================================
# STEP 7: circos_by_quadrant()  — links coloured by Pollock quadrant
# ==============================================================================

circos_by_quadrant <- function(cfg, pair_df, links = prepare_circos_links(cfg, pair_df)) {
  
  cat("--- [7] Circos: by quadrant (", cfg$label, ") ---\n", sep = "")
  d <- track_dirs(cfg)
  
  png(file.path(d$fig, paste0("circos_by_quadrant_", cfg$file_tag, ".png")),
      width = CANVAS_SIZE, height = CANVAS_SIZE, units = "in", res = PLOT_DPI)
  par(mar = c(10, 2, 4, 2))
  circos.clear()
  circos.par(gap.after = GAP_DEGREES, start.degree = 90, clock.wise = TRUE,
             track.margin = c(0.01, 0.05),
             canvas.xlim = c(-1.7, 1.7), canvas.ylim = c(-1.7, 1.7),
             points.overflow.warning = FALSE)
  chordDiagram(
    x                     = links$chord_df,
    grid.col              = links$sector_cols,
    col                   = links$quadrant_cols[links$class_links$dominant_quadrant],
    transparency          = 0, directional = 0,
    annotationTrack       = "grid", annotationTrackHeight = 0.04,
    link.border           = NA, order = links$all_classes,
    link.sort             = TRUE, link.decreasing = TRUE,
    preAllocateTracks     = list(track.height = 0.08)
  )
  draw_label_track()
  par(xpd = TRUE)
  legend(x = -1.35, y = -1.55,
         legend = c("Q1",
                    "Q2",
                    "Q3",
                    "Q4"),
         fill = unname(links$quadrant_cols), border = NA,
         title = "Pollock quadrant (link colour)", bty = "n", cex = 0.82, title.cex = 0.9)
  circos.clear()
  dev.off()
  cat("  saved.\n")
  invisible(links$class_links)
}



circos_by_quadrant(cfg_f, pair_df_f, links_f)
circos_by_quadrant(cfg_w, pair_df_w, links_w)


# ==============================================================================
# STEP 8: plot_synd_diatom_triangle()
# Triangular co-occurrence heatmap over the COMBINED set of Syndiniales +
# diatom (Bacillariophyceae + Mediophyceae) ASVs on a single shared axis.
#   Upper triangle = every pair among the set (fully tiled, no gaps)
#   fill = Pollock quadrant ; shade (alpha) = |residual r|
#   diagonal = group key (purple = Syndiniales, green = diatom)
#   axis labels (lowest taxonomy + ASV number) shown if the set is small enough
# ==============================================================================

plot_synd_diatom_triangle <- function(cfg, pair_df, loaded, links,
                                      lab_size = 9,
                                      show_labels = NA,
                                      strong_only = FALSE) {
  
  cat("--- [8] Syndiniales–diatom triangle (", cfg$label, ") ---\n", sep = "")
  
  d <- track_dirs(cfg)
  model_data <- loaded$model_data
  diatom_classes <- c("Bacillariophyceae", "Mediophyceae")
  
  # ── STEP 1: restrict to circos subset ─────────────────────────────
  circos_pairs <- pair_df %>%
    semi_join(links$class_links, by = c("Class_j", "Class_k"))
  
  target_otus <- unique(c(circos_pairs$OTU_j, circos_pairs$OTU_k))
  
  # ── STEP 2: build ASV table ───────────────────────────────────────
  asv_class <- bind_rows(
    pair_df %>% transmute(OTU = OTU_j, Class = Class_j),
    pair_df %>% transmute(OTU = OTU_k, Class = Class_k)
  ) %>%
    distinct() %>%
    mutate(group = case_when(
      Class == "Syndiniales" ~ "Syndiniales",
      Class == "Mediophyceae" ~ "Mediophyceae",
      Class == "Bacillariophyceae" ~ "Bacillariophyceae",
      TRUE ~ NA_character_
    )) %>%
    filter(
      OTU %in% target_otus,  
      !is.na(group)           
    )
  
  target_otus <- asv_class$OTU  
  
  # ── OPTIONAL: strong filtering ────────────────────────────────────
  if (strong_only) {
    keep <- circos_pairs %>%
      filter(strong_res,
             (Class_j == "Syndiniales" & Class_k %in% diatom_classes) |
               (Class_k == "Syndiniales" & Class_j %in% diatom_classes))
    
    target_otus <- intersect(target_otus, unique(c(keep$OTU_j, keep$OTU_k)))
    asv_class   <- asv_class %>% filter(OTU %in% target_otus)
  }
  
  if (length(target_otus) < 2) {
    cat("  Fewer than 2 valid ASVs — nothing to plot.\n")
    return(invisible(NULL))
  }
  
  # ── STEP 3: compute strength ordering ─────────────────────────────
  strength_tot <- circos_pairs %>%
    filter(OTU_j %in% target_otus, OTU_k %in% target_otus) %>%
    bind_rows(
      transmute(., OTU = OTU_j, a = abs(res_cor)),
      transmute(., OTU = OTU_k, a = abs(res_cor))
    ) %>%
    group_by(OTU) %>%
    summarise(tot = sum(a), .groups = "drop") %>%
    { setNames(.$tot, .$OTU) }
  
  ord <- function(v) v[order(-strength_tot[v])]
  
  synd_otus <- asv_class$OTU[asv_class$group == "Syndiniales"]
  medi_otus <- asv_class$OTU[asv_class$group == "Mediophyceae"]
  bacil_otus <- asv_class$OTU[asv_class$group == "Bacillariophyceae"]
  
  S <- c(ord(synd_otus), ord(medi_otus), ord(bacil_otus))
  group_of <- setNames(asv_class$group[match(S, asv_class$OTU)], S)
  
  cat("  ",
      length(synd_otus), " Syndiniales + ",
      length(medi_otus), " Mediophyceae + ",
      length(bacil_otus), " Bacillariophyceae = ",
      length(S), " total ASVs\n", sep = "")
  
  # ── STEP 4: labels ────────────────────────────────────────────────
  label_to_id <- setNames(names(model_data$otu_id_to_label),
                          model_data$otu_id_to_label)
  
  tax_full <- model_data$tax_mat[label_to_id[S], , drop = FALSE] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("OTU_ID") %>%
    mutate(OTU_Label = model_data$otu_id_to_label[OTU_ID])
  
  for (col in c("Order", "Family", "Genus", "Species"))
    if (!col %in% names(tax_full)) tax_full[[col]] <- NA_character_
  
  lab_lookup <- tax_full %>%
    mutate(label = case_when(
      !is.na(Species) & Species != "" ~ Species,
      !is.na(Genus)   & Genus   != "" ~ paste0(Genus, " sp."),
      !is.na(Family)  & Family  != "" ~ paste0(Family, " (fam.)"),
      !is.na(Order)   & Order   != "" ~ paste0(Order,  " (ord.)"),
      TRUE ~ Class
    ),
    label = paste0(label, " (ASV", sub(".*_ASV_", "", OTU_Label), ")"),
    label = gsub("Dino-Group", "DG", label)) %>%
    select(OTU_Label, label)
  
  S_labels <- lab_lookup$label[match(S, lab_lookup$OTU_Label)]
  
  # ── STEP 5: triangle data ─────────────────────────────────────────
  pos <- setNames(seq_along(S), S)
  
  P <- pair_df %>%
    filter(OTU_j %in% S, OTU_k %in% S) %>%
    mutate(
      i = pos[OTU_j],
      j = pos[OTU_k],
      x = pmin(i, j),
      y = pmax(i, j),
      quadrant = factor(quadrant, levels = c("Q1", "Q2", "Q3", "Q4"))
    )
  
  # ── STEP 6: diagonal colors ───────────────────────────────────────
  diag_df <- tibble(
    p = seq_along(S),
    col = case_when(
      group_of[S] == "Syndiniales" ~ "#7B2D8B",
      group_of[S] == "Mediophyceae" ~ "yellow4",
      group_of[S] == "Bacillariophyceae" ~ "darkgreen"
    )
  )
  
  quad_cols <- c(
    "Q1" = "#2166AC",
    "Q2" = "#F4A582",
    "Q3" = "#2AB5A0",
    "Q4" = "#B2182B"
  )
  
  # ── STEP 7: plot ──────────────────────────────────────────────────
  p_tri <- ggplot() +
    
    annotate("tile",
             x = diag_df$p, y = diag_df$p,
             fill = diag_df$col,
             width = 1, height = 1, colour = NA) +
    
    geom_text(
      data = diag_df,
      aes(x = p, y = p, label = S_labels),
      colour = "black",
      size = lab_size * 0.5,
      fontface = "bold",
      hjust = 0,
      nudge_x = 0.7
    ) +
    
    geom_tile(data = P,
              aes(x = x, y = y, fill = quadrant),
              colour = "grey95", linewidth = 0.1) +
    
    geom_point(data = tibble(group = c("Syndiniales", "Mediophyceae", "Bacillariophyceae")),
               aes(x = 1, y = 1, colour = group),
               size = 0.01,    
               alpha = 1) +   
    
    guides(colour = guide_legend(
      override.aes = list(size = 5, shape = 15, alpha = 1)
    )) +
  
    
    scale_fill_manual(values = quad_cols, name = "Decomposition quadrants") +
    
    scale_colour_manual(
      values = c(
        "Syndiniales" = "#7B2D8B",
        "Mediophyceae" = "yellow4",
        "Bacillariophyceae" = "darkgreen"
      ),
      name = "Group (diagonal)"
    ) +
    
    coord_fixed(clip = "off") +
    scale_x_continuous(expand = expansion(add = c(1, 8))) +   # extra room on the right
    scale_y_continuous(expand = expansion(add = c(1, 3))) +   # extra room on top/bottom
    
    labs(
      # title = paste0("Syndiniales–diatom co-occurrence — ", cfg$label),
      x = NULL, y = NULL
    ) +
    
    theme_minimal(base_size = 20) +
    theme(
      plot.title    = element_text(face = "bold"),
      panel.grid    = element_blank(),
      axis.text     = element_blank(),
      axis.ticks    = element_blank(),
      legend.position = "bottom",
      legend.text   = element_text(size = 18),
      legend.title  = element_text(size = 18 , face = "bold"),
      legend.key.size = unit(1, "cm"),
      plot.margin = margin(120, 40, 40, 40)   # top, right, bottom, left
    )
  
  print(p_tri)
  
  side <- max(8, length(S) * 0.22 + 4)
  
  ggsave(file.path(d$fig,
                   paste0("synd_diatom_triangle_", cfg$file_tag, ".png")),
         p_tri,
         width  = side + 6,    # extra room for labels spilling off the right
         height = side + 3,    # extra room for the bottom legend
         dpi    = 200,
         limitsize = FALSE)
  
  cat("  saved.\n")
  
  invisible(p_tri)
}



plot_synd_diatom_triangle(cfg_f, pair_df_f, loaded_f, links_f)
plot_synd_diatom_triangle(cfg_w, pair_df_w, loaded_w, links_w)


