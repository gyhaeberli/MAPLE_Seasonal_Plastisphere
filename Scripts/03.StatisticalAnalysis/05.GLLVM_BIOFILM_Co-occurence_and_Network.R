################################################################################
# GLLVM PRESENCE-ABSENCE ANALYSIS PIPELINE — SEASONAL FOULING COMMUNITY
#
# PRIMARY MODEL: best_inc_f_v2  (PA, M2 all-2-way interactions, fouling)
# NULL MODEL:    uncons_inc_f_v2 (PA, M7 unconstrained, fouling)
################################################################################

#  ENVIRONMENTAL AND RESIDUAL CO-OCCURRENCE 
################################################################################

library(gllvm)
library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(patchwork)
library(knitr)
library(scales)
library(vegan)
library(RColorBrewer)
library(tibble)
library(ggsci)
library(knitr)
library(kableExtra)
library(ggrepel)

cat("=== GLLVM ANALYSIS SCRIPT ===\n\n")

setwd("~/Github/MAPLE_Seasonal_Plastisphere/Scripts/03.StatisticalAnalysis")

SAVE_DIR_BASE  <- "~/Github/MAPLE_Seasonal_Plastisphere/Processed_data/gllvm_models"

OUT_DIR   <- "~/Github/MAPLE_Seasonal_Plastisphere/Results/gllvm_results"
TABLE_DIR <- file.path(OUT_DIR, "Tables")
FIG_DIR   <- file.path(OUT_DIR, "Figures")
for (d in c(OUT_DIR, FIG_DIR, TABLE_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# SECTION 0: LOAD DATA AND MODELS
# ==============================================================================


model_data <- readRDS(file.path(SAVE_DIR_BASE, "data_for_gllvm_FINAL.rds"))
metadata   <- model_data$metadata

patch_model <- function(fit, num_lv, include_X = TRUE) {
  fit$call$y           <- fit$y
  fit$call$num.lv      <- num_lv
  fit$call$studyDesign <- fit$studyDesign
  if (!is.null(fit$offset)) fit$call$offset <- fit$offset
  if (include_X && !is.null(fit$X)) fit$call$X <- fit$X
  fit
}

best_model <- patch_model(
  readRDS(file.path(SAVE_DIR_BASE, "GLLVM_final_incM2.rds")),
  num_lv = 2, include_X = TRUE
)

null_model <- patch_model(
  readRDS(file.path(SAVE_DIR_BASE, "GLLVM_null_incM7.rds")),
  num_lv = 2, include_X = FALSE
)

stopifnot(!is.null(best_model), !is.null(null_model))
cat("Models loaded successfully.\n\n")


################################################################################
# Pipeline:
#   3-A  Pollock decomposition — env vs residual correlation per OTU pair
#   3-B  Clustering of residual correlation matrix
#   3-C  Network analysis of strong residual associations
#   3-D  Hub / keystone taxa identification
################################################################################

library(gllvm)
library(igraph)
library(ggraph)
library(tidygraph)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(RColorBrewer)
library(pheatmap)
library(cluster)
library(scales)

# ── Global thresholds (adjust if needed) ─────────────────────────────────────
R_STRONG   <- 0.5   # |residual r| above this = "strong" association
N_CLUSTERS <- 5      # number of co-occurrence clusters (justified below)

cat("=== SECTION 3: RESIDUAL CO-OCCURRENCE ANALYSIS ===\n\n")


# ==============================================================================
# 3-A: POLLOCK DECOMPOSITION
#
# Background:
#   Raw OTU co-occurrence has two sources:
#     (1) Environmental filtering — two OTUs appear together because they both
#         prefer the same site/season/substrate. This is already captured by
#         the fixed effects in the GLLVM.
#     (2) Residual co-occurrence — two OTUs appear together (or avoid each
#         other) beyond what the environment predicts. This is what the latent
#         variables capture and what we interpret as potential biotic signal.
#
#   The Pollock (2014) framework makes this decomposition explicit by computing
#   both correlation types and plotting them against each other.
#
# Environmental correlation:
#   We multiply the coefficient matrix B (OTUs × predictors) by the design
#   matrix X (samples × predictors) to get fitted linear predictors — the
#   environmentally-predicted occurrence of each OTU at each site.
#   Correlating those fitted values across OTUs tells us how similarly two
#   OTUs respond to the MEASURED environment.
#
# Residual correlation:
#   getResidualCor() extracts this from the latent variables — the correlation
#   structure that remains AFTER the fixed effects are removed.
# ==============================================================================

cat("--- 3-A: Pollock decomposition ---\n\n")

# ── Compute environmental correlation ─────────────────────────────────────────

B     <- coef(best_model)$Xcoef   # OTUs × predictors
X_raw <- best_model$X

# Build the same dummy-coded design matrix the model used internally.
# The -1 drops the intercept because gllvm absorbs it into species intercepts.
X_numeric <- model.matrix(
  ~ site + season + substrate +
    site:season + site:substrate + season:substrate,
  data = X_raw
)[, -1]

# Safety check: column names must match B exactly
if (!identical(colnames(X_numeric), colnames(B))) {
  cat("Mismatch — only in X_numeric:", setdiff(colnames(X_numeric), colnames(B)), "\n")
  cat("Only in B:",                    setdiff(colnames(B), colnames(X_numeric)), "\n")
  stop("Design matrix columns do not match coefficient matrix. Check factor levels.")
}

# fitted_env: samples × OTUs — environmentally-predicted linear predictor
fitted_env <- X_numeric %*% t(B)

# env_cor: OTUs × OTUs — how similarly two OTUs respond to the measured env.
env_cor <- cor(fitted_env)
diag(env_cor) <- NA

# ── Compute residual correlation ───────────────────────────────────────────────

res_cor <- getResidualCor(best_model)
diag(res_cor) <- NA


# ── Compute trace reduction — how much of the residual covariance is removed
#    by the predictors? The trace is the sum of variances of all LVs.
#    A large reduction means environment explains most co-occurrence.
#    A small reduction means strong residual structure remains.

rcov_best <- getResidualCov(best_model, adjust = 0)
rcov_null <- getResidualCov(null_model, adjust = 0)
trace_reduction_pct <- round((1 - rcov_best$trace / rcov_null$trace) * 100, 1)

cat("Residual covariance trace — null model:", round(rcov_null$trace, 3), "\n")
cat("Residual covariance trace — best model:", round(rcov_best$trace, 3), "\n")
cat("Reduction by site + season + substrate: ", trace_reduction_pct, "%\n\n")
cat("Interpretation: the predictors explain", trace_reduction_pct,
    "% of the co-occurrence structure.\n",
    "The remaining", 100 - trace_reduction_pct,
    "% is residual and could reflect biotic interactions,\n",
    "unmeasured environmental gradients, or stochastic processes.\n\n")

# ── Build OTU-pair data frame ─────────────────────────────────────────────────

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

# Upper triangle only — each pair appears once
pairs <- which(upper.tri(res_cor), arr.ind = TRUE)

pair_df <- data.frame(
  OTU_j   = otu_labels[pairs[, 1]],
  OTU_k   = otu_labels[pairs[, 2]],
  res_cor  = res_cor[pairs],
  env_cor  = env_cor[pairs],
  stringsAsFactors = FALSE
) %>%
  left_join(tax_lookup %>% select(OTU_Label, Class_j = Class, Division_j = Division),
            by = c("OTU_j" = "OTU_Label")) %>%
  left_join(tax_lookup %>% select(OTU_Label, Class_k = Class, Division_k = Division),
            by = c("OTU_k" = "OTU_Label")) %>%
  mutate(
    same_class    = !is.na(Class_j) & !is.na(Class_k) & Class_j == Class_k,
    same_division = !is.na(Division_j) & !is.na(Division_k) & Division_j == Division_k,
    
    # Pollock quadrant: sign of env_cor × sign of res_cor
    # Q1 (++): both share habitat AND co-occur more → likely shared niche
    # Q2 (+-): share habitat BUT co-occur less → potential competition
    # Q3 (-+): different habitats BUT co-occur more → potential facilitation
    # Q4 (--): different habitats AND co-occur less → segregated, niche differentiation
    quadrant = case_when(
      env_cor >= 0 & res_cor >= 0 ~ "Q1",
      env_cor >= 0 & res_cor <  0 ~ "Q2",
      env_cor <  0 & res_cor >= 0 ~ "Q3",
      env_cor <  0 & res_cor <  0 ~ "Q4"
    ),
    
    # Surprise score: how much does residual exceed environmental expectation?
    # High surprise = the pair co-occurs (or avoids) far more than the
    # measured environment would predict → strongest candidate biotic signal.
    surprise   = abs(res_cor) - abs(env_cor),
    strong_res = abs(res_cor) > R_STRONG,
    concordant = sign(res_cor) == sign(env_cor)
  ) %>%
  filter(!is.na(res_cor), !is.na(env_cor))

cat("Total OTU pairs:", nrow(pair_df), "\n")
cat("Strong residual pairs (|r| >", R_STRONG, "):", sum(pair_df$strong_res), "\n\n")
cat("Strong residual pairs:", sum(pair_df$strong_res), 
    sprintf("(%.1f%%)\n", sum(pair_df$strong_res) / nrow(pair_df) * 100))


cat("Quadrant breakdown:\n")
print(table(pair_df$quadrant))
cat("\n")

write.csv(pair_df, file.path(TABLE_DIR, "S3A_pollock_pairs.csv"), row.names = FALSE)


# ── Pollock scatterplot ────────────────────────────────────────────────────────
# Each dot is one OTU pair. Blue = same taxonomic class, grey = cross-class.
# The diagonal line (slope = 1) marks where env_cor == res_cor.
# Points above the diagonal have MORE residual than environmental correlation —
# these are biologically the most interesting (biotic signal candidates).

# strong pairs color code
all_pollock <- ggplot(pair_df, aes(x = env_cor, y = res_cor)) +
  geom_hline(yintercept = 0, linewidth = 0.5, colour = "grey40") +
  geom_vline(xintercept = 0, linewidth = 0.5, colour = "grey40") +
  # All weak pairs — grey background
  geom_point(data = filter(pair_df, !strong_res),
             colour = "grey85", alpha = 0.3, size = 1.2) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  # Strong pairs coloured by quadrant
  geom_point(data = filter(pair_df, strong_res),
             aes(colour = quadrant),
             alpha = 0.75, size = 1.2) +
  annotate("text", x =  0.65, y =  0.9, label = "Q1",
           size = 2.7, colour = "grey25") +
  annotate("text", x =  0.65, y = -0.9, label = "Q2",
           size = 2.7, colour = "grey25") +
  annotate("text", x = -0.65, y =  0.9, label = "Q3",
           size = 2.7, colour = "grey25") +
  annotate("text", x = -0.65, y = -0.9, label = "Q4",
           size = 2.7, colour = "grey25") +
  scale_colour_manual(
    name   = "Quadrant",
    values = c(
      "Q1"      = "#2166AC",
      "Q2"    = "#F4A582",
      "Q3"    = "#2AB5A0",
      "Q4" = "#B2182B"
    ) ) +
  coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1)) +
  labs(
    # title    = "Pollock decomposition: environmental vs residual co-occurrence",
    # subtitle = paste0(
    #   "All pairs (n = ", nrow(pair_df), ") | ",
    #   "Coloured points: strong residual (|r| > ", R_STRONG,
    #   ", n = ", sum(pair_df$strong_res), ")\n",
    #    "Dashed diagonal = env_cor == res_cor"
    # ),
    x = "Design-based correlation\n", #(shared response to site + season + substrate)",
    y = "Residual correlation\n"#(co-occurrence beyond environmental filtering)"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "none")

print(all_pollock)
ggsave(file.path(FIG_DIR, "S3A_pollock_all_highlighted.png"),
       all_pollock, width = 8, height = 8, dpi = 200)
cat("  Pollock scatterplot all pairs highlighted saved.\n\n")


#all paris - association class level

p_pollock <- ggplot(pair_df, aes(x = env_cor, y = res_cor)) +
  geom_hline(yintercept = 0, linewidth = 0.5, colour = "grey40") +
  geom_vline(xintercept = 0, linewidth = 0.5, colour = "grey40") +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_point(data = filter(pair_df, !same_class),
             colour = "grey70", alpha = 0.25, size = 0.7) +
  geom_point(data = filter(pair_df, same_class),
             colour = "blue", alpha = 0.55, size = 0.8) +
  annotate("text", x =  0.65, y =  0.9, label = "Q1",
           size = 2.7, colour = "grey25") +
  annotate("text", x =  0.65, y = -0.9, label = "Q2",
           size = 2.7, colour = "grey25") +
  annotate("text", x = -0.65, y =  0.9, label = "Q3",
           size = 2.7, colour = "grey25") +
  annotate("text", x = -0.65, y = -0.9, label = "Q4",
           size = 2.7, colour = "grey25") +
  scale_colour_manual(name = "Pair type",
                      values = c("Same class" = "blue", "Cross-class" = "grey70")) +
  coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1)) +
  labs(
    title    = "Pollock decomposition: environmental vs residual co-occurrence",
    subtitle = paste0(
      "Each point = one OTU pair (n = ", nrow(pair_df), ")\n",
      "Blue = same taxonomic class | Dashed diagonal = env_cor == res_cor\n"
    ),
    x = "Desing-based correlation\n(shared response to site + season + substrate)",
    y = "Residual correlation\n(co-occurrence beyond environmental filtering)"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

print(p_pollock)

ggsave(file.path(FIG_DIR, "S3A_pollock_scatterplot.png"),
       p_pollock, width = 8, height = 8, dpi = 200)
cat("  Pollock scatterplot all pairs saved.\n\n")



# ── Quadrant summary table ─────────────────────────────────────────────────────

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

cat("Quadrant summary:\n")
print(quadrant_summary, n = Inf)
write.csv(quadrant_summary, file.path(TABLE_DIR, "S3A_quadrant_summary.csv"), row.names = FALSE)
cat("\n")

kable(quadrant_summary)%>%
  kable_styling(latex_options = "striped")


# ── Build quadrant summary table ─────────────────────────────────────────────

library(flextable)

quad_sum = flextable(quadrant_summary) %>%
  set_header_labels(
    quadrant       = "Quadrant",
    n_pairs        = "N pairs",
    pct_pairs      = "% pairs",
    n_strong       = "N strong",
    pct_strong     = "% strong",
    mean_res_cor   = "Mean res. cor.",
    mean_env_cor   = "Mean env. cor.",
    mean_surprise  = "Mean CoOc Score",
    n_same_class   = "N same class",
    pct_same_class = "% same class"
  ) %>%
  theme_vanilla() %>%
  autofit()
save_as_docx(quad_sum, path = file.path(TABLE_DIR, "quadrant_summary_biof.docx"))


# ── Top surprising pairs ───────────────────────────────────────────────────────

top_surprise <- pair_df %>%
  filter(strong_res) %>%
  arrange(desc(surprise)) %>%
  slice_head(n = 100) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

cat("Top 30 most surprising pairs (high residual, low environmental):\n")
print(top_surprise %>% select(OTU_j, OTU_k, Class_j, Class_k,
                              env_cor, res_cor, surprise, quadrant))
write.csv(top_surprise, file.path(TABLE_DIR, "S3A_top_surprising_pairs.csv"), row.names = FALSE)
cat("\n")



# high surprise score: res_corr >> env_corr

p_pollock_surprise <- ggplot() +
  geom_point(data = pair_df,
             aes(x = env_cor, y = res_cor),
             colour = "grey80", alpha = 0.25, size = 0.6) +
  geom_point(data = filter(pair_df, strong_res),
             aes(x = env_cor, y = res_cor),
             colour = "grey50", alpha = 0.4, size = 0.9) +
  # Top 30 — single highlight colour, no interaction typing
  geom_point(data = top_surprise,
             aes(x = env_cor, y = res_cor),
             colour = "orange", size = 3, alpha = 0.9) +
  ggrepel::geom_text_repel(
    data = top_surprise,
    aes(x = env_cor, y = res_cor,
        label = paste0(sub("_ASV_.*", "", OTU_j),
                       " × ",
                       sub("_ASV_.*", "", OTU_k))),
    colour       = "darkred",
    size         = 4,
    max.overlaps = 50,
    segment.size = 0.3,
    segment.alpha = 0.5,
    box.padding  = 0.4
  ) +
  geom_hline(yintercept = 0, linewidth = 0.5, colour = "grey30") +
  geom_vline(xintercept = 0, linewidth = 0.5, colour = "grey30") +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_hline(yintercept =  R_STRONG, linetype = "dotted",
             colour = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = -R_STRONG, linetype = "dotted",
             colour = "grey40", linewidth = 0.4) +
  annotate("text", x =  0.65, y =  0.95,
           label = "Q1",
           size = 2.8, colour = "grey30", hjust = 0.5) +
  annotate("text", x =  0.65, y = -0.95,
           label = "Q2",
           size = 2.8, colour = "grey30", hjust = 0.5) +
  annotate("text", x = -0.65, y =  0.95,
           label = "Q3",
           size = 2.8, colour = "grey30", hjust = 0.5) +
  annotate("text", x = -0.65, y = -0.95,
           label = "Q4",
           size = 2.8, colour = "grey30", hjust = 0.5) +
  coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1)) +
  labs(
  #   title    = "Pollock decomposition: environmental vs residual co-occurrence",
  #   subtitle = paste0(
  #     "Grey = all pairs (n = ", nrow(pair_df), ") | ",
  #     "Dark grey = strong residual |r| > ", R_STRONG,
  #     " (n = ", sum(pair_df$strong_res), ")\n",
  #     "Red = top 30 most surprising pairs (highest |res_cor| − |env_cor|)\n",
  #     "Dotted lines = ±", R_STRONG, " threshold | ",
  #     "Dashed diagonal = env_cor == res_cor"
  #   ),
    x = "Desing-based correlation\n(shared response to site + season + substrate)",
    y = "Residual correlation\n(co-occurrence beyond environmental filtering)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(size = 8)
  )

ggsave(file.path(FIG_DIR, "S3A_pollock_surprise_highlighted.png"),
       p_pollock_surprise, width = 10, height = 10, dpi = 200)

print(p_pollock_surprise)


#circos plots ################################################################

library(circlize)
library(dplyr)
library(RColorBrewer)

# STEP 1: PREPARE CLASS-LEVEL AGGREGATION FROM TOP 50 PAIRS
#
# Each pair has a Class_j and Class_k. We want to know, for each unique
# class-to-class combination, the mean residual correlation and whether the
# association is positive or negative. We also count how many pairs contribute
# to each class-class link so we can optionally scale by n_pairs.


#remove unclassified at calss level:
top_surprise <- pair_df %>%
  filter(strong_res) %>%
  filter(!is.na(Class_j), Class_j != "",
         !is.na(Class_k), Class_k != "") %>%
  arrange(desc(surprise)) %>%
  slice_head(n = 50) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

# Replace NA class with "Unclassified" so they appear in the plot
top30_clean <- top_surprise %>%
  mutate(
    Class_j = if_else(is.na(Class_j) | Class_j == "", "Unclassified", Class_j),
    Class_k = if_else(is.na(Class_k) | Class_k == "", "Unclassified", Class_k),
    association = if_else(res_cor > 0, "positive", "negative")
  )

# Aggregate to class level
# For pairs within the same class (Class_j == Class_k), we still include them
# — they will appear as self-links on the circos plot.
class_links <- top30_clean %>%
  group_by(Class_j, Class_k, association) %>%
  summarise(
    mean_res_cor  = round(mean(abs(res_cor)), 3),
    mean_surprise = round(mean(surprise), 3),
    n_pairs       = n(),
    # Dominant quadrant for this class pair
    dominant_quadrant = names(sort(table(quadrant), decreasing = TRUE))[1],
    .groups = "drop"
  ) %>%
  # Link width will be proportional to mean |res_cor| × n_pairs
  # This weights strong AND frequent class-level associations higher
  mutate(link_weight = mean_res_cor * n_pairs)

cat("Class-level links:\n")
print(class_links, n = Inf)
cat("\n")


# STEP 2: DEFINE COLOUR PALETTES
#
# Each unique class gets a sector colour.
# Links are coloured by association type (positive = warm, negative = cool)
# with transparency so overlapping links are readable.

all_classes <- sort(unique(c(class_links$Class_j, class_links$Class_k)))
n_classes   <- length(all_classes)

# Use a qualitative palette — Set1 + Set2 combined for up to ~16 classes
if (n_classes <= 8) {
  sector_cols <- setNames(brewer.pal(n_classes, "Set2"), all_classes)
} else if (n_classes <= 12) {
  sector_cols <- setNames(brewer.pal(n_classes, "Paired"), all_classes)
} else {
  sector_cols <- setNames(
    colorRampPalette(brewer.pal(12, "Paired"))(n_classes),
    all_classes
  )
}

# Link colours: positive = red family, negative = blue family
# add2 = transparency (0 = opaque, 255 = fully transparent)
link_col_positive <- adjustcolor("darkblue" , alpha.f = 0.6)
link_col_negative <- adjustcolor("red2" , alpha.f = 0.7)

# Quadrant colours for Plot 2
quadrant_cols <- c(
  "Q1"        = adjustcolor( "#2166AC", alpha.f = 0.65),
  "Q2"      = adjustcolor("#F4A582", alpha.f = 0.65),
  "Q3"      = adjustcolor("#2AB5A0", alpha.f = 0.65),
  "Q4"   = adjustcolor("#B2182B", alpha.f = 0.65)
)

# STEP 3: PLOT 1 — CIRCOS BY ASSOCIATION TYPE (POSITIVE / NEGATIVE)

CANVAS_SIZE  <- 12      # inches — increase further if still cramped
PLOT_DPI     <- 300
GAP_DEGREES  <- 7     # degrees of gap between sectors
LABEL_CEX    <- 1.5  # font size for sector labels
LABEL_OFFSET <- 0.5     # how far labels sit outside the track (in mm units)

# Helper that draws the label track identically for both plots
draw_label_track <- function() {
  circos.trackPlotRegion(
    track.index = 1,
    track.height = 0.08,          # taller track = more room for labels
    panel.fun = function(x, y) {
      xlim   <- get.cell.meta.data("xlim")
      ylim   <- get.cell.meta.data("ylim")
      sector <- get.cell.meta.data("sector.index")
      # Place text at outer edge of track
      circos.text(
        mean(xlim),
        ylim[2] + LABEL_OFFSET,   # push further outside the circle
        sector,
        facing     = "clockwise",
        niceFacing = TRUE,
        adj        = c(0, 0.5),
        cex        = LABEL_CEX,
        font       = 2
      )
    },
    bg.border = NA
  )
}

#
# How to read this plot:
#   - Each arc (sector) on the circle = one taxonomic class
#   - Sector width = total link weight involving that class
#     (how much residual co-occurrence involves this class)
#   - Each ribbon (link) = one class-to-class association
#   - Link width at each end = proportional to mean |res_cor| × n_pairs
#   - Red links = positive residual co-occurrence (taxa appear together
#     more than environment predicts — facilitation / parasitism candidates)
#   - Blue links = negative residual co-occurrence (taxa avoid each other
#     more than environment predicts — competition / predation candidates)

png(file.path(FIG_DIR, "S3_circos_top50_by_association.png"),
    width  = CANVAS_SIZE,
    height = CANVAS_SIZE,
    units  = "in",
    res    = PLOT_DPI)

# Large margins: top/right/bottom/left
# Bottom margin is large so legend text sits below the circle
par(mar = c(8, 2, 4, 2))

circos.clear()
circos.par(
  gap.after               = GAP_DEGREES,
  start.degree            = 90,
  clock.wise              = TRUE,
  track.margin            = c(0.01, 0.05),   # inner, outer margin of tracks
  # Shrink the plot region so labels have room outside
  canvas.xlim             = c(-1.4, 1.4),
  canvas.ylim             = c(-1.4, 1.4),
  points.overflow.warning = FALSE
)

chordDiagram(
  x = class_links %>%
    select(from = Class_j, to = Class_k, value = link_weight) %>%
    as.data.frame(),
  grid.col              = sector_cols,
  col                   = ifelse(class_links$association == "positive",
                                 link_col_positive, link_col_negative),
  transparency          = 0,
  directional           = 0,
  annotationTrack       = "grid",
  annotationTrackHeight = 0.04,
  link.border           = NA,
  order                 = all_classes,
  link.sort             = TRUE,
  link.decreasing       = TRUE,
  # Reduce the circle radius slightly to leave label space
  preAllocateTracks     = list(track.height = 0.08)
)

draw_label_track()

# Legends placed BELOW the plot using par coordinates
# xpd = TRUE allows drawing outside the plot region
par(xpd = TRUE)

legend(
  x      = -1.35,           # left side
  y      = -1.55,           # below circle
  legend = c("Positive co-occurrence", "Negative co-occurence"),
  fill   = c(link_col_positive, link_col_negative),
  border = NA,
  title  = "Association type (link colour)",
  bty    = "n",
  cex    = 0.9,
  title.cex = 0.9
)

# Add width explanation below
legend(
  x      = -0.1,           # right side
  y      = -1.55,
  legend = "Link width = mean |residual cor.| × N OTU pairs",
  bty    = "n",
  cex    = 0.9,
  title  = "Link width",
  title.cex = 0.9
)

circos.clear()
dev.off()
cat("Plot 1 saved: S3_circos_top50_by_association.png\n")




# STEP 4: PLOT 2 — CIRCOS COLOURED BY QUADRANT
#
# Same structure as Plot 1 but links are coloured by their dominant Pollock
# quadrant. 

png(file.path(FIG_DIR, "S3_circos_top50_by_quadrant.png"),
    width  = CANVAS_SIZE,
    height = CANVAS_SIZE,
    units  = "in",
    res    = PLOT_DPI)

par(mar = c(8, 2, 4, 2))

circos.clear()
circos.par(
  gap.after               = GAP_DEGREES,
  start.degree            = 90,
  clock.wise              = TRUE,
  track.margin            = c(0.01, 0.05),
  canvas.xlim             = c(-1.4, 1.4),
  canvas.ylim             = c(-1.4, 1.4),
  points.overflow.warning = FALSE
)

chordDiagram(
  x = class_links %>%
    select(from = Class_j, to = Class_k, value = link_weight) %>%
    as.data.frame(),
  grid.col              = sector_cols,
  col                   = quadrant_cols[class_links$dominant_quadrant],
  transparency          = 0,
  directional           = 0,
  annotationTrack       = "grid",
  annotationTrackHeight = 0.04,
  link.border           = NA,
  order                 = all_classes,
  link.sort             = TRUE,
  link.decreasing       = TRUE,
  preAllocateTracks     = list(track.height = 0.08)
)

draw_label_track()

par(xpd = TRUE)

legend(
  x      = -1.35,
  y      = -1.55,
  legend = c(
    "Q1",
    "Q2",
    "Q3",
    "Q4"
  ),
  fill   = unname(quadrant_cols),
  border = NA,
  title  = "Pollock quadrant (link colour)",
  bty    = "n",
  cex    = 0.82,
  title.cex = 0.9
)


circos.clear()
dev.off()
cat("Plot 2 saved: S3_circos_top50_by_quadrant.png\n")



# STEP 5: SUMMARY TABLE — what drives each class-level link
# Printed to console so you can cross-reference with the circos plots

cat("\n=== Class-level link summary ===\n\n")

class_links %>%
  arrange(desc(link_weight)) %>%
  select(
    From          = Class_j,
    To            = Class_k,
    Association   = association,
    Mean_res_cor  = mean_res_cor,
    Mean_surprise = mean_surprise,
    N_pairs       = n_pairs,
    Link_weight   = link_weight,
    Dominant_quadrant = dominant_quadrant
  ) %>%
  print(n = Inf)


kable(class_links) %>%
  kable_styling(latex_options = "striped")

class_links_dedup <- class_links %>%
  rowwise() %>%
  mutate(
    pair_key = paste(sort(c(Class_j, Class_k)), collapse = "___")
  ) %>%
  ungroup() %>%
  distinct(pair_key, .keep_all = TRUE) %>%
  select(-pair_key)


kable(class_links_dedup) %>%
  kable_styling(latex_options = "striped")


#with lowest taxonomic assignment

# All taxonomy columns available in your taxonomy matrix, from broadest to finest
tax_cols_full <- intersect(
  c("Division", "Class", "Order", "Family", "Genus", "Species"),
  colnames(model_data$tax_mat)
)

label_to_id <- setNames(
  names(model_data$otu_id_to_label),
  model_data$otu_id_to_label
)

# Get all unique OTUs appearing in the top 30 pairs
all_otus_top30 <- unique(c(top_surprise$OTU_j, top_surprise$OTU_k))
otu_ids_top30  <- label_to_id[all_otus_top30]

# Build full taxonomy table for these OTUs
tax_full <- model_data$tax_mat[otu_ids_top30, tax_cols_full, drop = FALSE] %>%
  as.data.frame() %>%
  rownames_to_column("OTU_ID") %>%
  mutate(OTU_Label = model_data$otu_id_to_label[OTU_ID])

annotated_pairs <- top_surprise %>%
  left_join(tax_full %>% select(OTU_Label, 
                                Div_j = Division, Cls_j = Class,
                                Ord_j = Order,    Fam_j = Family,
                                Gen_j = Genus,    Spp_j = Species),
            by = c("OTU_j" = "OTU_Label")) %>%
  left_join(tax_full %>% select(OTU_Label,
                                Div_k = Division, Cls_k = Class,
                                Ord_k = Order,    Fam_k = Family,
                                Gen_k = Genus,    Spp_k = Species),
            by = c("OTU_k" = "OTU_Label")) %>%
  mutate(
    lowest_j = case_when(
      !is.na(Spp_j) & Spp_j != "" ~ Spp_j,
      !is.na(Gen_j) & Gen_j != "" ~ paste0(Gen_j, " sp."),
      !is.na(Fam_j) & Fam_j != "" ~ paste0(Fam_j, " (fam.)"),
      !is.na(Ord_j) & Ord_j != "" ~ paste0(Ord_j, " (ord.)"),
      TRUE                        ~ Cls_j
    ),
    lowest_k = case_when(
      !is.na(Spp_k) & Spp_k != "" ~ Spp_k,
      !is.na(Gen_k) & Gen_k != "" ~ paste0(Gen_k, " sp."),
      !is.na(Fam_k) & Fam_k != "" ~ paste0(Fam_k, " (fam.)"),
      !is.na(Ord_k) & Ord_k != "" ~ paste0(Ord_k, " (ord.)"),
      TRUE                        ~ Cls_k
    )
  ) %>%
  select(OTU_j, lowest_j, Cls_j, OTU_k, lowest_k, Cls_k,
         env_cor, res_cor, surprise, quadrant)

print(annotated_pairs)


kable(annotated_pairs) %>%
  kable_styling(latex_options = "striped")


## save table

csv_table = annotated_pairs %>%
  rename(CoOc_score = surprise)

write.csv(annotated_pairs, file.path(TABLE_DIR, "S3A_annotated_top50_pairs.csv"), row.names = FALSE)



library(flextable)
library(dplyr)

quadrant_colors <- c(
  "Q1" = "#CCE5FF",
  "Q2" = "#FFE0B2",
  "Q3" = "#D5F5E3",
  "Q4" = "#FFCCCC"
)

q_labels <- c(
  "Q1",
  "Q2",
  "Q3",
  "Q4"
)

table_data <- class_links_dedup %>%
  mutate(
    `Class pair`        = paste(Class_j, "—", Class_k),
    `Association`       = association,
    `N OTU pairs`       = n_pairs,
    `Link weight`       = round(link_weight, 3),
    `Mean resid. corr.` = round(mean_res_cor, 3),
    `Mean score`        = round(mean_surprise, 3),
    Q                   = substr(trimws(as.character(dominant_quadrant)), 1, 2)
  ) %>%
  arrange(Q, desc(`Link weight`)) %>%
  select(`Class pair`, `Association`, `N OTU pairs`, `Link weight`,
         `Mean resid. corr.`, `Mean score`, Q)

ft <- flextable(table_data) %>%
  theme_vanilla() %>%
  bold(j = "Class pair") %>%
  align(j = c("N OTU pairs", "Link weight", "Mean resid. corr.", "Mean score"),
        align = "right", part = "body") %>%
  align(j = c("N OTU pairs", "Link weight", "Mean resid. corr.", "Mean score"),
        align = "right", part = "header")

# Colour rows by quadrant
for (q in c("Q1", "Q2", "Q3", "Q4")) {
  rows <- which(table_data$Q == q)
  if (length(rows) == 0) next
  ft <- ft %>%
    bg(i = rows, bg = quadrant_colors[q], part = "body")
}

# Add quadrant group separator rows
# flextable doesn't have pack_rows, so we add a bold label row before each group
ft <- ft %>%
  autofit()

ft






# ==============================================================================
# DIAGNOSTIC: ARE STRONG RESIDUAL PAIRS DISTRIBUTED ACROSS CONDITIONS?
#
# This is an informal check — not a formal statistical test.
# We take the top N residual pairs from the GLLVM decomposition and ask:
# within each level of site / season / substrate, do these pairs still
# co-occur in the raw PA data?
#
# What this tells you:
#   - A pair that co-occurs in ALL or MOST condition levels is likely a
#     genuine cross-condition association.
#   - A pair that only co-occurs in ONE level is a candidate for being
#     condition-specific (even if the GLLVM residual is dataset-wide).
#
# What this does NOT tell you:
#   - Whether the association is statistically significant within each level.
#   - Whether the pattern reflects biotic interaction vs. unmeasured gradient.
#   - This uses raw PA co-occurrence, not GLLVM-corrected residuals.
# ==============================================================================

cat("--- DIAGNOSTIC: Condition distribution of top residual pairs ---\n\n")

N_DIAG_PAIRS <- 50   # how many top pairs to inspect

# ── Select top pairs by residual correlation strength ─────────────────────────

diag_pairs <- pair_df %>%
  filter(strong_res) %>%
  arrange(desc(abs(res_cor))) %>%
  slice_head(n = N_DIAG_PAIRS)

cat("Checking top", N_DIAG_PAIRS, "pairs by |residual r|\n\n")

# ── Build full PA matrix with metadata attached ────────────────────────────────
# We use the full otu_filt (pre top-200 subsetting) to maximise sample coverage,
# but restrict columns to OTUs actually in the model.

pa_full <- (model_data$otu_filt[, model_data$top_otus] > 0) * 1L
colnames(pa_full) <- model_data$otu_id_to_label[model_data$top_otus]

meta_diag <- model_data$metadata   # already filtered and factor-levelled

stopifnot(identical(rownames(pa_full), rownames(meta_diag)))

# ── For each condition variable, compute within-level co-occurrence rate ───────

condition_vars <- c("site", "season", "substrate")

diag_results <- lapply(condition_vars, function(cvar) {
  
  levels_cvar <- levels(meta_diag[[cvar]])
  
  pair_rows <- lapply(seq_len(nrow(diag_pairs)), function(i) {
    
    otu_j <- diag_pairs$OTU_j[i]
    otu_k <- diag_pairs$OTU_k[i]
    
    # Check both OTUs are present in the PA matrix columns
    if (!otu_j %in% colnames(pa_full) || !otu_k %in% colnames(pa_full)) {
      return(NULL)
    }
    
    level_rows <- lapply(levels_cvar, function(lv) {
      
      idx <- which(meta_diag[[cvar]] == lv)
      if (length(idx) < 3) return(NULL)   # skip tiny groups
      
      pj <- pa_full[idx, otu_j]
      pk <- pa_full[idx, otu_k]
      
      n_samples      <- length(idx)
      n_both_present <- sum(pj == 1 & pk == 1)
      n_either       <- sum(pj == 1 | pk == 1)
      
      # Jaccard: co-occurrence relative to how often either appears
      # Ranges 0 (never together) to 1 (always together when either present)
      jaccard <- if (n_either > 0) n_both_present / n_either else NA_real_
      
      # Phi coefficient: correlation between two binary vectors
      # Equivalent to Pearson r on binary data
      # Positive = co-occur more than chance, negative = avoid
      phi <- if (sd(pj) > 0 && sd(pk) > 0) cor(pj, pk) else NA_real_
      
      data.frame(
        condition_var  = cvar,
        condition_level = lv,
        OTU_j          = otu_j,
        OTU_k          = otu_k,
        res_cor        = diag_pairs$res_cor[i],
        n_samples      = n_samples,
        n_both         = n_both_present,
        n_either       = n_either,
        jaccard        = round(jaccard, 3),
        phi            = round(phi, 3),
        stringsAsFactors = FALSE
      )
    })
    
    do.call(rbind, level_rows)
  })
  
  do.call(rbind, pair_rows)
})

diag_df <- do.call(rbind, diag_results)

# ── Summary: for each pair, how many condition levels show co-occurrence? ──────
# A pair is "present" in a level if Jaccard > 0 (they co-occur at least once).
# We also flag pairs where phi is consistently positive across levels.

pair_summary <- diag_df %>%
  group_by(condition_var, OTU_j, OTU_k, res_cor) %>%
  summarise(
    n_levels_total    = n(),
    n_levels_cooccur  = sum(jaccard > 0,    na.rm = TRUE),
    n_levels_positive = sum(phi > 0,        na.rm = TRUE),
    n_levels_strong   = sum(phi > 0.2,      na.rm = TRUE),
    mean_jaccard      = round(mean(jaccard, na.rm = TRUE), 3),
    mean_phi          = round(mean(phi,     na.rm = TRUE), 3),
    min_phi           = round(min(phi,      na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  mutate(
    coverage = round(n_levels_cooccur / n_levels_total, 2),
    # Flag pairs that appear in all levels (fully distributed)
    # vs. those concentrated in one level
    distribution = case_when(
      coverage == 1.00             ~ "All levels",
      coverage >= 0.50             ~ "Most levels",
      coverage >  0                ~ "Few levels",
      TRUE                         ~ "Absent"
    )
  ) %>%
  arrange(condition_var, desc(mean_phi))

cat("=== Pair distribution summary by condition variable ===\n\n")
print(pair_summary, n = Inf)



# ── Console summary ───────────────────────────────────────────────────────────

cat("=== DISTRIBUTION FLAGS ===\n\n")

for (cv in condition_vars) {
  cat("Condition:", cv, "\n")
  tbl <- pair_summary %>%
    filter(condition_var == cv) %>%
    count(distribution) %>%
    mutate(pct = round(n / sum(n) * 100, 1))
  print(tbl)
  cat("\n")
}

# cat("NOTE: 'All levels' means the pair co-occurs at least once in every\n")
# cat("condition level in the raw PA data. It does NOT mean the association\n")
# cat("strength is equal across levels — inspect phi values per level for that.\n\n")












# ==============================================================================
# 3-B: CLUSTERING OF RESIDUAL CORRELATION MATRIX
#
# We cluster OTUs by their residual co-occurrence profiles — how each OTU
# co-varies with all other OTUs after environment is accounted for.
# OTUs in the same cluster tend to co-occur together (positive within-cluster
# residual correlations) and may share ecological roles or dependencies.
#
# Method: Ward.D2 hierarchical clustering on distance = 1 - residual_r.
# Cluster number selected by comparing silhouette width and cluster balance.
# ==============================================================================

cat("--- 3-B: Residual correlation clustering ---\n\n")

dist_mat   <- as.dist(1 - getResidualCor(best_model))
hclust_res <- hclust(dist_mat, method = "ward.D2")

# ── Diagnostic 1: silhouette width across k = 2..10 ──────────────────────────
# Silhouette width measures how well each OTU fits its assigned cluster.
# It ranges from -1 (wrong cluster) to +1 (perfectly clustered).
# Pick the k that maximises average silhouette across all OTUs.

sil_widths <- sapply(2:10, function(k) {
  cl  <- cutree(hclust_res, k = k)
  sil <- silhouette(cl, dist_mat)
  mean(sil[, "sil_width"])
})

best_k_sil <- which.max(sil_widths) + 1
cat("Silhouette-optimal k:", best_k_sil, "\n")

sil_df <- data.frame(k = 2:10, avg_silhouette = round(sil_widths, 4))
print(sil_df)

# ── Diagnostic 2: cluster balance at each k ───────────────────────────────────
# Silhouette can favour k=2 even when one cluster has 95% of OTUs.
# Check balance before committing.

cat("\nCluster size balance (k = 2 to 8):\n")
for (k in 2:8) {
  cl    <- cutree(hclust_res, k = k)
  sizes <- sort(table(cl), decreasing = TRUE)
  pcts  <- round(sizes / sum(sizes) * 100, 1)
  cat("  k =", k, ":", paste(pcts, collapse = "% / "), "%\n")
}
cat("\n")
cat("Selected N_CLUSTERS =", N_CLUSTERS,
    "(set at top of script — adjust based on diagnostics above)\n\n")

# ── Cut tree and annotate clusters ────────────────────────────────────────────

clusters   <- cutree(hclust_res, k = N_CLUSTERS)
cluster_pal <- setNames(
  RColorBrewer::brewer.pal(max(N_CLUSTERS, 3), "Set2")[1:N_CLUSTERS],
  paste0("C", seq_len(N_CLUSTERS))
)

cluster_df <- data.frame(
  OTU_Label = otu_labels,
  OTU_ID    = otu_ids,
  Cluster   = paste0("C", clusters),
  stringsAsFactors = FALSE
) %>%
  left_join(
    model_data$tax_mat[otu_ids, tax_cols_use, drop = FALSE] %>%
      as.data.frame() %>%
      rownames_to_column("OTU_ID"),
    by = "OTU_ID"
  )

cat("Cluster sizes:\n")
print(table(cluster_df$Cluster))
cat("\n")

cat("Top 5 taxonomic classes per cluster:\n")
cluster_df %>%
  mutate(Class = replace_na(Class, "Unclassified")) %>%
  count(Cluster, Class) %>%
  group_by(Cluster) %>%
  slice_max(n, n = 5, with_ties = FALSE) %>%
  mutate(prop = round(n / sum(n), 2)) %>%
  print(n = Inf)
cat("\n")

write.csv(cluster_df, file.path(TABLE_DIR, "S3B_cluster_taxonomy.csv"), row.names = FALSE)

# ── Residual correlation heatmap ──────────────────────────────────────────────

cluster_order <- order(clusters)
res_cor_plot  <- getResidualCor(best_model)

annot_heatmap <- data.frame(
  Cluster  = factor(cluster_df$Cluster[cluster_order],
                    levels = paste0("C", seq_len(N_CLUSTERS))),
  Division = replace_na(cluster_df$Division[cluster_order], "Unknown"),
  row.names = otu_labels[cluster_order]
)

n_div <- length(unique(annot_heatmap$Division))
annot_colours <- list(
  Cluster  = cluster_pal,
  Division = setNames(
    colorRampPalette(brewer.pal(9, "Set1"))(n_div),
    unique(annot_heatmap$Division)
  )
)

png(file.path(FIG_DIR, "S3B_residual_heatmap.png"),
    width = 16, height = 14, units = "in", res = 300)
pheatmap(
  res_cor_plot[cluster_order, cluster_order],
  color             = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
  breaks            = seq(-1, 1, length.out = 102),
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  gaps_row          = cumsum(table(sort(clusters))),
  gaps_col          = cumsum(table(sort(clusters))),
  annotation_row    = annot_heatmap,
  annotation_col    = annot_heatmap,
  annotation_colors = annot_colours,
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  treeheight_row    = 0,
  treeheight_col    = 0,
  main              = paste0(
    "Residual OTU co-occurrence | ", trace_reduction_pct,
    "% of co-occurrence variation explained by design | k = ", N_CLUSTERS, " clusters"
  ),
  fontsize     = 9,
  border_color = NA,
  legend_breaks  = c(-1, -0.5, 0, 0.5, 1),
  legend_labels  = c("-1", "-0.5", "0", "0.5", "1")
)
dev.off()
cat("Residual heatmap saved.\n\n")







# ==============================================================================
# 3-C: NETWORK ANALYSIS OF STRONG RESIDUAL ASSOCIATIONS
#
# Background:
#   We build a network where each OTU is a node and each strong residual
#   co-occurrence (|r| > R_STRONG) is an edge. Positive edges = OTUs that
#   co-occur more than the environment predicts (facilitation / shared
#   unmeasured niche). Negative edges = OTUs that avoid each other beyond
#   environmental expectation (competition / antagonism candidates).
#
# Network metrics computed per node:
#   - Degree: number of strong residual associations. High-degree OTUs are
#     broadly connected and potentially important in structuring the community.
#   - Strength: sum of |r| across all edges. Weights how strong the connections
#     are, not just how many.
#   - Betweenness centrality: how often an OTU lies on shortest paths between
#     other OTUs. High betweenness = "bridge" taxa linking otherwise separate
#     groups — ecologically interpreted as potential keystone species.
#   - Positive degree / negative degree: separately counts positive (potential
#     facilitation) and negative (potential competition/predation) associations.
# ==============================================================================

cat("--- 3-C: Network analysis ---\n\n")


# new threshold for network
R_STRONG = 0.5

# ── Build edge list ───────────────────────────────────────────────────────────

edges <- pair_df %>%
  filter(abs(res_cor) > R_STRONG) %>%    # filter directly on value, not the flag
  select(OTU_j, OTU_k, res_cor, env_cor, surprise, quadrant,
         same_class, Class_j, Class_k) %>%
  mutate(
    association = if_else(res_cor > 0, "positive", "negative"),
    abs_r       = abs(res_cor)
  )

cat("Edges in network (|r| >", R_STRONG, "):", nrow(edges), "\n")
cat("  Positive associations:", sum(edges$association == "positive"), "\n")
cat("  Negative associations:", sum(edges$association == "negative"), "\n\n")

# if (nrow(edges) == 0) {
#   cat("No strong residual pairs found. Lower R_STRONG threshold and re-run.\n\n")
# } else {

# ── Build node metadata ─────────────────────────────────────────────────────

node_otus <- unique(c(edges$OTU_j, edges$OTU_k))

node_df <- tax_lookup %>%
  filter(OTU_Label %in% node_otus) %>%
  mutate(
    Class_plot    = if_else(is.na(Class) | Class == "", "Unclassified", Class),
    Division_plot = if_else(is.na(Division) | Division == "", "Unclassified", Division)
  ) %>%
  left_join(
    cluster_df %>% select(OTU_Label, Cluster),
    by = "OTU_Label") 

# ── Build igraph object ─────────────────────────────────────────────────────

g <- graph_from_data_frame(
  d        = edges %>% select(from = OTU_j, to = OTU_k,
                              weight = abs_r, association, res_cor),
  directed = FALSE,
  vertices = node_df %>% select(name = OTU_Label, everything())
)


# ── Compute node-level metrics ──────────────────────────────────────────────

V(g)$node_degree      <- as.numeric(igraph::degree(g))
V(g)$node_strength    <- as.numeric(igraph::strength(g, weights = E(g)$weight))
V(g)$node_betweenness <- as.numeric(igraph::betweenness(g, 
                                                        weights = 1 / E(g)$weight,
                                                        normalized = TRUE))
V(g)$pos_degree <- sapply(V(g)$name, function(otu) {
  sum(edges$association == "positive" &
        (edges$OTU_j == otu | edges$OTU_k == otu))
})
V(g)$neg_degree <- sapply(V(g)$name, function(otu) {
  sum(edges$association == "negative" &
        (edges$OTU_j == otu | edges$OTU_k == otu))
})

cat("Network summary:\n")
cat("  Nodes:", vcount(g), "\n")
cat("  Edges:", ecount(g), "\n")
cat("  Density:", round(edge_density(g), 4), "\n")
cat("  Mean degree:", round(mean(V(g)$node_degree), 2), "\n")
cat("  Components:", components(g)$no, "\n\n")

# ── Node metrics table ──────────────────────────────────────────────────────

# Build a lowest-available taxonomy label for each OTU
tax_priority <- tax_lookup %>%
  select(OTU_Label, Species, Genus, Family, Order, Class) %>%
  mutate(
    lowest_tax = case_when(
      !is.na(Species) & Species != "" ~ paste0(Species, " (species)"),
      !is.na(Genus)   & Genus   != "" ~ paste0(Genus,   " (genus)"),
      !is.na(Family)  & Family  != "" ~ paste0(Family,  " (family)"),
      !is.na(Order)   & Order   != "" ~ paste0(Order,   " (order)"),
      !is.na(Class)   & Class   != "" ~ paste0(Class,   " (class)"),
      TRUE                            ~ "Unclassified"
    )
  ) %>%
  select(OTU_Label, lowest_tax, Class)

node_metrics <- as_data_frame(g, what = "vertices") %>%
  rename(OTU_Label = name) %>%
  select(OTU_Label, Class_plot, Division_plot, Cluster,
          node_degree, node_strength, node_betweenness,
         pos_degree, neg_degree) %>%
  left_join(tax_priority, by = "OTU_Label") %>%
  mutate(
    signed_dominance = (pos_degree - neg_degree) / (pos_degree + neg_degree),
    across(c(node_strength, node_betweenness), ~ round(.x, 4))
  ) %>%
  select(OTU_Label, lowest_tax, Class, Cluster,
          node_degree, node_strength, node_betweenness,
         pos_degree, neg_degree, signed_dominance) %>%
  arrange(desc(node_betweenness))

cat("Top 20 OTUs by betweenness centrality:\n")
print(node_metrics %>% head(50) %>%
        select(OTU_Label, lowest_tax, Class, Cluster,
               node_degree, node_strength, node_betweenness,
               pos_degree, neg_degree))

write.csv(node_metrics, file.path(TABLE_DIR, "S3C_node_metrics.csv"), row.names = FALSE)
write.csv(as_data_frame(g, what = "edges"),
          file.path(TABLE_DIR, "S3C_edge_list.csv"), row.names = FALSE)


# ==============================================================================
# 3-D: Important TAXA IDENTIFICATION
#
# Hub taxa are identified using two complementary criteria:
#
#   (1) High betweenness: OTUs that bridge otherwise disconnected parts of
#       the network. Removing these would fragment the network most.
#       Ecologically: potential keystones or "ecological engineers".
#
#   (2) High strength with predominantly negative edges: OTUs that strongly
#       and consistently avoid many other OTUs. This pattern is consistent
#       with competitive dominants.
#
#   (3) High strength with predominantly positive edges: OTUs that strongly
#       co-occur with many others. Consistent with facilitative taxa
#       (e.g., biofilm matrix producers that support other species).
# ==============================================================================

cat("--- 3-D: Important taxa ---\n\n")

# Betweenness: taxa above this sit on a disproportionate share of shortest paths
# Strength: taxa above this have consistently strong residual associations

cat("Node metric distributions (to guide threshold choice):\n")
print(quantile(node_metrics$node_betweenness, probs = c(0.5, 0.75, 0.9, 0.95, 0.99)))
print(quantile(node_metrics$node_strength,    probs = c(0.5, 0.75, 0.9, 0.95, 0.99)))
cat("\n")

BETWEEN_THRESH <- 0.05   # ~95th percentile — only taxa that genuinely bridge network sections
STRENGTH_THRESH <- 42    # ~90th percentile — top 10% of connectors by association strength


metrics <- node_metrics %>%
  mutate(
    high_betweenness = node_betweenness > BETWEEN_THRESH,
    high_strength    = node_strength    > STRENGTH_THRESH
  )

cat("Taxa with high betweenness (>", BETWEEN_THRESH, "):\n")
metrics %>% filter(high_betweenness) %>%
  arrange(desc(node_betweenness)) %>%
  select(OTU_Label, Class, lowest_tax, node_betweenness, node_strength,
         node_degree, pos_degree, neg_degree) %>%
  print()
cat("\n")

cat("Taxa with high strength (>", STRENGTH_THRESH, "):\n")
metrics %>% filter(high_strength) %>%
  arrange(desc(node_strength)) %>%
  select(OTU_Label, Class, node_betweenness, node_strength,
         node_degree, pos_degree, neg_degree) %>%
  print()
cat("\n")



cat("All nodes with predominantly negative associations (>50% negative edges):\n")
node_metrics %>%
  filter((pos_degree + neg_degree) > 0) %>%
  mutate(
    pct_neg = neg_degree / (pos_degree + neg_degree),
    neg_dominated = pct_neg > 0.5
  ) %>%
  filter(neg_dominated) %>%
  arrange(desc(neg_degree), desc(pct_neg)) %>%
  # Join full taxonomy
  # left_join(
  #   model_data$tax_mat[otu_ids, tax_cols_use, drop = FALSE] %>%
  #     as.data.frame() %>%
  #     rownames_to_column("OTU_ID") %>%
  #     mutate(OTU_Label = model_data$otu_id_to_label[OTU_ID]),
  #   by = "OTU_Label"
  # ) %>%
  # select(OTU_Label,
  #        Division, Class = Class.x, Order, Family, Genus, Species,
  #        node_degree, pos_degree, neg_degree, pct_neg,
  #        node_betweenness, node_strength) %>%
  print()
cat("\n")

# Assocaitions summary

edges %>%
  mutate(
    class_pair = paste(
      pmin(Class_j, Class_k),
      pmax(Class_j, Class_k),
      sep = " × "
    )
  ) %>%
  count(class_pair, association, sort = TRUE) %>%
  print()




## FULL TABLE
# ── Combined important taxa table ─────────────────────────────────────────────

# 1. High betweenness taxa
high_between_tbl <- metrics %>%
  filter(high_betweenness) %>%
  arrange(desc(node_betweenness)) %>%
  mutate(category = "High betweenness")

# 2. High strength taxa (not already in betweenness set)
high_strength_tbl <- metrics %>%
  filter(high_strength & !high_betweenness) %>%
  arrange(desc(node_strength)) %>%
  mutate(category = "High strength")

# 3. Negative-dominated taxa
neg_dom_tbl <- node_metrics %>%
  filter((pos_degree + neg_degree) > 0) %>%
  mutate(
    pct_neg      = neg_degree / (pos_degree + neg_degree),
    neg_dominated = pct_neg > 0.5,
    category     = "Negative-dominated"
  ) %>%
  filter(neg_dominated) %>%
  arrange(desc(neg_degree))

# ── Bind and format ───────────────────────────────────────────────────────────

important_taxa_table <- bind_rows(
  high_between_tbl,
  high_strength_tbl,
  neg_dom_tbl
) %>%
  mutate(
    pct_pos = round(pos_degree / (pos_degree + neg_degree) * 100, 1),
    pct_neg = round(neg_degree / (pos_degree + neg_degree) * 100, 1)
  ) %>%
  select(
    Category         = category,
    Taxon            = lowest_tax,
    Class,
    Betweenness      = node_betweenness,
    Degree           = node_degree,
    `Node Strength`   = node_strength,
    `Positive edges (%)` = pct_pos,
    `Negative edges (%)` = pct_neg
  )


# ── For the network plot: label taxa meeting either threshold ─────────────────

hub_labels <- metrics %>%
  filter(high_betweenness | high_strength) %>%
  arrange(desc(node_betweenness), desc(node_strength)) %>%
  mutate(
    hub_number = row_number(),
    flags = case_when(
      high_betweenness & high_strength ~ "betweenness + strength",
      high_betweenness                 ~ "betweenness",
      high_strength                    ~ "strength"
    )
  ) %>%
  select(OTU_Label, flags, hub_number,
         node_betweenness, node_strength, node_degree,
         Class, pos_degree, neg_degree)

cat("Labelled taxa in network plot:\n")
print(hub_labels)
cat("\n")

tg <- as_tbl_graph(g) %>%
  left_join(
    metrics %>%
      select(name = OTU_Label, high_betweenness, high_strength,
             node_betweenness, node_strength),
    by = "name"
  ) %>%
  left_join(
    hub_labels %>% select(name = OTU_Label, hub_number),
    by = "name"
  ) %>%
  mutate(
    node_label = if_else(!is.na(hub_number),
                         as.character(hub_number),
                         ""),
    is_hub     = !is.na(hub_number)
  )


# ── Build class colour palette ────────────────────────────────────────────────
# Get all unique classes present in the network nodes
all_classes_net <- unique(as_data_frame(tg, what = "vertices")$Class_plot)
all_classes_net <- all_classes_net[!is.na(all_classes_net)]

n_cls <- length(all_classes_net)

# Generate enough colours — Paired + Set1 + Set3 combined handles up to ~36
class_cols <- setNames(
  colorRampPalette(
    c(brewer.pal(12, "Paired"),
      brewer.pal(9,  "Set1"),
      brewer.pal(12, "Set3"))
  )(n_cls),
  all_classes_net
)



set.seed(42)
p_net_numbered <- ggraph(tg, layout = "stress") +
  
  # ── Edges ────────────────────────────────────────────────────────────────────
  geom_edge_link(
    aes(colour = association,
        width  = weight,
        alpha  = weight),
    show.legend = TRUE
  ) +
  
  # ── Non-hub nodes: small filled circles ──────────────────────────────────────
  geom_node_point(
    data  = . %>% filter(!is_hub),
    aes(colour = Class),
    size  = 5,
    alpha = 0.75
  ) +
  
  # ── Hub nodes: coloured outer ring ───────────────────────────────────────────
  geom_node_point(
    data  = . %>% filter(is_hub),
    aes(colour = Class),
    size  = 14,
    alpha = 1
  ) +
  
  # ── Hub nodes: white inner fill ──────────────────────────────────────────────
  geom_node_point(
    data   = . %>% filter(is_hub),
    colour = "white",
    size   = 10,
    alpha  = 1
  ) +
  
  # ── Hub number labels ─────────────────────────────────────────────────────────
  geom_node_text(
    data     = . %>% filter(is_hub),
    aes(label  = node_label,
        colour = Class),
    size     = 4,
    fontface = "bold"
  ) +
  
  # ── Scales ────────────────────────────────────────────────────────────────────
  scale_edge_colour_manual(
    values = c(positive = "#4575b4", negative = "#d73027"),
    name   = "Association",
    labels = c(positive = "Positive (co-occurrence)",
               negative = "Negative (avoidance)")
  ) +
  scale_edge_width(range = c(0.2, 2.0), guide = "none") +
  scale_edge_alpha(range = c(0.2, 0.8), guide = "none") +
  scale_colour_manual(
    values = class_cols,
    name   = "Taxonomic class",
    guide  = guide_legend(
      override.aes = list(shape = 16, size = 4, alpha = 1)
    )
  ) +
  
  # labs(
  #   title    = "Residual co-occurrence network",
  #   subtitle = paste0(
  #     "Edges: |residual r| > ", R_STRONG,
  #     " after removing site + season + substrate effects\n",
  #     "Numbered circles = hub taxa (Bridge or High strength) | ",
  #     "Node colour = taxonomic class\n",
  #     "Blue edges = positive co-occurrence | Red edges = avoidance"
  #   )
  # ) +
  theme_graph(base_family = "sans", base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size  = 9),
    legend.position = "right"
  )

print(p_net_numbered)

ggsave(
  file.path(FIG_DIR, "S3C_network_numbered_stress.png"),
  p_net_numbered,
  width  = 14,
  height = 12,
  dpi    = 200
)
cat("Numbered network plot (stress layout) saved.\n\n")


##TABLE

library(flextable)
library(officer)   # needed to write to .docx

# ── Build the table data ──────────────────────────────────────────────────────

important_taxa_table <- bind_rows(
  high_between_tbl,
  high_strength_tbl,
  neg_dom_tbl
) %>%
  mutate(
    pct_pos = round(pos_degree / (pos_degree + neg_degree) * 100, 1),
    pct_neg = round(neg_degree / (pos_degree + neg_degree) * 100, 1)
  ) %>%
  select(
    Category             = category,
    Taxon                = lowest_tax,
    Class,
    Betweenness          = node_betweenness,
    Degree               = node_degree,
    `Node Strength`      = node_strength,
    `Positive edges (%)` = pct_pos,
    `Negative edges (%)` = pct_neg
  )

# ── Build the flextable ───────────────────────────────────────────────────────

ft <- flextable(important_taxa_table) %>%
  
  # Merge repeated Category values in the first column
  merge_v(j = "Category") %>%
  valign(j = "Category", valign = "top") %>%
  
  # Header styling
  bold(part = "header") %>%
  bg(part = "header", bg = "#2C3E50") %>%
  color(part = "header", color = "white") %>%
  align(part = "header", align = "center") %>%
  
  # Zebra-stripe body rows
  bg(i = ~ as.integer(rownames(important_taxa_table)) %% 2 == 0,
     bg = "#F2F2F2", part = "body") %>%
  
  # Highlight negative-dominated rows in a soft red
  bg(i = ~ Category == "Negative-dominated",
     bg = "#FDECEA", part = "body") %>%
  
  # Round numeric columns to 3 decimal places for Betweenness & Node Strength
  colformat_double(
    j      = c("Betweenness", "Node Strength"),
    digits = 3
  ) %>%
  colformat_double(
    j      = c("Positive edges (%)", "Negative edges (%)"),
    digits = 1
  ) %>%
  
  # Column widths (inches) — adjust to taste
  width(j = "Category",            width = 1.4) %>%
  width(j = "Taxon",               width = 1.6) %>%
  width(j = "Class",               width = 1.4) %>%
  width(j = "Betweenness",         width = 1.1) %>%
  width(j = "Degree",              width = 0.8) %>%
  width(j = "Node Strength",       width = 1.1) %>%
  width(j = "Positive edges (%)",  width = 1.2) %>%
  width(j = "Negative edges (%)",  width = 1.2) %>%
  
  # General alignment
  align(j = c("Betweenness", "Degree", "Node Strength",
              "Positive edges (%)", "Negative edges (%)"),
        align = "right", part = "body") %>%
  align(j = c("Category", "Taxon", "Class"),
        align = "left", part = "body") %>%
  
  # Font
  font(fontname = "Calibri", part = "all") %>%
  flextable::fontsize(size = 10, part = "all") %>%
  flextable::fontsize(size = 11, part = "header") %>%
  
  # Borders
  border_outer(part = "all",
               border = officer::fp_border(color = "#2C3E50", width = 1.5)) %>%
  border_inner_h(part = "body",
                 border = officer::fp_border(color = "#CCCCCC", width = 0.5)) %>%
  
  # Auto-fit row height
  autofit()

# ── Write to Word ─────────────────────────────────────────────────────────────

doc <- read_docx() %>%
  body_add_par("Important Taxa — Network Summary", style = "heading 1") %>%
  body_add_flextable(ft) %>%
  body_end_section_landscape()

print(doc, target = file.path(TABLE_DIR, "important_taxa_table.docx"))
cat("Table saved to Word document.\n")




# ==============================================================================
# 3-D (continued): TRANSITION SUMMARY — WHY SYNDINIALES × DIATOMS?
#
# Before zooming into the synd-diatom sub-network, we summarise which classes
# dominate the hub taxa. This makes the motivation for the targeted analysis
# explicit rather than implied.
# ==============================================================================

cat("--- 3-D: Hub taxa class summary (transition to sub-network) ---\n\n")

hub_class_summary <- hub_labels %>%
  group_by(Class) %>%
  summarise(
    n_hubs        = n(),
    hub_numbers   = paste(sort(hub_number), collapse = ", "),
    mean_between  = round(mean(node_betweenness), 4),
    mean_strength = round(mean(node_strength),    3),
    .groups = "drop"
  ) %>%
  arrange(desc(n_hubs))

cat("Classes represented among hub taxa:\n")
print(hub_class_summary, n = Inf)

# How many of the top hubs are Syndiniales or diatoms?
synd_diatom_hubs <- hub_labels %>%
  filter(Class %in% c("Syndiniales", "Bacillariophyceae", "Mediophyceae"))

cat(sprintf(
  "\n%d of %d hub taxa are Syndiniales or diatoms (Bacillariophyceae / Mediophyceae).\n",
  nrow(synd_diatom_hubs), nrow(hub_labels)
))
cat("These classes also dominate the Pollock circos plots.\n")
cat("Proceeding to ASV-level sub-network analysis.\n\n")


# ==============================================================================
# 3-D: SYNDINIALES–DIATOM SUB-NETWORK ANALYSIS
#
# The Pollock decomposition (3-A) identified Syndiniales and diatoms as the
# dominant class pair in the top surprise co-occurrences. Hub analysis (above)
# confirms that Syndiniales and diatom ASVs account for the majority of high-
# betweenness and high-strength nodes in the full community network.
#
# We now zoom to ASV level, asking:
#   (1) Which specific Syndiniales and diatom ASVs are associated?
#   (2) What is the Pollock quadrant of each pair — does the association
#       exceed environmental expectation (Q3) or reflect shared habitat (Q1)?
#   (3) Is the association consistent across sites, seasons, and substrates?
# ==============================================================================

cat("=== SECTION 3-D: SYNDINIALES–DIATOM SUB-NETWORK ===\n\n")


# ------------------------------------------------------------------------------
# STEP 1: IDENTIFY ALL SYNDINIALES–DIATOM PAIRS FROM THE FULL NETWORK
#
# We filter the edge list (already built in 3-C) to pairs where one partner
# is Syndiniales and the other is Bacillariophyceae or Mediophyceae.
# Both positive (co-occurrence) and negative (avoidance) edges are kept.
# ------------------------------------------------------------------------------

cat("--- Step 1: Syndiniales–diatom pairs from network edges ---\n\n")

synd_diatom_edges <- edges %>%
  left_join(
    node_metrics %>% select(OTU_Label, Class_from = Class),
    by = c("OTU_j" = "OTU_Label")
  ) %>%
  left_join(
    node_metrics %>% select(OTU_Label, Class_to = Class),
    by = c("OTU_k" = "OTU_Label")
  ) %>%
  filter(
    (Class_from == "Syndiniales" & Class_to %in% c("Bacillariophyceae", "Mediophyceae")) |
      (Class_from %in% c("Bacillariophyceae", "Mediophyceae") & Class_to == "Syndiniales")
  ) %>%
  mutate(
    synd_otu   = if_else(Class_from == "Syndiniales", OTU_j, OTU_k),
    diatom_otu = if_else(Class_from == "Syndiniales", OTU_k, OTU_j),
    diatom_class = if_else(Class_from == "Syndiniales", Class_to, Class_from)
  )

cat("Syndiniales–diatom pairs in network (|r| >", R_STRONG, "):", nrow(synd_diatom_edges), "\n")
cat("  Positive (co-occurrence):", sum(synd_diatom_edges$association == "positive"), "\n")
cat("  Negative (avoidance):",     sum(synd_diatom_edges$association == "negative"), "\n")
cat("  Unique Syndiniales ASVs:", n_distinct(synd_diatom_edges$synd_otu), "\n")
cat("  Unique diatom ASVs:     ", n_distinct(synd_diatom_edges$diatom_otu), "\n\n")

if (nrow(synd_diatom_edges) == 0) {
  cat("WARNING: No Syndiniales–diatom pairs found. Check Class labels in node_metrics.\n\n")
  stop("Stopping 3-D sub-network: no pairs to analyse.")
}

# Add Pollock quadrant — already present in edges via pair_df join in 3-C
synd_diatom_pollock <- synd_diatom_edges %>%
  mutate(pollock_signal = quadrant)

cat("Pollock quadrant breakdown:\n")
print(table(synd_diatom_pollock$association,
            synd_diatom_pollock$pollock_signal,
            useNA = "ifany"))
cat("\n")

write.csv(synd_diatom_pollock,
          file.path(TABLE_DIR, "S3D_synd_diatom_pairs.csv"),
          row.names = FALSE)


# ------------------------------------------------------------------------------
# STEP 2: FLAG HUB TAXA WITHIN THE SUB-NETWORK
#
# Cross-reference synd and diatom ASVs against hub_labels so the PA matrix
# and sub-network plot can mark which ASVs are also community-level hubs.
# ------------------------------------------------------------------------------

hub_otus <- hub_labels$OTU_Label

synd_diatom_pollock <- synd_diatom_pollock %>%
  mutate(
    synd_is_hub    = synd_otu   %in% hub_otus,
    diatom_is_hub  = diatom_otu %in% hub_otus,
    synd_hub_num   = hub_labels$hub_number[match(synd_otu,   hub_labels$OTU_Label)],
    diatom_hub_num = hub_labels$hub_number[match(diatom_otu, hub_labels$OTU_Label)]
  )

cat("Hub taxa appearing in synd-diatom pairs:\n")
synd_diatom_pollock %>%
  filter(synd_is_hub | diatom_is_hub) %>%
  select(synd_otu, synd_hub_num, diatom_otu, diatom_hub_num,
         association, res_cor, quadrant) %>%
  print()
cat("\n")


# ------------------------------------------------------------------------------
# STEP 3: BUILD READABLE LABELS FOR PLOTTING
# ------------------------------------------------------------------------------

clean_label <- function(x) {
  x %>%
    gsub("Dino-Group-I-Clade-",   "DGI C",   .) %>%
    gsub("Dino-Group-II-Clade-",  "DGII C",  .) %>%
    gsub("Dino-Group-III-Clade-", "DGIII C", .) %>%
    gsub("_X_ASV_([0-9]+)",       " (\\1)",  .) %>%
    gsub("_ASV_([0-9]+)$",        " (\\1)",  .)
}

# Diatom labels: lowest taxonomy + ASV number for uniqueness
diatom_labels <- node_metrics %>%
  filter(OTU_Label %in% synd_diatom_pollock$diatom_otu) %>%
  select(OTU_Label, lowest_tax, Class) %>%
  left_join(
    tax_lookup %>% select(OTU_Label, Species, Genus),
    by = "OTU_Label"
  ) %>%
  mutate(
    diatom_plot_label = case_when(
      !is.na(Species) & Species != "" ~ Species,
      !is.na(Genus)   & Genus   != "" ~ paste0(Genus, " sp."),
      TRUE                            ~ Class
    ),
    diatom_plot_label = paste0(
      diatom_plot_label, " (", sub(".*_ASV_", "", OTU_Label), ")"
    )
  ) %>%
  select(OTU_Label, diatom_plot_label)

pa_data <- synd_diatom_pollock %>%
  mutate(synd_label = clean_label(synd_otu)) %>%
  left_join(diatom_labels, by = c("diatom_otu" = "OTU_Label")) %>%
  rename(diatom_label = diatom_plot_label)


# ------------------------------------------------------------------------------
# STEP 4: PRESENCE/ABSENCE MATRIX — SYNDINIALES × DIATOM ASV PAIRS
#
# Rows = diatom ASVs, columns = Syndiniales ASVs.
# Filled cell = strong residual association exists (|r| > R_STRONG).
# Colour = Pollock quadrant.
# Bold axis label = community-level hub taxon.
# ------------------------------------------------------------------------------

cat("--- Step 4: Presence/absence matrix ---\n\n")

# Order axes by connectivity
synd_order <- pa_data %>%
  count(synd_label, sort = TRUE) %>%
  pull(synd_label)

diatom_order <- pa_data %>%
  group_by(diatom_label) %>%
  summarise(
    n_synd = n_distinct(synd_label),
    mean_r = mean(res_cor),
    .groups = "drop"
  ) %>%
  arrange(desc(n_synd), desc(mean_r)) %>%
  pull(diatom_label)

# Build hub lookup BEFORE factor conversion (columns still exist here)
synd_hub_lookup <- pa_data %>%
  distinct(synd_label, synd_is_hub)

diatom_hub_lookup <- pa_data %>%
  distinct(diatom_label, diatom_is_hub)

# Convert to factors
pa_data <- pa_data %>%
  mutate(
    synd_label   = factor(synd_label,   levels = rev(synd_order)),
    diatom_label = factor(diatom_label, levels = diatom_order),
    signal_type  = factor(pollock_signal)
  )

# Build face vectors aligned to factor levels
# Must match the order of levels() exactly for ggplot axis text
synd_face <- ifelse(
  levels(pa_data$synd_label) %in%
    synd_hub_lookup$synd_label[synd_hub_lookup$synd_is_hub],
  "bold", "plain"
)

diatom_face <- ifelse(
  levels(pa_data$diatom_label) %in%
    diatom_hub_lookup$diatom_label[diatom_hub_lookup$diatom_is_hub],
  "bold", "plain"
)

signal_cols <- c(
  "Q1" = "#2166AC",
  "Q2" = "#F4A582",
  "Q3" = "#2AB5A0",
  "Q4" = "#B2182B"
)

n_synd   <- n_distinct(pa_data$synd_label)
n_diatom <- n_distinct(pa_data$diatom_label)

p_pa <- ggplot(pa_data,
               aes(x = synd_label, y = diatom_label, fill = signal_type)) +
  
  geom_tile(colour = "white", linewidth = 0.4) +
  
  scale_fill_manual(
    values = signal_cols,
    name   = "Association type",
    labels = c(
      Q1 = "Q1",
      Q2 = "Q2",
      Q3 = "Q3",
      Q4 = "Q4"
    ),
    drop = TRUE
  ) +
  scale_x_discrete(position = "top") +
  labs(
    x       = "Syndiniales ASV",
    y       = "Diatom ASV",
    caption = paste0(
      "Filled cell = strong residual association |r| > ", R_STRONG,
      " | Bold labels = community-level hub taxon"
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x      = element_text(
      angle = 45,
      hjust = 0,
      size  = 9,
      face  = synd_face
    ),
    axis.text.y      = element_text(
      size  = 9,
      face  = diatom_face
    ),
    axis.title.x     = element_text(size = 10, face = "italic"),
    axis.title.y     = element_text(size = 10, face = "italic"),
    panel.grid       = element_blank(),
    legend.position  = "bottom",
    legend.title     = element_text(size = 10),
    legend.text      = element_text(size = 9),
    plot.caption     = element_text(size = 8, colour = "grey40"),
    plot.margin      = margin(5, 10, 5, 5)
  )

fig_w <- max(7,  n_synd   * 1.1 + 3)
fig_h <- max(6,  n_diatom * 0.45 + 3)

ggsave(
  file.path(FIG_DIR, "S3D_synd_diatom_PA_matrix.png"),
  p_pa,
  width  = fig_w,
  height = fig_h,
  dpi    = 200
)
print(p_pa)
cat("PA matrix saved: S3D_synd_diatom_PA_matrix.png\n\n")

# ------------------------------------------------------------------------------
# STEP 5: SEASONAL DETECTION PLOTS
#
# Are Syndiniales and diatoms detected in the same seasons?
# Point size = number of ASVs detected at least once in that group × condition.
# This is a raw PA check — not GLLVM-corrected — but gives biological context
# for interpreting the residual associations above.
# ------------------------------------------------------------------------------

cat("--- Step 5: Seasonal detection frequency ---\n\n")

# All Syndiniales and diatom ASVs in the full model (not just those with edges)
all_synd_otus <- node_metrics %>%
  filter(Class == "Syndiniales") %>%
  pull(OTU_Label)

all_diatom_otus <- node_metrics %>%
  filter(Class %in% c("Bacillariophyceae", "Mediophyceae")) %>%
  pull(OTU_Label)

cat("Total Syndiniales ASVs in network:", length(all_synd_otus), "\n")
cat("Total diatom ASVs in network:     ", length(all_diatom_otus), "\n\n")

# Factor level ordering — adjust to match your metadata
seasons          <- levels(meta_diag$season)
substrate_levels <- levels(meta_diag$substrate)

# Helper: mean detection frequency + n ASVs detected per condition group
compute_freq_full <- function(otu_list, group_vars) {
  valid <- otu_list[otu_list %in% colnames(pa_full)]
  
  meta_diag %>%
    select(all_of(group_vars)) %>%
    distinct() %>%
    rowwise() %>%
    mutate(
      idx = list({
        mask <- rep(TRUE, nrow(meta_diag))
        for (gv in group_vars) mask <- mask & (meta_diag[[gv]] == get(gv))
        which(mask)
      }),
      n_samples       = length(idx),
      freq_pct        = if (n_samples < 2 | length(valid) == 0) NA_real_
      else round(mean(pa_full[idx, valid]) * 100, 1),
      n_asvs_detected = if (n_samples < 2 | length(valid) == 0) NA_integer_
      else sum(colSums(pa_full[idx, valid]) > 0L)
    ) %>%
    select(-idx) %>%
    ungroup()
}

synd_season   <- compute_freq_full(all_synd_otus,   "season") %>%
  mutate(taxon = "Syndiniales",       n_asvs_total = length(all_synd_otus))
diatom_season <- compute_freq_full(all_diatom_otus, "season") %>%
  mutate(taxon = "Bacillariophyceae / Mediophyceae", n_asvs_total = length(all_diatom_otus))

synd_sub   <- compute_freq_full(all_synd_otus,   c("season", "substrate")) %>%
  mutate(taxon = "Syndiniales",       n_asvs_total = length(all_synd_otus))
diatom_sub <- compute_freq_full(all_diatom_otus, c("season", "substrate")) %>%
  mutate(taxon = "Bacillariophyceae / Mediophyceae", n_asvs_total = length(all_diatom_otus))

class_season <- bind_rows(synd_season, diatom_season) %>%
  mutate(season = factor(season, levels = seasons))

class_sub <- bind_rows(synd_sub, diatom_sub) %>%
  mutate(
    season    = factor(season,    levels = seasons),
    substrate = factor(substrate, levels = substrate_levels)
  )

taxon_cols <- c(
  "Syndiniales"                        = "#7B2D8B",
  "Bacillariophyceae / Mediophyceae"   = "#2E8B57"
)

shared_scales <- list(
  scale_colour_manual(values = taxon_cols, name = NULL),
  scale_size_continuous(
    name  = "ASVs detected in group",
    range = c(2, 8),
    guide = guide_legend(override.aes = list(colour = "grey40"))
  ),
  scale_y_continuous(
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100),
    labels = function(x) paste0(x, "%")
  ),
  labs(
    x       = NULL,
    y       = "Mean detection frequency (% of samples)",
    caption = "Point size = number of ASVs detected at least once in that group"
  ),
  theme_bw(base_size = 12),
  theme(
    strip.text       = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "grey92", colour = NA),
    axis.text.x      = element_text(angle = 30, hjust = 1, size = 9),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )
)

# Plot 1: all substrates collapsed
p_season <- ggplot(
  class_season,
  aes(x = season, y = freq_pct,
      colour = taxon, group = taxon,
      size   = n_asvs_detected)
) +
  geom_line(linewidth = 1, alpha = 0.85, na.rm = TRUE) +
  geom_point(alpha = 0.9, na.rm = TRUE) +
  shared_scales

ggsave(file.path(FIG_DIR, "S3D_seasonal_detection_class.png"),
       p_season, width = 6, height = 5, dpi = 200)
print(p_season)
cat("Figure saved: S3D_seasonal_detection_class.png\n\n")

# Plot 2: faceted by substrate
p_substrate <- ggplot(
  class_sub,
  aes(x = season, y = freq_pct,
      colour = taxon, group = taxon,
      size   = n_asvs_detected)
) +
  geom_line(linewidth = 1, alpha = 0.85, na.rm = TRUE) +
  geom_point(alpha = 0.9, na.rm = TRUE) +
  facet_wrap(~ substrate, ncol = 2) +
  shared_scales

ggsave(
  file.path(FIG_DIR, "S3D_seasonal_detection_by_substrate.png"),
  p_substrate,
  width  = 8,
  height = ceiling(length(substrate_levels) / 2) * 3.5,
  dpi    = 200
)
print(p_substrate)
cat("Figure saved: S3D_seasonal_detection_by_substrate.png\n\n")


# ------------------------------------------------------------------------------
# STEP 6: SAVE ANNOTATED PAIR TABLE TO WORD
# ------------------------------------------------------------------------------

cat("--- Step 6: Annotated pair table ---\n\n")

sd_table <- synd_diatom_pollock %>%
  left_join(
    node_metrics %>% select(OTU_Label, synd_tax = lowest_tax),
    by = c("synd_otu" = "OTU_Label")
  ) %>%
  left_join(diatom_labels, by = c("diatom_otu" = "OTU_Label")) %>%
  mutate(
    `Syndiniales ASV` = if_else(
      synd_is_hub,
      paste0(synd_tax, " [hub ", synd_hub_num, "]"),
      synd_tax
    ),
    `Diatom ASV` = if_else(
      diatom_is_hub,
      paste0(diatom_plot_label, " [hub ", diatom_hub_num, "]"),
      diatom_plot_label
    )
  ) %>%
  select(
    `Syndiniales ASV`,
    `Diatom ASV`,
    `Diatom class`   = diatom_class,
    `Residual r`     = res_cor,
    `Env. r`         = env_cor,
    Association      = association,
    Quadrant         = pollock_signal,
    `Synd. hub`      = synd_is_hub,
    `Diatom hub`     = diatom_is_hub
  ) %>%
  arrange(Association, desc(abs(`Residual r`)))

quadrant_fill <- c(
  "Q1" = "#CCE5FF",
  "Q2" = "#FFE0B2",
  "Q3" = "#D5F5E3",
  "Q4" = "#FFCCCC"
)

ft_sd <- flextable(sd_table) %>%
  theme_vanilla() %>%
  bold(j = c("Syndiniales ASV", "Diatom ASV")) %>%
  color(
    i = ~ `Synd. hub` == TRUE | `Diatom hub` == TRUE,
    j = c("Syndiniales ASV", "Diatom ASV"),
    color = "#8B0000"
  ) %>%
  align(
    j     = c("Residual r", "Env. r"),
    align = "right", part = "all"
  ) %>%
  colformat_double(j = c("Residual r", "Env. r"), digits = 3) %>%
  bg(part = "header", bg = "#2C3E50") %>%
  color(part = "header", color = "white") %>%
  bold(part = "header")

# Colour rows by quadrant
for (q in names(quadrant_fill)) {
  rows <- which(sd_table$Quadrant == q)
  if (length(rows) == 0) next
  ft_sd <- ft_sd %>% bg(i = rows, bg = quadrant_fill[q], part = "body")
}

ft_sd <- ft_sd %>%
  width(j = "Syndiniales ASV", width = 2.0) %>%
  width(j = "Diatom ASV",      width = 2.0) %>%
  width(j = "Diatom class",    width = 1.3) %>%
  width(j = "Residual r",      width = 0.9) %>%
  width(j = "Env. r",          width = 0.9) %>%
  width(j = "Association",     width = 1.0) %>%
  width(j = "Quadrant",        width = 0.8) %>%
  width(j = "Synd. hub",       width = 0.8) %>%
  width(j = "Diatom hub",      width = 0.9) %>%
  flextable::fontsize(size = 9,  part = "body") %>%
  flextable::fontsize(size = 10, part = "header") %>%
  font(fontname = "Calibri", part = "all") %>%
  border_outer(border = officer::fp_border(color = "#2C3E50", width = 1.5)) %>%
  border_inner_h(border = officer::fp_border(color = "#CCCCCC", width = 0.5))

doc_sd <- read_docx() %>%
  body_add_par("Syndiniales–Diatom Residual Associations", style = "heading 1") %>%
  body_add_par(paste0(
    "Strong residual pairs (|r| > ", R_STRONG, ") between Syndiniales and diatom ASVs. ",
    "Row colour = Pollock quadrant. Hub taxa (cross-referenced to network figure) ",
    "shown in dark red with hub number in brackets."
  ), style = "Normal") %>%
  body_add_flextable(ft_sd) %>%
  body_end_section_landscape()

print(doc_sd, target = file.path(TABLE_DIR, "S3D_synd_diatom_table.docx"))
cat("Syndiniales–diatom table saved to Word.\n\n")

cat("=== SECTION 3-D COMPLETE ===\n\n")

# ==============================================================================
# 3-D SUPPLEMENT: RELATIVE ABUNDANCE & FOO — SYNDINIALES vs BACILLARIOPHYCEAE
# Full dataset (ps_combined), Version A (denominator = all reads in sample)
#
# FIGURE 1 — Boxplot: RA distribution, Water vs Biofilm, jittered dots
# FIGURE 2 — Seasonal barplot: mean RA, facet grid taxon × sample_type
# FIGURE 3 — Substrate barplot: mean RA across substrate types
# ==============================================================================

library(phyloseq)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

# ── Load data ──────────────────────────────────────────────────────────────────
PS_DIR <- "~/Github/MAPLE_Seasonal_Plastisphere/Processed_data/phyloseq_objects"
ps_combined <- readRDS(file.path(PS_DIR, "ps_st_wf_allST_combined.rds"))

sample_data(ps_combined)$sample_type <- factor(
  sample_data(ps_combined)$sample_type,
  levels = c("Water", "Biofilm")
)
sample_data(ps_combined)$season <- factor(
  sample_data(ps_combined)$season,
  levels = c("Winter", "Spring", "Summer", "Fall", "Winter2")
)
sample_data(ps_combined)$substrate <- factor(
  sample_data(ps_combined)$substrate,
  levels = c("Filter", "Glass", "PE", "Weathered_PE", "PET", "Weathered_PET")
)

ps_raw               <- ps_combined
ps_synd_raw          <- subset_taxa(ps_raw, Class == "Syndiniales")
ps_bacill_raw        <- subset_taxa(ps_raw, Class == "Bacillariophyceae")

cat("Syndiniales ASVs:        ", ntaxa(ps_synd_raw),   "\n")
cat("Bacillariophyceae ASVs:  ", ntaxa(ps_bacill_raw), "\n\n")

# ── Colour palettes ────────────────────────────────────────────────────────────
# Using consistent label "Bacillariophyceae" everywhere — no trailing s

LABEL_SYND   <- "Syndiniales"
LABEL_BACILL <- "Bacillariophyceae"
taxon_levels <- c(LABEL_SYND, LABEL_BACILL)

taxon_cols <- c(
  "Syndiniales"       = "#7B2D8B",
  "Bacillariophyceae" = "#2E8B57"
)
taxon_fill <- c(
  "Syndiniales"       = adjustcolor("#7B2D8B", alpha.f = 0.45),
  "Bacillariophyceae" = adjustcolor("#2E8B57", alpha.f = 0.45)
)

substrate_levels <- c("Filter", "Glass", "PE", "Weathered_PE", "PET", "Weathered_PET")
season_levels    <- c("Winter", "Spring", "Summer", "Fall", "Winter2")

# ── Helper: per-sample RA (Version A — denominator = all reads) ────────────────
extract_ra <- function(ps_subset, ps_total, group_label) {
  
  num_mat <- as(otu_table(ps_subset), "matrix")
  if (taxa_are_rows(ps_subset)) num_mat <- t(num_mat)
  
  den_mat <- as(otu_table(ps_total), "matrix")
  if (taxa_are_rows(ps_total)) den_mat <- t(den_mat)
  
  total_reads <- rowSums(den_mat)
  group_reads <- rowSums(num_mat)
  ra_pct      <- ifelse(total_reads > 0, group_reads / total_reads * 100, NA_real_)
  
  meta <- data.frame(sample_data(ps_subset), stringsAsFactors = FALSE)
  
  data.frame(
    sample_id   = rownames(meta),
    ra_pct      = ra_pct,
    sample_type = meta$sample_type,
    season      = meta$season,
    substrate   = meta$substrate,
    taxon_group = group_label,
    stringsAsFactors = FALSE
  )
}

# ── Build per-sample RA data frame ─────────────────────────────────────────────
ra_all <- bind_rows(
  extract_ra(ps_synd_raw,   ps_raw, LABEL_SYND),
  extract_ra(ps_bacill_raw, ps_raw, LABEL_BACILL)
) %>%
  mutate(
    sample_type = factor(sample_type, levels = c("Water", "Biofilm")),
    season      = factor(season,      levels = season_levels),
    substrate   = factor(substrate,   levels = substrate_levels),
    taxon_group = factor(taxon_group, levels = taxon_levels)
  )

cat("Total per-sample observations:", nrow(ra_all), "\n\n")

# ==============================================================================
# FIGURE 1: BOXPLOT — overall RA distribution, Water vs Biofilm
# ==============================================================================

p_box <- ggplot(
  ra_all %>% filter(!is.na(ra_pct)),
  aes(x      = taxon_group,
      y      = ra_pct,
      fill   = taxon_group,
      colour = taxon_group)
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha  = 0.45,
    width  = 0.5,
    colour = "grey30"
  ) +
  geom_jitter(
    width = 0.18,
    size  = 1.0,
    alpha = 0.40
  ) +
  facet_wrap(~ sample_type, ncol = 2) +
  scale_fill_manual(values   = taxon_fill, guide = "none") +
  scale_colour_manual(values = taxon_cols, guide = "none") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(
    x = NULL,
    y = "Relative abundance"
  ) +
  theme_classic(base_size = 14) +
  theme(
    strip.text.x       = element_text(size = 12, face = "bold"),
    strip.background   = element_rect(fill = "grey92", colour = NA),
    axis.title.y       = element_text(size = 13, face = "bold"),
    axis.text.x        = element_text(size = 12, face = "bold.italic"),
    axis.text.y        = element_text(size = 12, face = "bold"),
    axis.line          = element_line(linewidth = 0.7, colour = "black"),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.4),
    plot.margin        = margin(5, 10, 5, 5)
  )

print(p_box)
ggsave(
  file.path(FIG_DIR, "S3D_ra_boxplot.png"),
  p_box, width = 8, height = 5, dpi = 200
)
cat("Figure 1 saved: S3D_ra_boxplot.png\n\n")


# ==============================================================================
# FIGURE 2: SEASONAL BARPLOT — mean RA, facet grid taxon × sample_type
# ==============================================================================

season_summary <- ra_all %>%
  filter(!is.na(ra_pct), !is.na(season)) %>%
  group_by(taxon_group, season, sample_type) %>%
  summarise(
    mean_ra   = mean(ra_pct, na.rm = TRUE),
    se_ra     = sd(ra_pct,   na.rm = TRUE) / sqrt(n()),
    n_samples = n(),
    .groups   = "drop"
  )

cat("Seasonal RA summary:\n")
print(season_summary, n = Inf)
cat("\n")

p_ra_season <- ggplot(
  season_summary,
  aes(x = season, y = mean_ra, fill = taxon_group)
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width    = 0.65,
    alpha    = 0.85
  ) +
  geom_errorbar(
    aes(ymin = pmax(mean_ra - se_ra, 0),
        ymax = mean_ra + se_ra),
    position  = position_dodge(width = 0.75),
    width     = 0.25,
    colour    = "grey30",
    linewidth = 0.5
  ) +
  facet_grid(taxon_group ~ sample_type, scales = "free_y") +
  scale_fill_manual(values = taxon_cols, name = NULL) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(
    x = NULL,
    y = "Mean relative abundance"
  ) +
  theme_classic(base_size = 14) +
  theme(
    strip.text.x       = element_text(size = 12, face = "bold"),
    strip.text.y       = element_text(size = 11, face = "bold.italic"),
    strip.background   = element_rect(fill = "grey92", colour = NA),
    axis.title.y       = element_text(size = 13, face = "bold"),
    axis.text.x        = element_text(angle = 30, hjust = 1,
                                      size = 11, face = "bold"),
    axis.text.y        = element_text(size = 11, face = "bold"),
    axis.line          = element_line(linewidth = 0.7, colour = "black"),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.4),
    legend.position    = "none",
    plot.margin        = margin(5, 10, 5, 5)
  )

print(p_ra_season)
ggsave(
  file.path(FIG_DIR, "S3D_ra_seasonal_bar.png"),
  p_ra_season, width = 9, height = 7, dpi = 200
)
cat("Figure 2 saved: S3D_ra_seasonal_bar.png\n\n")


# ==============================================================================
# FIGURE 3: SUBSTRATE BARPLOT — mean RA across substrate types
# Filter = seawater, all others = biofilm substrates
# No sample_type facet — substrate already encodes habitat
# ==============================================================================

substrate_summary <- ra_all %>%
  filter(!is.na(ra_pct), !is.na(substrate)) %>%
  group_by(taxon_group, substrate) %>%
  summarise(
    mean_ra   = mean(ra_pct, na.rm = TRUE),
    se_ra     = sd(ra_pct,   na.rm = TRUE) / sqrt(n()),
    n_samples = n(),
    .groups   = "drop"
  )

cat("Substrate RA summary:\n")
print(substrate_summary, n = Inf)
cat("\n")

p_ra_sub <- ggplot(
  substrate_summary,
  aes(x = substrate, y = mean_ra, fill = taxon_group)
) +
  geom_col(
    position = position_dodge(width = 0.75),
    width    = 0.65,
    alpha    = 0.85
  ) +
  geom_errorbar(
    aes(ymin = pmax(mean_ra - se_ra, 0),
        ymax = mean_ra + se_ra),
    position  = position_dodge(width = 0.75),
    width     = 0.25,
    colour    = "grey30",
    linewidth = 0.5
  ) +
  scale_fill_manual(values = taxon_cols, name = NULL) +
  scale_x_discrete(labels = c(
    "Filter"        = "Seawater",
    "Glass"         = "Glass",
    "PE"            = "PE",
    "Weathered_PE"  = "Weathered\nPE",
    "PET"           = "PET",
    "Weathered_PET" = "Weathered\nPET"
  ))+
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(
    x = NULL,
    y = "Mean relative read abundance"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title.y       = element_text(size = 13, face = "bold"),
    axis.text.x        = element_text(size = 12, face = "bold"),
    axis.text.y        = element_text(size = 12, face = "bold"),
    axis.line          = element_line(linewidth = 0.7, colour = "black"),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.4),
    legend.position    = "bottom",
    legend.text        = element_text(size = 12, face = "bold.italic"),
    plot.margin        = margin(5, 10, 5, 5)
  )

print(p_ra_sub)
ggsave(
  file.path(FIG_DIR, "S3D_ra_substrate_bar.png"),
  p_ra_sub, width = 8, height = 5, dpi = 200
)
cat("Figure 3 saved: S3D_ra_substrate_bar.png\n\n")

cat("=== RELATIVE ABUNDANCE SUPPLEMENT COMPLETE ===\n\n")