#################################################################################
# TITLE: Anatomical Paper
# PURPOSE: Craniometric Analyses- Refined
# AUTHOR: KRK, MCS, AMP, RRD, JTH
# LAST UPDATED ON 3/17/2021
# UPDATES: 
# 1. Clean Code
# 2. Refined analysis

# NOTES: 
################################################################################

library(HDMD)
library(mnormt)
library(psych)
library(MASS)
library(readr)
library(tibble)
library(fansi)
library(mice)
library(ggplot2)
library(ellipse)
library(FactoMineR)
library(factoextra)
library(missMDA)

set.seed(1234)

read_or_get_data <- function(object_name, file_name = paste0(object_name, ".csv")) {
  if (file.exists(file_name)) {
    return(as.data.frame(readr::read_csv(file_name, show_col_types = FALSE)))
  }
  if (exists(object_name, inherits = TRUE)) {
    return(as.data.frame(get(object_name, inherits = TRUE)))
  }
  stop(paste("Could not find", object_name, "or", file_name))
}

rename_nlw <- function(x) {
  if ("NLW" %in% names(x) && !"NLB" %in% names(x)) {
    names(x)[names(x) == "NLW"] <- "NLB"
  }
  x
}

make_ellipse_data <- function(dat_pred) {
  ellipse_data <- lapply(unique(dat_pred$Group), function(g) {
    d <- dat_pred[dat_pred$Group == g, ]
    if (nrow(d) < 3 || any(!is.finite(d$LD1)) || any(!is.finite(d$LD2)) || sd(d$LD1) == 0 || sd(d$LD2) == 0) {
      return(NULL)
    }
    cbind(
      as.data.frame(with(d, ellipse(cor(LD1, LD2), scale = c(sd(LD1), sd(LD2)), centre = c(mean(LD1), mean(LD2))))),
      Group = g
    )
  })
  do.call(rbind, ellipse_data)
}

plot_lda_results <- function(data, model_formula, output_file) {
  lda_plot <- lda(model_formula, data = data)
  pred <- predict(lda_plot)
  dat_pred <- data.frame(Group = pred$class, pred$x)
  dat_ell <- make_ellipse_data(dat_pred)
  p <- ggplot(dat_pred, aes(x = LD1, y = LD2, col = Group)) +
    geom_point(size = 4.5, shape = 20, alpha = 0.75) +
    theme_classic()
  if (!is.null(dat_ell)) {
    p <- p + geom_path(data = dat_ell, aes(x = x, y = y, color = Group), size = 1, linetype = 2)
  }
  ggsave(filename = output_file, plot = p, width = 16, height = 12, units = "cm", scale = 1.6, quality = 100)
  p
}

run_mahalanobis <- function(data, model_formula, output_file) {
  variables <- all.vars(model_formula)[-1]
  grouping <- data$Group
  data_md <- as.matrix(data[, variables])
  mahala_sq <- pairwise.mahalanobis(data_md, grouping = grouping, center = TRUE, inverted = FALSE)
  d2 <- sqrt(mahala_sq$distance)
  group_names <- rownames(mahala_sq$means)
  rownames(d2) <- group_names
  colnames(d2) <- group_names
  write.csv(d2, output_file)
  invisible(d2)
}

plot_individual_mahalanobis <- function(data, model_formula, title = "Squared Mahalanobis distances") {
  variables <- all.vars(model_formula)[-1]
  data_md <- as.matrix(data[, variables])
  sx <- cov(data_md)
  d2 <- mahalanobis(data_md, colMeans(data_md), sx)
  plot(density(d2, bw = 0.5), main = title)
  rug(d2)
  invisible(d2)
}

run_lda_analysis <- function(reference_data, anatomical_data, combined_data, model_formula, cv_file, ccr_file, indiv_file, post_file, md_file, plot_file) {
  reference_data <- rename_nlw(as.data.frame(reference_data))
  anatomical_data <- rename_nlw(as.data.frame(anatomical_data))
  combined_data <- rename_nlw(as.data.frame(combined_data))
  reference_data$Group <- as.factor(reference_data$Group)
  anatomical_data$Group <- as.factor(anatomical_data$Group)
  combined_data$Group <- as.factor(combined_data$Group)
  lda_cv <- lda(model_formula, data = reference_data, CV = TRUE)
  cv_table <- table(reference_data$Group, lda_cv$class)
  write.csv(cv_table, cv_file)
  ccr <- diag(prop.table(cv_table, 1))
  write.csv(ccr, ccr_file)
  lda_fit <- lda(model_formula, data = reference_data)
  pred <- predict(lda_fit, newdata = anatomical_data)
  lda_indiv <- data.frame(original = anatomical_data$ID, pred = pred$class)
  write.csv(lda_indiv, indiv_file, row.names = FALSE)
  write.csv(pred$posterior, post_file)
  plot_obj <- plot_lda_results(combined_data, model_formula, plot_file)
  md <- run_mahalanobis(combined_data, model_formula, md_file)
  individual_md <- plot_individual_mahalanobis(combined_data, model_formula)
  list(
    lda_cv = lda_cv,
    cv_table = cv_table,
    ccr = ccr,
    overall_ccr = sum(diag(prop.table(cv_table))),
    lda_fit = lda_fit,
    predictions = pred,
    lda_indiv = lda_indiv,
    plot = plot_obj,
    mahalanobis = md,
    individual_mahalanobis = individual_md
  )
}

prepare_anatomical_data <- function() {
  data_anat <- read_or_get_data("MSUANAT", "MSUANAT.csv")
  data_anat <- as.data.frame(complete(mice(data_anat, m = 5)))
  data_anat_003 <- read_or_get_data("ANAT003", "ANAT003.csv")
  data_anat_combined <- as.data.frame(rbind(data_anat, data_anat_003))
  write.csv(data_anat_combined, "data.anat.comb.csv", row.names = FALSE)
  data_anat_combined
}

data.anat.comb <- prepare_anatomical_data()

broad_formula <- Group ~ GOL + BNL + BBH + XCB + XFB + ZYB + ASB + BPL + NLH + NLB + MAB + OBH + OBB + STB + FOL + UFHT
refined_formula <- Group ~ GOL + BNL + BBH + XCB + ZYB + ASB + BPL + NLH + NLB + MAB + OBH + OBB + STB + FOL + UFHT

data.ref.broad <- rename_nlw(read_or_get_data("Anat_Imputed_Data", "Anat_Imputed_Data.csv"))
data.comp.broad <- rbind(data.anat.comb, data.ref.broad)

broad.results <- run_lda_analysis(
  reference_data = data.ref.broad,
  anatomical_data = data.anat.comb,
  combined_data = data.comp.broad,
  model_formula = broad_formula,
  cv_file = "lda.ref.csv",
  ccr_file = "lda.ref.ccr.csv",
  indiv_file = "lda.indiv.csv",
  post_file = "lda.post.csv",
  md_file = "md.csv",
  plot_file = "LDFA-CODED-broad.jpg"
)

data.ref.refined.ind.only <- rename_nlw(read_or_get_data("Anat_Imputed_Refined_Data_Ind_Only", "Anat_Imputed_Refined_Data_Ind_Only.csv"))
data.comp.refined.ind.only <- rbind(data.anat.comb, data.ref.refined.ind.only)

refined.ind.only.results <- run_lda_analysis(
  reference_data = data.ref.refined.ind.only,
  anatomical_data = data.anat.comb,
  combined_data = data.comp.refined.ind.only,
  model_formula = refined_formula,
  cv_file = "lda.ref.no.xfb.csv",
  ccr_file = "lda.ref.ccr.no.xfb.csv",
  indiv_file = "lda.indiv.no.xfb.csv",
  post_file = "lda.post.no.xfb.csv",
  md_file = "md.no.xfb.csv",
  plot_file = "LDFA-CODED-refined-ind-only.jpg"
)

data.ref.refined.all <- rename_nlw(read_or_get_data("Anat_Imputed_Refined_Data", "Anat_Imputed_Refined_Data.csv"))
data.comp.refined.all <- rbind(data.anat.comb, data.ref.refined.all)

refined.all.results <- run_lda_analysis(
  reference_data = data.ref.refined.all,
  anatomical_data = data.anat.comb,
  combined_data = data.comp.refined.all,
  model_formula = refined_formula,
  cv_file = "lda.ref.2.no.xfb.csv",
  ccr_file = "lda.ref.ccr.2.no.xfb.csv",
  indiv_file = "lda.indiv.3.no.xfb.csv",
  post_file = "lda.post.3.no.xfb.csv",
  md_file = "md.2.no.xfb.csv",
  plot_file = "LDFA-CODED-refined-all.jpg"
)

data.ref.howells.only <- rename_nlw(read_or_get_data("Anat_Imputed_Refined_Data_How_Only", "Anat_Imputed_Refined_Data_How_Only.csv"))
howells.only.plot <- plot_lda_results(data.ref.howells.only, refined_formula, "LDFA-CODED-howells-only.jpg")

if (file.exists("data.csv")) {
  famd.data <- read.csv("data.csv", head = TRUE, check.names = FALSE)
  anat.ids <- c("ANAT-001", "ANAT-002", "ANAT-003", "ANAT-004", "ANAT-005", "ANAT-006", "ANAT-007", "ANAT-008", "ANAT-009", "ANAT-010", "ANAT-011", "ANAT-012", "ANAT-013", "ANAT-014", "ANAT-019", "ANAT-020")
  rownames(famd.data)[seq_along(anat.ids)] <- anat.ids
  famd.df <- famd.data[1:16, 2:65]
  numeric.cols <- 4:min(64, ncol(famd.df))
  famd.df[numeric.cols] <- lapply(famd.df[numeric.cols], as.numeric)
  famd.ncp <- min(64, ncol(famd.df), nrow(famd.df) - 1)
  famd.impute.ncp <- min(5, famd.ncp)
  famd.impute <- imputeFAMD(famd.df, ncp = famd.impute.ncp)
  res.famd <- FAMD(famd.df, tab.disj = famd.impute$tab.disj, ncp = famd.ncp, graph = FALSE)
  eig.val <- get_eigenvalue(res.famd)
  famd.var <- get_famd_var(res.famd)
  famd.quanti.var <- get_famd_var(res.famd, "quanti.var")
  famd.quali.var <- get_famd_var(res.famd, "quali.var")
  famd.ind <- get_famd_ind(res.famd)
  famd.variable.contributions <- as.data.frame(res.famd$var$contrib)
  write.csv(eig.val, "famd.eigenvalues.csv")
  write.csv(famd.var$coord, "famd.variable.coordinates.csv")
  write.csv(famd.var$cos2, "famd.variable.cos2.csv")
  write.csv(famd.var$contrib, "famd.variable.contributions.csv")
  famd.scree.plot <- fviz_screeplot(res.famd)
  famd.variable.plot <- fviz_famd_var(res.famd, col.var = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, title = "Variable Associations")
  famd.dim1.plot <- fviz_contrib(res.famd, "var", axes = 1, top = 35, color = "#00AFBB", fill = "#00AFBB", title = "Contribution of Variables to Dimension 1")
  famd.dim2.plot <- fviz_contrib(res.famd, "var", axes = 2, top = 25, color = "#E7B800", fill = "#E7B800", title = "Contribution of Variables to Dimension 2")
  famd.quanti.plot <- fviz_famd_var(res.famd, "quanti.var", col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, title = "Quantitative Variables Correlations")
  famd.quanti.cos2.plot <- fviz_famd_var(res.famd, "quanti.var", col.var = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE)
  famd.quali.plot <- fviz_famd_var(res.famd, "quali.var", col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))
  famd.individual.plot <- fviz_famd_ind(res.famd, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), repel = TRUE, title = "Individuals & Qualitative Biplot")
}
