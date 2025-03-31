# Loading required libraries
library(glmnet)
library(survival)
library(dplyr)
library(survminer)
library(timeROC)
library(Hmisc)

# ─────────────────────────────────────────────────────────────
# MODEL: Elastic Net Cox Regression Model
# ─────────────────────────────────────────────────────────────
# This model uses:
# - Model Type: Cox proportional hazards model (for survival analysis)
# - Penalty Applied: Elastic Net regularization
# - Alpha = 0.5 (equal blend of LASSO and Ridge)
# - Lambda = determined by 10-fold cross-validation (lambda.min)
# - Features: Top 500 survival-associated genes + Age, Sex, BCG

# Loading the datasets
uromol <- readRDS("UROMOL_TaLG.teachingcohort (1).rds")
knowles <- readRDS("knowles_matched_TaLG_final.rds")

# First, identifying shared genes
shared_genes <- intersect(colnames(uromol$exprs), colnames(knowles$exprs))

# Normalizing the Z-scores for gene expression within UROMOL
exprs_u <- uromol$exprs[, shared_genes]
exprs_u_z <- scale(exprs_u)

# Applying univariate Cox regression to filter only those genes that are associated with recurrence free survival
gene_pvals <- apply(exprs_u_z[, shared_genes], 2, function(g) {
  coxph_fit <- try(coxph(Surv(uromol$RFS_time, uromol$Recurrence) ~ g), silent = TRUE)
  if (inherits(coxph_fit, "try-error")) return(1)
  summary(coxph_fit)$coefficients[5]
})

# Selecting the top 500 most survival-associated genes
top500_cox_genes <- names(sort(gene_pvals))[1:500]
exprs_u_cox <- exprs_u_z[, top500_cox_genes]

# Subsetting clinical features: Age, Sex, BCG
clinical_u <- uromol %>%
  select(RFS_time, Recurrence, Age, Sex, BCG) %>%
  filter(complete.cases(.)) %>%
  mutate(Sex = factor(Sex))

# Filtering expression to match clinical data
exprs_u_cox_filtered <- exprs_u_cox[rownames(clinical_u), ]

# Designing the matrix and outcome for the model
y_train <- with(clinical_u, Surv(RFS_time, Recurrence))
X_clinical <- model.matrix(~ Age + Sex + BCG, data = clinical_u)[, -1]
X_train_cox <- cbind(X_clinical, exprs_u_cox_filtered)

# Fitting the Elastic Net Cox model
set.seed(123)
cv_fit_cox <- cv.glmnet(x = X_train_cox, y = y_train, family = "cox", alpha = 0.5)

# Predicting the risk scores and stratifying patients
risk_cox <- predict(cv_fit_cox, newx = X_train_cox, s = "lambda.min", type = "link")
clinical_u$risk_score <- risk_cox
clinical_u$risk_group <- ifelse(risk_cox > median(risk_cox), "High", "Low")

# Kaplan-Meier plot
fit_cox <- survfit(Surv(RFS_time, Recurrence) ~ risk_group, data = clinical_u)
ggsurvplot(fit_cox, data = clinical_u, pval = TRUE, risk.table = TRUE,
           title = "KM: Top 500 Survival-Associated Genes")

# Evaluating AUC on UROMOL
# ─────────────────────────────────────────────────────────────
# Diagnostic check to ensure consistent training data
cat("Training samples (X_train_cox):", nrow(X_train_cox), "
")
cat("Training samples (y_train):", length(y_train), "
")
cat("Events (Recurrence = 1):", sum(clinical_u$Recurrence == 1), "
")
roc_cox <- timeROC(T = clinical_u$RFS_time, delta = clinical_u$Recurrence,
                   marker = clinical_u$risk_score, cause = 1,
                   times = c(12, 24, 36), iid = TRUE)

# Integrated AUC across time points and mean of AUCs
print(roc_cox$AUC)

auc_mean <- mean(roc_cox$AUC)
cat("
Mean AUC (12/24/36 months):", round(auc_mean, 3), "
")

# ─────────────────────────────────────────────────────────────
# VALIDATION ON KNOWLES USING ONLY SHARED COLUMNS
# ─────────────────────────────────────────────────────────────

# Subsetting and normalizing z-score Knowles gene expression
exprs_k <- knowles$exprs[, top500_cox_genes]
exprs_k_z <- scale(exprs_k)

# Preparing the Knowles clinical data
knowles_clinical <- knowles %>%
  select(Progression, PFS_time., Recurrence, RFS_time, Age, Sex, BCG) %>%
  mutate(
    Sex = factor(Sex, levels = levels(clinical_u$Sex)),
    BCG = as.numeric(BCG),
    RFS_time = as.numeric(as.character(RFS_time)),
    Recurrence = as.numeric(as.character(Recurrence))
  )
knowles_clinical <- knowles_clinical[rownames(exprs_k_z), ]

# Designing the matrix for validation
Xk_clinical <- model.matrix(~ Age + Sex + BCG, data = knowles_clinical)[, -1]
X_valid_cox <- cbind(Xk_clinical, exprs_k_z)

# Matching columns with training matrix
common_cols <- intersect(colnames(X_train_cox), colnames(X_valid_cox))
X_valid_cox <- X_valid_cox[, common_cols]
X_train_cox <- X_train_cox[, common_cols]

# Predicting the risk scores in validation cohort
knowles_clinical$risk_score <- predict(cv_fit_cox, newx = X_valid_cox, s = "lambda.min", type = "link")
knowles_clinical$risk_group <- ifelse(knowles_clinical$risk_score > median(knowles_clinical$risk_score), "High", "Low")

# Filtering only complete cases for RFS
knowles_clinical_complete <- knowles_clinical %>%
  filter(!is.na(RFS_time) & !is.na(Recurrence) & is.finite(RFS_time))

# Validation with  KM plot
if (nrow(knowles_clinical_complete) > 0) {
  fit_val <- survfit(Surv(RFS_time, Recurrence) ~ risk_group, data = knowles_clinical_complete)
  
  km_val_plot <- ggsurvplot(
    fit_val,
    data = knowles_clinical_complete,
    pval = TRUE,
    risk.table = TRUE,
    title = "Validation Cohort: Kaplan-Meier Plot"
  )
  print(km_val_plot)
  
  # Time-dependent AUC
  roc_val <- timeROC(
    T = knowles_clinical_complete$RFS_time,
    delta = knowles_clinical_complete$Recurrence,
    marker = knowles_clinical_complete$risk_score,
    cause = 1,
    times = c(12, 24, 36),
    iid = TRUE
  )
  print(roc_val$AUC)
  
 
