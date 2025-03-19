scoring_and_predict_auc = function(expr, pheno, genesets, sample_colname, group_colname,
                                   dataset_title=""){
  require(ROCR)
  require(preprocessCore)
  require(e1071)
  require(GSVA)
  require(pROC)
  # Step 1: Prepare the data using GSVA
  gsvapar = ssgseaParam(as.matrix(expr), geneSets = genesets, normalize = T)
  gsva.exp = gsva(gsvapar, verbose = TRUE, BPPARAM = BiocParallel::SerialParam(progressbar = TRUE))
  gsva.exp = as.data.frame(t(gsva.exp))
  gsva.exp[,(ncol(gsva.exp) + 1)] = rownames(gsva.exp)
  names(gsva.exp)[ncol(gsva.exp)] = sample_colname
  gsva.exp = merge(gsva.exp, pheno, by = sample_colname)
  # Step 2: Fit a model (logistic regression in this example)
  # Ensure group_colname is a binary factor
  gsva.exp[[group_colname]] <- factor(gsva.exp[[group_colname]], levels = unique(gsva.exp[[group_colname]]))
  levels(gsva.exp[[group_colname]]) <- c(0, 1)
  formula = as.formula(paste(group_colname, "~", names(genesets)))
  model = glm(formula = formula, data = gsva.exp, family = "binomial")
  # Step 3: Calculate AUC
  actual_responses = gsva.exp[[group_colname]]
  predicted_probs = predict(model, type = "response")
  par(pty = "s")
  p.roc = plot.roc(actual_responses, predicted_probs,
                   legacy.axes = T, percent = T, col = "#A30023",
                   cex.lab = 2,
                   cex.axis = 1.7, cex.main = 2,
                   lwd = 4, main = paste0(dataset_title," ROC Curve"))
  auc_obj = p.roc$auc
  auc_decimal <- auc_obj / 100
  text(30, 20, paste("AUC =", round(auc_decimal, 4)), cex = 2.5, col = "#A30023")
  par(pty = "m")
  return(gsva.exp, p.roc)
}