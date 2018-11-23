# This function requires:
#    i: index (e.g. the position (raw) of the CpG to be analyzed)
#    model: a formula. Left side the outcome, right side covariates
#    molecular data: a matrix having features in rows and samples in columns
#    pheno.data: a data.frame having samples in rows and variables in columns

lmFeature <- function(i, model, molecular.data, pheno.data){
  sel <- rownames(molecular.data)[i]
  mt <- as.formula(model)
  
  vars <- all.vars(mt)
  condition <- vars[1]
  covars <- vars[-1]
  
  data.sel <- data.frame(feature=molecular.data[sel,], 
                         pheno.data[, vars])
  
  dep <- paste(vars, collapse="+")
  mm <- as.formula(paste("feature ~", dep))
  mod <- lm(mm, data=data.sel)
  
  ans <- matrix(summary(mod)$coefficients[c(1,2,4)], nrow=1)
  colnames(ans) <- c("beta", "s.e.", "p.value")
  rownames(ans) <- sel
  return(ans)
}