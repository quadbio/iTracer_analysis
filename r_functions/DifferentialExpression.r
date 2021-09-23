# ANOVA-based test for expression changed linearly related to a given factor/numeric variable
ancova_group_test <- function(expr, # expression matrix, rows as genes
                              group, # group values, factor or numeric
                              covar = NULL, # covariate data frame
                              num_threads = 20, # number of threads to use
                              return_coef_group = NULL) # the level of factor in the group to return the coefficient
{
  require(doParallel, quietly = T)
  idx <- which(!is.na(group))
  expr <- expr[,idx]
  
  group <- group[idx]
  if (is.factor(group) | is.character(group)){
    group <- factor(group)
    g_lev <- levels(group)
    if (!is.null(return_coef_group) & sum(g_lev == return_coef_group)>0){
      g_lev <- c(setdiff(g_lev, return_coef_group), return_coef_group)
    } else{
      return_coef_group <- NULL
    }
    group <- factor(group, levels = g_lev)
  }
  
  if(!is.null(covar)) covar <- setNames(data.frame(covar[idx,]), colnames(covar))
  
  registerDoParallel(num_threads)
  res <- foreach(i = 1:nrow(expr), .combine = rbind) %dopar%{
    e <- as.numeric(expr[i,])
    
    if (is.null(covar)){
      m0 <- lm(e ~ 1)
      m1 <- lm(e ~ group)
    } else{
      data <- data.frame(y = e, covar, group = group)
      m0 <- lm(y ~ . - group, data = data)
      m1 <- lm(y ~ ., data = data)
    }
    
    a0 <- anova(m0)
    a1 <- anova(m1)
    a01 <- anova(m0, m1)
    
    p1 <- a1$Pr[nrow(a1)-1]
    p2 <- a01$Pr[2]
    p3 <- pf(a0["Residuals","Mean Sq"] / a1["Residuals","Mean Sq"], df1 = a0["Residuals","Df"], df2 = a1["Residuals","Df"], lower.tail = F)
    
    var_tot <- sum(a0[,2])
    var_covar <- sum(a0[-nrow(a0),2])
    var_g <- a1[nrow(a1)-1,2]
    
    if (is.factor(group)){
      coef_g <- ifelse(is.null(return_coef_group), NA, coef(m1)[paste0("group",return_coef_group)])
    } else{
      coef_g <- coef(m1)["group"]
    }
    
    return(c(var_g, var_covar, var_tot, p1, p2, p3, coef_g))
  }
  stopImplicitCluster()
  res <- data.frame(res)
  rownames(res) <- rownames(expr)
  colnames(res) <- c("var_group","var_covar","var_total","p_ANOVA", "p_ANCOVA", "p_Resi", "coef")
  return(res)
}
