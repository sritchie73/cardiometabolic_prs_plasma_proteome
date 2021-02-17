# Aggregate aptamer P-values into a single p-value for the protein,
# handling directionality of effect sizes
prot_pvalue <- function(pval, beta, avg.func=mean) {
  if(length(pval) == 1) {
    # When there is only 1 aptamer, simply return the p-value.
    return(pval)
  } else {
    # First, we need to convert two-sided p-values to one-sided p-values,
    # which we'll orient to the direction of effect of the smallest PRS to
    # aptamer p-value for this protein.
    one.sided <- pval/2
    to.invert <- sign(beta) != sign(beta[which.min(pval)])
    one.sided[to.invert] <- 1 - one.sided[to.invert]

    # Convert to Z-scores
    zs <- qnorm(one.sided, lower.tail=FALSE)

    # Average Z-scores
    avg.z <- avg.func(zs)

    # Sometimes the avg.z may be NaN, e.g. where taking the mean where there are P-values of 0.
    if (is.nan(avg.z)) {
      if (length(unique(sign(zs)))) {
        warning("Avg Z-score is NaN due to input P-values of 0. Beta effects do not have consistent direction so returning P=0.5")
        return(0.5)
      } else {
        warning("Avg Z-score is NaN due to input P-values of 0. Beta effects have consistent direction so returning P=0")
        return(0)
      }
    }

    # Convert back to p-value
    meta.one.sided <- pnorm(avg.func(zs), lower.tail=FALSE)

    # Convert to two-sided p-value
    meta.two.sided <- ifelse(meta.one.sided > 0.5, (1-meta.one.sided)*2, meta.one.sided*2)

    return(meta.two.sided)
  }
}


