library(MendelianRandomization)
library(data.table)
library(scales)
library(ggplot2)
library(ggrepel)

# Try all MR methods
#
# Equivalent to to MendelianRandomization::mr_allmethods()
# but is robust to failure of assumptions (e.g. too few variants).
# All methods that can be run, will be run, rather than erroring
# out if one method fails.
#
# @param mri object of class MRInput()
#
# @return a data.table of results equivalent to whats shown in
#         show(mr_allmethods()).
mr_tryall <- function(mri) {
  # Set up failure data.tables
  mr_fail_dt <- data.table(mr_estimate = NA_real_)
  mr_fail_dt[, c("mr_se", "mr_L95", "mr_U95", "mr_pval") := NA_real_]
  
	mr_ivw_dt <- tryCatch({
    # MendelianRandomization functions behave badly, issuing messages with cat() 
    # instead of message() making for awkward message suppression:
		msg <- capture.output({ mro <- mr_ivw(mri) })
		mr_ivw_dt <- data.table(method = "IVW", 
														mr_estimate = mro@Estimate,
														mr_se = mro@StdError, 
														mr_L95 = mro@CILower, 
														mr_U95 = mro@CIUpper,
														mr_pval = mro@Pvalue)
	}, error = function(e) { return(cbind(data.table(method = "IVW"), mr_fail_dt)) })
	mro_dt <- mr_ivw_dt

	mr_ivw_penalised <- tryCatch({
		msg <- capture.output({ mro <- mr_ivw(mri, penalized = TRUE) })
		mr_ivw_penalised <- data.table(method = "Penalized IVW", 
																	 mr_estimate = mro@Estimate,
																	 mr_se = mro@StdError, 
																	 mr_L95 = mro@CILower, 
																	 mr_U95 = mro@CIUpper,
																	 mr_pval = mro@Pvalue)
	}, error = function(e) { return(cbind(data.table(method = "Penalized IVW"), mr_fail_dt)) })
	mro_dt <- rbind(mro_dt, mr_ivw_penalised)

	mr_ivw_robust <- tryCatch({
		msg <- capture.output({ mro <- mr_ivw(mri, robust=TRUE) })
		mr_ivw_robust <- data.table(method = "Robust IVW", 
																mr_estimate = mro@Estimate,
																mr_se = mro@StdError, 
															  mr_L95 = mro@CILower, 
																mr_U95 = mro@CIUpper,
																mr_pval = mro@Pvalue)
	}, # Warning can be generated when the robust lm fails to converge - results are
		 # returned but its not clear to me how reliable these are. Given we're doing an
		 # MR scan, better to just throw these out.
	   warning = function(e) { return(cbind(data.table(method = "Robust IVW"), mr_fail_dt)) },
	   error   = function(e) { return(cbind(data.table(method = "Robust IVW"), mr_fail_dt)) })
	mro_dt <- rbind(mro_dt, mr_ivw_robust)

	mr_ivw_robust_penalised <- tryCatch({
		msg <- capture.output({ mro <- mr_ivw(mri, robust = TRUE, penalized = TRUE) })
		mr_ivw_robust_penalised <- data.table(method = "Penalized robust IVW", 
																					mr_estimate = mro@Estimate,
																					mr_se = mro@StdError, 
																					mr_L95 = mro@CILower, 
																					mr_U95 = mro@CIUpper,
																					mr_pval = mro@Pvalue)
  }, warning = function(e) { return(cbind(data.table(method = "Penalized robust IVW"), mr_fail_dt)) },
     error   = function(e) { return(cbind(data.table(method = "Penalized robust IVW"), mr_fail_dt)) })
	mro_dt <- rbind(mro_dt, mr_ivw_robust_penalised)

	mr_median_simple <- tryCatch({
		msg <- capture.output({ mro <- mr_median(mri, weighting = "simple") })
		mr_median_simple <- data.table(method = "Simple median", 
																	 mr_estimate = mro@Estimate,
																	 mr_se = mro@StdError, 
																	 mr_L95 = mro@CILower, 
																	 mr_U95 = mro@CIUpper,
																	 mr_pval = mro@Pvalue)
	}, error = function(e) { return(cbind(data.table(method = "Simple median"), mr_fail_dt)) })
	mro_dt <- rbind(mro_dt, mr_median_simple)

	mr_median_weighted <- tryCatch({
		msg <- capture.output({ mro <- mr_median(mri, weighting = "weighted") })
		mr_median_weighted <- data.table(method = "Weighted median", 
																		 mr_estimate = mro@Estimate,
																		 mr_se = mro@StdError, 
																		 mr_L95 = mro@CILower, 
																		 mr_U95 = mro@CIUpper,
																		 mr_pval = mro@Pvalue)
	}, error = function(e) { return(cbind(data.table(method = "Weighted median"), mr_fail_dt)) })
	mro_dt <- rbind(mro_dt, mr_median_weighted)

	mr_median_penalized <- tryCatch({
		msg <- capture.output({ mro <- mr_median(mri, weighting = "penalized") })
		mr_median_penalized <- data.table(method = "Penalized weighted median", 
																			mr_estimate = mro@Estimate,
																			mr_se = mro@StdError, 
																			mr_L95 = mro@CILower, 
																			mr_U95 = mro@CIUpper,
																			mr_pval = mro@Pvalue)
	}, error = function(e) { return(cbind(data.table(method = "Penalized weighted median"), mr_fail_dt)) })
	mro_dt <- rbind(mro_dt, mr_median_penalized)

	mr_mode_weighted_simple <- tryCatch({
		msg <- capture.output({ mro <- mr_mbe(mri, weighting = "weighted", stderror = "simple") })
		mr_mode_weighted_simple <- data.table(method = "Weighted mode (simple SE)", 
																				  mr_estimate = mro@Estimate,
																				  mr_se = mro@StdError, 
																				  mr_L95 = mro@CILower, 
																				  mr_U95 = mro@CIUpper,
																				  mr_pval = mro@Pvalue)
	}, error = function(e) { return(cbind(data.table(method = "Weighted mode (simple SE)"), mr_fail_dt)) })
	mro_dt <- rbind(mro_dt, mr_mode_weighted_simple)

  mr_mode_weighted_delta <- tryCatch({
    msg <- capture.output({ mro <- mr_mbe(mri, weighting = "weighted", stderror = "delta") })
    mr_mode_weighted_delta <- data.table(method = "Weighted mode (delta SE)",
                                         mr_estimate = mro@Estimate,
                                         mr_se = mro@StdError,
                                         mr_L95 = mro@CILower,
                                         mr_U95 = mro@CIUpper,
                                         mr_pval = mro@Pvalue)
  }, error = function(e) { return(cbind(data.table(method = "Weighted mode (delta SE)"), mr_fail_dt)) })
  mro_dt <- rbind(mro_dt, mr_mode_weighted_delta)

	mr_mode_unweighted_simple <- tryCatch({
		msg <- capture.output({ mro <- mr_mbe(mri, weighting = "unweighted", stderror = "simple") })
		mr_mode_unweighted_simple <- data.table(method = "Unweighted mode (simple SE)", 
																						mr_estimate = mro@Estimate,
																						mr_se = mro@StdError, 
																						mr_L95 = mro@CILower, 
																						mr_U95 = mro@CIUpper,
																						mr_pval = mro@Pvalue)
	}, error = function(e) { return(cbind(data.table(method = "Unweighted mode (simple SE)"), mr_fail_dt)) })
	mro_dt <- rbind(mro_dt, mr_mode_unweighted_simple)

  mr_mode_unweighted_delta <- tryCatch({
    msg <- capture.output({ mro <- mr_mbe(mri, weighting = "unweighted", stderror = "delta") })
    mr_mode_unweighted_delta <- data.table(method = "Unweighted mode (delta SE)",
																					 mr_estimate = mro@Estimate,
																					 mr_se = mro@StdError,
																					 mr_L95 = mro@CILower,
																					 mr_U95 = mro@CIUpper,
																					 mr_pval = mro@Pvalue)
  }, error = function(e) { return(cbind(data.table(method = "Unweighted mode (delta SE)"), mr_fail_dt)) })
  mro_dt <- rbind(mro_dt, mr_mode_unweighted_delta)

##   mr_hetpen_dt <- tryCatch({
##     msg <- capture.output({ mro <- mr_hetpen(mri) }) 
##     # note can have multiple CIs for different groups of variants, 
##     # however only the biggest group/best estimate will have a value in
##     # the Estimate slot. There are no P-values for this method
##     mr_hetpen_dt <- data.table(method = "Heterogeneity penalized mode",
## 															 mr_estimate = NA_real_,
##                                mr_se = NA_real_,
## 															 mr_L95 = mro@CILower,
## 															 mr_U95 = mro@CIUpper,
## 															 mr_pval = NA_real_)
##     mr_hetpen_dt[1, mr_estimate := mro@Estimate]
##     if (nrow(mr_hetpen_dt) > 1) {
##       mr_hetpen_dt[, method := sprintf("%s (%s)", method, .I)]
##     }
##     mr_hetpen_dt
##   }, error = function(e) { return(cbind(data.table(method = "Heterogeneity penalized mode"), mr_fail_dt)) })
##   mro_dt <- rbind(mro_dt, mr_hetpen_dt)

	mr_egger_dt <- tryCatch({
		msg <- capture.output({ mro <- mr_egger(mri) })
		mr_egger_dt <- data.table(method = c("MR-Egger", "(intercept)"), 
														  mr_estimate = c(mro@Estimate, mro@Intercept),
															mr_se = c(mro@StdError.Est, mro@StdError.Int),
															mr_L95 = c(mro@CILower.Est, mro@CILower.Int),
															mr_U95 = c(mro@CIUpper.Est, mro@CIUpper.Int),
															mr_pval = c(mro@Causal.pval, mro@Pleio.pval)) # also mro@Pval.Est, mro@Pval.Int
	}, error = function(e) { return(cbind(data.table(method = c("MR-Egger", "(intercept)")), mr_fail_dt)) })
	mro_dt <- rbind(mro_dt, mr_egger_dt)

	mr_egger_penalised <- tryCatch({
		msg <- capture.output({ mro <- mr_egger(mri, penalized = TRUE) })
		mr_egger_penalised <- data.table(method = c("Penalized MR-Egger", "(intercept)"), 
																		 mr_estimate = c(mro@Estimate, mro@Intercept),
																		 mr_se = c(mro@StdError.Est, mro@StdError.Int),
																		 mr_L95 = c(mro@CILower.Est, mro@CILower.Int),
																		 mr_U95 = c(mro@CIUpper.Est, mro@CIUpper.Int),
																		 mr_pval = c(mro@Causal.pval, mro@Pleio.pval)) # also mro@Pval.Est, mro@Pval.Int
	}, error = function(e) { return(cbind(data.table(method = c("Penalized MR-Egger", "(intercept)")), mr_fail_dt)) })
	mro_dt <- rbind(mro_dt, mr_egger_penalised)

	mr_egger_robust <- tryCatch({
		msg <- capture.output({ mro <- mr_egger(mri, robust = TRUE) })
		mr_egger_robust <- data.table(method = c("Robust MR-Egger", "(intercept)"), 
																  mr_estimate = c(mro@Estimate, mro@Intercept),
																	mr_se = c(mro@StdError.Est, mro@StdError.Int),
																	mr_L95 = c(mro@CILower.Est, mro@CILower.Int),
																	mr_U95 = c(mro@CIUpper.Est, mro@CIUpper.Int),
																	mr_pval = c(mro@Causal.pval, mro@Pleio.pval)) # also mro@Pval.Est, mro@Pval.Int
	}, warning = function(w) { return(cbind(data.table(method = c("Robust MR-Egger", "(intercept)")), mr_fail_dt)) },
		 error   = function(e) { return(cbind(data.table(method = c("Robust MR-Egger", "(intercept)")), mr_fail_dt)) })
	mro_dt <- rbind(mro_dt, mr_egger_robust)

	mr_egger_penalised_robust <- tryCatch({
		msg <- capture.output({ mro <- mr_egger(mri, robust = TRUE, penalized = TRUE) })
		mr_egger_penalised_robust <- data.table(method = c("Penalized robust MR-Egger", "(intercept)"),
																						mr_estimate = c(mro@Estimate, mro@Intercept),
																						mr_se = c(mro@StdError.Est, mro@StdError.Int),
																						mr_L95 = c(mro@CILower.Est, mro@CILower.Int),
																						mr_U95 = c(mro@CIUpper.Est, mro@CIUpper.Int),
																						mr_pval = c(mro@Causal.pval, mro@Pleio.pval)) # also mro@Pval.Est, mro@Pval.Int
	}, warning = function(w) { return(cbind(data.table(method = c("Penalized robust MR-Egger", "(intercept)")), mr_fail_dt)) },
		 error   = function(e) { return(cbind(data.table(method = c("Penalized robust MR-Egger", "(intercept)")), mr_fail_dt)) })
	mro_dt <- rbind(mro_dt, mr_egger_penalised_robust)

  return(mro_dt)
}

# Extract dose response fits
#
# @param mro_dt data.table returned by mr_tryall()
# @return a data.table containing the intercept and slope for each mr-estimate
dose_response_fits <- function(mro_dt) {
  eggers <- data.table(
    method = mro_dt[which(method == "(intercept)")-1, method],
    intercept = mro_dt[which(method == "(intercept)"), mr_estimate],
    mr_L95 = mro_dt[which(method == "(intercept)"), mr_L95],
    mr_U95 = mro_dt[which(method == "(intercept)"), mr_U95]
  )

  plot_dt <- mro_dt[method != "(intercept)", .(method, intercept=0, slope=mr_estimate)]
  plot_dt[eggers, on = .(method), intercept := i.intercept]
  plot_dt[, method := factor(method, levels=unique(method))]
  plot_dt <- plot_dt[!is.na(slope)]
  return(plot_dt)
}

# Extract the 95% Confidence Intervals for the dose response curves
#
# @param mro_dt data.table returned by mr_tryall()
# @return a data.table containing the 95% confidence intervals for the
#         intercept and slope for each mr-estimate
dose_response_ci95 <- function(mro_dt) {
  eggers <- data.table(
    method = mro_dt[which(method == "(intercept)")-1, method],
    intercept = mro_dt[which(method == "(intercept)"), mr_estimate],
    mr_L95 = mro_dt[which(method == "(intercept)"), mr_L95],
    mr_U95 = mro_dt[which(method == "(intercept)"), mr_U95]
  )

  ci95_dt <- mro_dt[method != "(intercept)", .(method,
    intercept.l95=0, intercept.u95=0,
    slope.l95=mr_L95, slope.u95=mr_U95)]
  ci95_dt[eggers, on = .(method), c("intercept.l95", "intercept.u95") := .(mr_L95, mr_U95)]
  ci95_dt[, method := factor(method, levels=unique(method))]
  ci95_dt <- ci95_dt[!is.na(slope.l95)]
  return(ci95_dt)
}

# Determine the X and Y limits for the dose response curve plots
#
# Plot will show all instruments, their standard errors, the x/y axis
# intercepts, and if applicable will include the 95% confidence intervals
# for the MR-Egger intercepts.
#
# @param instruments data.table containing information about each instrument
# @param ci95_dt data.table returned by dose_response_ci95().
# @return a list containing the xlim and ylim vectors 
#				  along with the plot width and plot height.
dose_response_plot_lim <- function(instruments, ci95_dt) {
  xlim <- c(min(instruments$Prot.Effect.pQTL - instruments$Prot.SE.pQTL),
            max(instruments$Prot.Effect.pQTL + instruments$Prot.SE.pQTL))
  ylim <- c(min(instruments$effect.gwas - instruments$se.gwas),
            max(instruments$effect.gwas + instruments$se.gwas))
  ylim <- c(min(c(ylim[1], ci95_dt$intercept.l95)), 
            max(c(ylim[2], ci95_dt$intercept.u95)))
	if (xlim[1] > 0) xlim[1] <- 0
	if (xlim[2] < 0) xlim[2] <- 0
	if (ylim[1] > 0) ylim[1] <- 0
	if (ylim[2] < 0) ylim[2] <- 0
	pw <- xlim[2] - xlim[1]
	ph <- ylim[2] - ylim[1]
	xlim[1] <- xlim[1] - pw * 0.2
	xlim[2] <- xlim[2] + pw * 0.2
	ylim[1] <- ylim[1] - ph * 0.2
	ylim[2] <- ylim[2] + ph * 0.2
  return(list(xlim=xlim, ylim=ylim, pw=pw, ph=ph))
}

#' Polygon points for showing line uncertainty.
#'
#' Given a line of best fit defined as y = a + bx, e.g. an estimate of 
#' causality / dose response curve, uncertainty around the intercept 
#' and slope of that line (e.g. 95% confidence intervals), generates a 
#' series of points that can be used to draw a polygon showing the range 
#' of possible lines of fit (e.g. via geom_polygon()).
#'
#' Warning: currently buggy if the intercept.l95 or intercept.u95 are not
## within the plot window.
#'
#' @param intercept.l95 lower 95% confidence interval or lower limit for the y-intercept
#' @param intercept.u95 upper 95% confidence interval or upper limit for the y-intercept
#' @param slope.l95 lower 95% condience interval or lower limit for the slope
#' @param slope.u95 upper 95% condience interval or upper limit for the slope
#' @param xlim vector of length two giving the x-axis plot limits.
#' @param ylim vector of length two giving the y-ayis plot limits.
#' 
#' @return a data.table of x and y coordinates used to define the polygon
line_ci95_poly <- function(
	intercept.l95, intercept.u95, 
  slope.l95, slope.u95,
  xlim, ylim
) {
  intersect_hline <- function(yintercept, slope, h) {
    # y = a + bx, y=c       becomes
    # c = a + bx
    # (c-a)/b = x
    if (slope == 0) ifelse(h > 1, Inf, -Inf)
    if (is.infinite(slope)) return(0)
    (h - yintercept)/slope
  }
  intersect_vline <- function(yintercept, slope, v) {
    # y = a + bx, x=c       becomes
    # y = a + b*c
    if (slope == 0) return(yintercept)
    if (is.infinite(slope)) ifelse(v > 1, Inf, -Inf)
    yintercept + slope * v
  }
  
  # define plot corners
  top_right <- data.table(x=xlim[2], y=ylim[2])
  bottom_right <- data.table(x=xlim[2], y=ylim[1])
  bottom_left <- data.table(x=xlim[1], y=ylim[1])
  top_left <- data.table(x=xlim[1], y=ylim[2])
  
  # Starting at the intercept.u95, we'll go around the plot 
  # finding the intersection points and adding them to the 
  # polygon data.table
  poly_dt <- data.table(x=0, y=intercept.u95)
  
  # (1) From the intercept.u95 to the bottom/right/top border based on the slope.u95
  top_intersect <- intersect_hline(intercept.u95, slope.u95, ylim[2]) 
  right_intersect <- intersect_vline(intercept.u95, slope.u95, xlim[2])
  bottom_intersect <- intersect_hline(intercept.u95, slope.u95, ylim[1])
  if (slope.u95 >= 0) {
    if (right_intersect > ylim[2]) {
      new_point <- data.table(x=top_intersect, y=ylim[2])
      this_axis <- "top"
    } else {
      new_point <- data.table(x=xlim[2], y=right_intersect)
      this_axis <- "right"
    }
  } else if (slope.u95 < 0) {
    if (right_intersect < ylim[1]) {
      new_point <- data.table(x=bottom_intersect, y=ylim[1])
      this_axis <- "bottom"
    } else {
      new_point <- data.table(x=xlim[2], y=right_intersect)
      this_axis <- "right"
    }
  }
  poly_dt <- rbind(poly_dt, new_point)
  
  # (2) From the intercept.l95 to the bottom/right/top border based on the slope.l95
  #     - we also have to add the plot corners if applicable.
  prev_axis <- this_axis
  top_intersect <- intersect_hline(intercept.l95, slope.l95, ylim[2]) 
  right_intersect <- intersect_vline(intercept.l95, slope.l95, xlim[2])
  bottom_intersect <- intersect_hline(intercept.l95, slope.l95, ylim[1])
  if (slope.l95 >= 0) {
    if (right_intersect > ylim[2]) {
      new_point <- data.table(x=top_intersect, y=ylim[2])
      this_axis <- "top"
    } else {
      new_point <- data.table(x=xlim[2], y=right_intersect)
      this_axis <- "right"
    }
  } else if (slope.l95 < 0) {
    if (right_intersect < ylim[1]) {
      new_point <- data.table(x=bottom_intersect, y=ylim[1])
      this_axis <- "bottom"
    } else {
      new_point <- data.table(x=xlim[2], y=right_intersect)
      this_axis <- "right"
    }
  }
  # add corners if applicable
  if (prev_axis == "top" && this_axis == "right") {
    poly_dt <- rbind(poly_dt, top_right)
  } else if (prev_axis == "top" && this_axis == "bottom") {
    poly_dt <- rbind(poly_dt, top_right, bottom_right)
  } else if (prev_axis == "right" && this_axis == "bottom") {
    poly_dt <- rbind(poly_dt, bottom_right)
  }
  # add new point
  poly_dt <- rbind(poly_dt, new_point)
  
  # back to y-intercept
  poly_dt <- rbind(poly_dt, data.table(x=0, y=intercept.l95))
  
  # (3) From the intercept.l95 to the bottom/left/top border based on the slope.u95
  bottom_intersect <- intersect_hline(intercept.l95, slope.u95, ylim[1])
  left_intersect <- intersect_vline(intercept.l95, slope.u95, xlim[1])
  top_intersect <- intersect_hline(intercept.l95, slope.u95, ylim[2]) 
  if (slope.u95 >= 0) {
    if (left_intersect < ylim[1]) {
      new_point <- data.table(x=bottom_intersect, y=ylim[1])
      this_axis <- "bottom"
    } else {
      new_point <- data.table(x=xlim[1], y=left_intersect)
      this_axis <- "left"
    }
  } else if (slope.u95 < 0) {
    if (left_intersect > ylim[2]) {
      new_point <- data.table(x=top_intersect, y=ylim[2])
      this_axis <- "top"
    } else {
      new_point <- data.table(x=xlim[1], y=left_intersect)
      this_axis <- "left"
    }
  }
  poly_dt <- rbind(poly_dt, new_point)

  # (4) From the intercept.u95 to the bottom/left/top border based on the slope.l95
  #     - we also have to add the plot corners if applicable.
  prev_axis <- this_axis
  bottom_intersect <- intersect_hline(intercept.u95, slope.l95, ylim[1])
  left_intersect <- intersect_vline(intercept.u95, slope.l95, xlim[1])
  top_intersect <- intersect_hline(intercept.u95, slope.l95, ylim[2])
  if (slope.l95 >= 0) {
    if (left_intersect < ylim[1]) {
      new_point <- data.table(x=bottom_intersect, y=ylim[1])
      this_axis <- "bottom"
    } else {
      new_point <- data.table(x=xlim[1], y=left_intersect)
      this_axis <- "left"
    }
  } else if (slope.l95 < 0) {
    if (left_intersect > ylim[2]) {
      new_point <- data.table(x=top_intersect, y=ylim[2])
      this_axis <- "top"
    } else {
      new_point <- data.table(x=xlim[1], y=left_intersect)
      this_axis <- "left"
    }
  }
  # add corners if applicable
  if (prev_axis == "bottom" && this_axis == "left") {
    poly_dt <- rbind(poly_dt, bottom_left)
  } else if (prev_axis == "bottom" && this_axis == "top") {
    poly_dt <- rbind(poly_dt, bottom_left, top_left)
  } else if (prev_axis == "left" && this_axis == "top") {
    poly_dt <- rbind(poly_dt, top_left)
  }
  # add new point
  poly_dt <- rbind(poly_dt, new_point)

  # Add starting point 
  poly_dt <- rbind(poly_dt, poly_dt[1])
  
  return(poly_dt)
}

# Create ggplot object showing dose response curve
#
# @param instruments data.table containing information about each instrument
# @param mro_dt data.table returned by mr_tryall().
# @param label logical; should snp labels be shown?
#
gg_dose_response <- function(instruments, mro_dt, label=TRUE,
                             point.size=2, bar.width=0.5,
                             point.color="black", bar.color="black") {
	lines_dt <- dose_response_fits(mro_dt)
	ci95_dt <- dose_response_ci95(mro_dt)
	plot_lim <- dose_response_plot_lim(instruments, ci95_dt)

  points_dt <- data.table(
    exposure.beta = instruments$Prot.Effect.pQTL,
		exposure.lSE = instruments$Prot.Effect.pQTL - instruments$Prot.SE.pQTL,
		exposure.uSE = instruments$Prot.Effect.pQTL + instruments$Prot.SE.pQTL,
    outcome.beta = instruments$effect.gwas,
		outcome.lSE = instruments$effect.gwas - instruments$se.gwas,
		outcome.uSE = instruments$effect.gwas + instruments$se.gwas,
    variants = instruments$ALT_ID,
    MAF = instruments$EAF.pqtl,
    type = instruments$type)

  # Construct base plot elements 
  g <- ggplot(points_dt) + 
       # Aesthetic mapping used for showing instrument effects on exposure/outcome
			 aes(x = exposure.beta, xmin = exposure.lSE, xmax = exposure.uSE,
           y = outcome.beta, ymin = outcome.lSE, ymax = outcome.uSE) +
       # plot limits 
       scale_x_continuous(limits = plot_lim$xlim, expand = c(0,0)) +
       scale_y_continuous(limits = plot_lim$ylim, expand = c(0,0)) +
       # lines showing x and y axis intercepts
			 geom_hline(yintercept=0, color="#525252", linetype="dotted") +
       geom_vline(xintercept=0, color="#525252", linetype="dotted")
 
  # Add polygons showing dose-response curve uncertainty (if applicable)
  if (nrow(ci95_dt) > 0) {
    # First need to define polygon using function above
		ci95_poly <- ci95_dt[!is.na(intercept.l95) & !is.na(intercept.u95) &
										 !is.na(slope.l95) & !is.na(slope.u95),
										 line_ci95_poly(intercept.l95, intercept.u95,
																		slope.l95, slope.u95,
																		plot_lim$xlim, plot_lim$ylim),
										 by=.(method)]

    # A path is added around the polygon separately so we can control the 
    # alpha level of the path.
	  g <- g + geom_polygon(data=ci95_poly, inherit.aes=FALSE, show.legend=FALSE,
													aes(x=x, y=y, fill=method), alpha=0.1) +
						 geom_path(data=ci95_poly, inherit.aes=FALSE, show.legend=FALSE,
											 aes(x=x, y=y, color=method), size=0.05, alpha=0.4)
  }
  # Add dose-response lines of best fit (if applicable)
  if (nrow(lines_dt) > 0) {
    g <- g + geom_abline(data=lines_dt, linetype=2, size=0.5,
              					 aes(intercept=intercept, slope=slope, color=method))
  }
   Add dose response fills and colors. Ignored if above not on plot.
   color_map <- c("Simple median" = "#004529", 
		  						 "Weighted median" = "#238443",
                  "Penalized weighted median" = "#78c679",
                  "IVW" = "#810f7c", 
                  "Penalized IVW" = "#88419d",
                  "Robust IVW" = "#8c6bb1", 
                  "Penalized robust IVW" = "#8c96c6",
                  "MR-Egger" = "#e31a1c", 
                  "Penalized MR-Egger" = "#fc4e2a",
                  "Robust MR-Egger" = "#fd8d3c", 
                  "Penalized robust MR-Egger" = "#feb24c")
  g <- g + scale_color_manual(name = "Causal estimate (95% CI)", drop = TRUE, values = color_map) +
           scale_fill_manual(name = "Causal estimate (95% CI)", drop = TRUE, values = color_map)

  # Show instrument betas on exposure and outcome (and SEs)
  g <- g + geom_errorbarh(height=0, color=bar.color, size=bar.width) + 
           geom_errorbar(width=0, color=bar.color, size=bar.width) + 
           geom_point(data=points_dt[type == "trans"], shape=21, size=point.size, color="black", fill="black") +
           geom_point(data=points_dt[type == "cis"], shape=21, size=point.size, color="black", fill="#feb24c") + 

  # Add labels to instruments
  if (label) {
		g <- g + geom_label_repel(aes(label = instruments),
															size = 2, force = 10, min.segment.length = 0, segment.size=0.3,
															nudge_x = plot_lim$pw * 0.01, nudge_y = -plot_lim$ph * 0.01)
  }

  # Add theme parameters (these can be overridden after return by adding new theme() element)
  g <- g + theme_bw() + 
					 theme(panel.grid=element_blank(),
             axis.title=element_text(size = 10),
             axis.text=element_text(size = 8),
             legend.title=element_text(size = 10, face = "bold"),
             legend.text=element_text(size = 8))
  return(g)
}
