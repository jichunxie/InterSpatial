#' Compute Nearest Outgoing Transport Distance (NOTD)
#'
#' For each receiver meta-cell, finds the minimum Wasserstein distance from any source meta-cell.
#'
#' @param meta_cells_distance A list returned by `compute_wasserstein_distance()` or `compute_wasserstein_distance_between_cell_type()`, containing a data frame named `wasserstein_distance`.
#'
#' @return A data frame with one row per receiver meta-cell, containing:
#' \itemize{
#'   \item \code{clust_ind_receptor}: Receiver meta-cell ID
#'   \item \code{min_wass_dist}: Minimum Wasserstein distance from any source meta-cell
#' }
#'
#' @export

compute_NOTD <- function(meta_cells_distance)
{
  results <- meta_cells_distance$wasserstein_distance
  NOTD_df <- results %>%
    group_by(clust_ind_receptor) %>%
    summarise(min_wass_dist = min(wass_dist, na.rm = TRUE), .groups = "drop")
  as.data.frame(NOTD_df)
}



#' Compute Kendall's Tau Between NOTD from Source Meta-cells and Receptor Expression in Receiving Meta-cells
#'
#' Computes Kendall's correlation between NOTD (minimum Wasserstein distance to source meta-cells)
#' and expression of a specified receptor gene across receiver meta-cells.
#'
#' @param meta_cells_distance A list returned by `compute_wasserstein_distance()` or `compute_wasserstein_distance_between_cell_type()`, containing \code{wasserstein_distance}, \code{source_name}, and \code{receiver_name}.
#' @param gene_count_matrix_meta_cells A matrix of gene expression values where rows are meta-cells and columns are genes.
#' @param receptor_gene A character string indicating the name of the receptor gene to test.
#'
#' @return A named list with the following entries:
#' \itemize{
#'   \item \code{tau}: Raw Kendall's tau statistic
#'   \item \code{tau_round}: Rounded tau (2 decimals)
#'   \item \code{p_value}: Raw p-value
#'   \item \code{P_val}: Formatted string for plotting (e.g. "p-value < 0.001")
#' }
#'
#' @export

compute_tau <- function(meta_cells_distance, gene_count_matrix_meta_cells, receptor_gene) {
  NOTD_df <- compute_NOTD(meta_cells_distance)
  NOTD_df$receptor_val <- gene_count_matrix_meta_cells[NOTD_df$clust_ind_receptor, receptor_gene]

  plot_data <- data.frame(
    distance = NOTD_df$min_wass_dist,
    receptor_val = NOTD_df$receptor_val
  )

  corr_val <- cor(plot_data$distance, plot_data$receptor_val, method = "kendall")
  corr_val_round <- round(corr_val, 2)
  test_cor <- cor.test(plot_data$distance, plot_data$receptor_val, method = "kendall")
  p_val <- test_cor$p.value
  formatted_p_val <- ifelse(p_val < 0.001, "p-value < 0.001", paste0("p-value = ", round(p_val, 3)))

  return (list("tau" = corr_val, "tau_round" = corr_val_round, "p_value" = p_val, "P_val" = formatted_p_val))
}


#' Plot Receptor Expression vs. NOTD With Piecewise Linear Fit
#'
#' Generates a scatterplot of receptor gene expression versus NOTD (minimum Wasserstein distance),
#' fitted with a piecewise linear regression and confidence band.
#'
#' @param meta_cells_distance A list returned by `compute_wasserstein_distance()` or `compute_wasserstein_distance_between_cell_type()`, containing \code{wasserstein_distance}, \code{source_name}, and \code{receiver_name}.
#' @param gene_count_matrix_meta_cells A matrix of gene expression values where rows are meta-cells and columns are genes.
#' @param receptor_gene Character string specifying the receptor gene to plot.
#' @param breakpoints Optional numeric vector of breakpoints for piecewise linear regression. If \code{NULL}, uses the 25\% and 50\% quantiles of NOTD as defaults.
#' @param ggtitle_custom Optional character string to override the default ggplot title.
#'
#' @return A \code{ggplot2} object showing the scatter plot with regression line and confidence interval.
#'
#' @importFrom ggplot2 ggplot geom_point geom_line geom_ribbon xlab ylab ggtitle theme annotate aes element_text element_rect element_line scale_y_continuous
#' @importFrom stats cor cor.test lm predict
#' @export

plot_receptor_expression <- function(meta_cells_distance,
                                     gene_count_matrix_meta_cells,
                                     receptor_gene,
                                     breakpoints = NULL,
                                     ggtitle_custom = NULL) {
  NOTD_df <- compute_NOTD(meta_cells_distance)
  NOTD_df$receptor_val <- gene_count_matrix_meta_cells[NOTD_df$clust_ind_receptor, receptor_gene]

  plot_data <- data.frame(
    distance = NOTD_df$min_wass_dist,
    receptor_val = NOTD_df$receptor_val
  )

  corr_val <- cor(plot_data$distance, plot_data$receptor_val, method = "kendall")
  corr_val_round <- round(corr_val, 2)
  test_cor <- cor.test(plot_data$distance, plot_data$receptor_val, method = "kendall")
  p_val <- test_cor$p.value
  formatted_p_val <- ifelse(p_val < 0.001, "p-value < 0.001", paste0("p-value = ", round(p_val, 3)))

  # Set breakpoints if not provided
  if (is.null(breakpoints)) {
    breakpoints <- quantile(plot_data$distance, probs = c(0.25, 0.5, 0.75))
    fit <- lm(receptor_val ~ distance + I((distance - breakpoints[1]) * (distance >= breakpoints[1]))+
                I((distance - breakpoints[2]) * (distance >= breakpoints[2])),
              data = plot_data)
  }else{
    if(breakpoints=="Fit 1")
    {
      breakpoints <- quantile(plot_data$distance, probs = c(0.25, 0.5, 0.75))
      fit <- lm(receptor_val ~ distance +
                  I((distance - breakpoints[2]) * (distance >= breakpoints[2])),
                data = plot_data)
    }
    # Fit piecewise linear model

  }



  # Predict with confidence intervals
  pred_with_conf <- predict(fit, interval = "confidence")
  plot_data$pred <- pred_with_conf[, "fit"]
  plot_data$lwr <- pred_with_conf[, "lwr"]
  plot_data$upr <- pred_with_conf[, "upr"]

  source_name <- meta_cells_distance$source_name
  receiver_name <- meta_cells_distance$receiver_name

  title_text <- if (!is.null(ggtitle_custom)) ggtitle_custom else paste0(receptor_gene, " Receptor in : ", receiver_name)

  pl <- ggplot(plot_data, aes(distance, receptor_val)) +
    geom_point(size = 3, color = "blue") +
    geom_line(aes(y = pred), size = 3, col = "red") +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4, fill = "#979696") +
    xlab(paste0("NOTD from ", source_name)) + ylab("Expression") +
    ggtitle(title_text) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      plot.title = element_text(size = 35, face = 'bold'),
      legend.title = element_text(size = 20, face = 'bold'),
      legend.text = element_text(size = 20),
      legend.position = "none",
      legend.key.width = unit(2, "cm"),
      legend.spacing.x = unit(1.0, 'cm'),
      axis.title = element_text(face = "bold", size = 35),
      axis.text = element_text(size = 30)
    ) +
    annotate("text", x = Inf, y = Inf,
             label = bquote(tau == .(corr_val_round)),
             size = 15, color = "red", hjust = 1, vjust = 1) +
    annotate("text", x = Inf, y = Inf,
             label = formatted_p_val,
             size = 15, color = "black", hjust = 1, vjust = 2.5)

  return(pl)
}
