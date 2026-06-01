######################################################
###plot
library(ggplot2)
library(patchwork)
library(RColorBrewer)


methods <- c("HGCGMA", "NPN-Glasso", "Glasso", "MAEGG","Homo-GCGMA")
Ns <- c(150, 300, 500, 800)
ps <- c(30, 50, 100)

df <- expand.grid(N = Ns, p = ps, method = methods)
df$mean <- NA_real_
df$sd   <- NA_real_

for (p in ps) {
  for (N in Ns) {
    var_mean <- paste0("KL_ave_N", N, "_p", p)
    var_sd   <- paste0("KL_sd_N", N, "_p", p)
    if (exists(var_mean) && exists(var_sd)) {
      mean_vec <- get(var_mean)
      sd_vec   <- get(var_sd)
      idx <- df$N == N & df$p == p
      df$mean[idx] <- mean_vec
      df$sd[idx]   <- sd_vec
    } else {
      warning("var ", var_mean, " or ", var_sd, " not exist")
    }
  }
}
df <- na.omit(df)
df$method <- factor(df$method, levels = methods)


cols <- brewer.pal(5, "Set1")
names(cols) <- methods
shade_cols <- scales::alpha(cols, 0.25)


auto_ylim <- list()
for (p_val in ps) {
  d <- subset(df, p == p_val)
  max_val <- max(d$mean + d$sd, na.rm = TRUE)
  min_val <- min(d$mean - d$sd, na.rm = TRUE)
  auto_ylim[[as.character(p_val)]] <- c(min_val, max_val)
}


adjusted_upper_30 <- auto_ylim[["30"]][2] + 0.5

plots <- list()
for (p_val in ps) {
  d <- subset(df, p == p_val)
  
  
  if (p_val == 30) {
    y_upper <- auto_ylim[["30"]][2]+2
    y_lower <- min(d$mean - d$sd, na.rm = TRUE)  
  } else if (p_val == 50) {
    
    y_upper <- auto_ylim[["50"]][2]+1
    y_lower <- auto_ylim[["50"]][1]
  } else if (p_val == 100) {
    
    y_upper <- auto_ylim[["100"]][2]
    y_lower <- auto_ylim[["100"]][1]
  }
  
  p <- ggplot(d, aes(x = N, y = mean, color = method, fill = method, group = method)) +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.25, linetype = 0) +
    geom_line(size = 1) +
    geom_point(aes(shape = method), size = 3) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = shade_cols) +
    scale_shape_manual(values = c(16, 17, 15, 18,8)) +
    scale_x_continuous(breaks = Ns) +
    coord_cartesian(ylim = c(y_lower, y_upper)) +   
    labs(title = paste0("p = ", p_val), x = "N", y = NULL) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          legend.position = c(1, 1),
          legend.justification = c("right", "top"),
          legend.background = element_rect(fill = "white", color = "black", size = 0.2),
          legend.title = element_blank(),
          legend.key = element_blank())
  plots[[as.character(p_val)]] <- p
}


combined_plot <- wrap_plots(plots, nrow = 1)
print(combined_plot)
#ggsave("KL_non_Gaussian.jpeg", combined_plot, width = 12, height = 4, dpi = 300)

