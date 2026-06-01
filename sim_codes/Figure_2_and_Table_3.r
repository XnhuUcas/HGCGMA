##########
##plot
# install.packages(c("igraph", "ggplot2", "patchwork", "ggnetwork"))
rm(list = ls())
library(MASS)
library(igraph)
library(ggplot2)
library(ggnetwork)
library(patchwork)
library(Matrix)
library(kernlab)
source("candidate.r")
source("cov_est.r")
source("Omega_est.r")
source("estimate_Omega.r")
source("algorithms.r")
source("data_generate.r")
p=dim(Omega_est_method[,,1])[1]
Omega_positive=array(0,dim = c(p,p,5))
for (k in 1:5) {
  Omega=ensure_positive_definite(Omega_est_method[,,k])
  Omega_positive[,,k]=cov2cor(Omega)
}

library(org.Hs.eg.db)
original_ensembl <- colnames(brain_data_used[[1]])
hgnc_symbols <- mapIds(
  org.Hs.eg.db,
  keys = original_ensembl,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first" 
)
hgnc_symbols[is.na(hgnc_symbols)] <- original_ensembl[is.na(hgnc_symbols)]

dimnames(Omega_positive) <- list(hgnc_symbols, hgnc_symbols, NULL)

gene_100=dimnames(Omega_positive)[[1]]

plot_top50_compact_concentric <- function(mat, title_name, target_genes = NULL, layout_coords = NULL) {
  
  g_full <- graph_from_adjacency_matrix(mat, mode = "undirected", diag = FALSE, weighted = TRUE)
  V(g_full)$label_name <- rownames(mat)
  
  if (is.null(target_genes)) {
    full_degrees <- degree(g_full)
    top_nodes <- order(full_degrees, decreasing = TRUE)[1:min(50, vcount(g_full))]
    target_genes <- V(g_full)$label_name[top_nodes]
  }
  
  
  mat_sub <- mat[target_genes, target_genes, drop = FALSE]
  g <- graph_from_adjacency_matrix(mat_sub, mode = "undirected", diag = FALSE, weighted = TRUE)
  V(g)$label_name <- rownames(mat_sub)
  
  
  node_degrees <- degree(g)
  V(g)$degree <- node_degrees
  
  
  if (is.null(layout_coords)) {
    breaks <- quantile(node_degrees, probs = c(0, 0.5, 0.85, 1))
    if(length(unique(breaks)) < 4) {
      shells <- cut(rank(node_degrees, ties.method = "first"), breaks = 3, labels = FALSE)
    } else {
      shells <- cut(node_degrees, breaks = breaks, include.lowest = TRUE, labels = FALSE)
    }
    
    coords <- matrix(0, nrow = vcount(g), ncol = 2)
    for (layer in 1:3) {
      idx <- which(shells == layer)
      if (length(idx) > 0) {
        radius <- ifelse(layer == 1, 2.0, ifelse(layer == 2, 1.1, 0.4))
        angles <- seq(0, 2 * pi, length.out = length(idx) + 1)[1:length(idx)]
        coords[idx, 1] <- radius * cos(angles)
        coords[idx, 2] <- radius * sin(angles)
      }
    }
    layout_coords <- coords
  }
  
  if (ecount(g) == 0) {
    g_net <- ggnetwork(g, layout = layout_coords)
  } else {
    E(g)$abs_weight <- abs(E(g)$weight)
    g_net <- ggnetwork(g, layout = layout_coords)
  }
  
  if ("name" %in% colnames(g_net)) {
    g_net$gene_label <- g_net$name
  } else if ("vertex.names" %in% colnames(g_net)) {
    g_net$gene_label <- g_net$vertex.names
  } else {
    g_net$gene_label <- rownames(mat_sub)[g_net$vertex.id]
  }
  
  p_net <- ggplot(g_net, aes(x = x, y = y, xend = xend, yend = yend)) +
    
    geom_edges(aes(linewidth = abs_weight), color = "grey30", alpha = 0.5, show.legend = TRUE) +
    
    geom_nodes(aes(size = degree, color = degree), show.legend = TRUE) +
    
    geom_nodetext_repel(aes(label = gene_label), fontface = "bold", size = 2.8, color = "black") +
    
    scale_linewidth_continuous(name = "Edge Weight", range = c(0.2, 1.8)) +
    scale_size_continuous(name = "Degree", range = c(3, 8)) + 
    scale_color_gradient(name = "Degree", low = "cyan2", high = "red") + 
    
   
  guides(
    color = guide_colourbar(title = "Degree", barwidth = 5, barheight = 0.5), 
    size = "none", 
    linewidth = guide_legend(title = "Edge Weight", override.aes = list(color = "grey30")) 
  ) +
    
    labs(title = title_name) +
    theme_blank() + 
    
    
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    panel.border = element_blank(),      
    strip.background = element_blank(),   
    
    
    legend.position = "bottom",          
    legend.direction = "horizontal",     
    legend.box = "horizontal",           
    legend.background = element_blank(),  
    legend.box.background = element_blank(), 
    legend.key = element_blank(),         
    
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 7),
    legend.spacing.x = unit(0.5, "cm")    
  )
  
  return(list(plot = p_net, genes = target_genes, layout = layout_coords))
}

base_result <- plot_top50_compact_concentric(Omega_positive[,,1], "HGCGMA")
fixed_genes <- base_result$genes
fixed_layout <- base_result$layout

res=plot_top50_compact_concentric(Omega_positive[,,1], " ", target_genes = fixed_genes, layout_coords = fixed_layout)

#ggsave("HGCGMA.pdf",plot = res$plot,width = 6, height = 6, dpi = 1200)
################################
g_base <- graph_from_adjacency_matrix(Omega_positive[,,1], mode = "undirected", diag = FALSE, weighted = TRUE)

all_degrees <- degree(g_base)

top50_genes_with_degrees <- data.frame(
  Gene_Name = base_result$genes,
  Degree = all_degrees[base_result$genes]
)

top50_genes_with_degrees <- top50_genes_with_degrees[order(top50_genes_with_degrees$Degree, decreasing = TRUE), ]


