#### ---- PLOT FUNCTIONS ---- ####

# standard theme for umap
theme_umap = function(...){
  theme_bw() +
    theme(
      ...,
      aspect.ratio = 3/3,
      panel.grid = element_blank(),
      strip.background = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank()
    ) 
}

fix_legend = function(...){
  guides(
    color = guide_legend(override.aes = list(size = 2)),
    nrow = 2
  )
}