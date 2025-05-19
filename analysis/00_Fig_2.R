###################################################################################################
#' Color gradient for Figure 2
#'
#'Produces:
#'        
#'        - outputs/Fig_2_gradient.tiff
#'        
#' @author Nicolas Mouquet, \email{nicolas.mouquet@@cnrs.fr},
#' @date 2025/04/10 first created
##################################################################################################


library(ggplot2)

n <- 300
df <- expand.grid(x = seq(0, 1, length.out = n),
                  y = seq(0, 1, length.out = n))

df$z <- (df$x + (1 - df$y)) / 2 

gradient_colors <- colorRampPalette(c("#8DD694", "#FFF2CC", "#FF5E3A"))(100)

a <- ggplot(df, aes(x = x, y = y, fill = z)) +
  geom_tile() +
  scale_fill_gradientn(colors = gradient_colors) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "none")

ggsave(file=here::here("outputs","Fig_2_gradient.tiff"), a,width = 15, height = 15, dpi = 300, units = "cm", device='tiff') 

