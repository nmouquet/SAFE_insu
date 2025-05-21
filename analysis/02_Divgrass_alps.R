###################################################################################################
#' Spatial Insurance, analysis for alpine herbaceous communities
#'
#'Uses:
#'
#'      - upload the DIVGRASS dataset data/Gauzere_2023/divgrass_coms.Rdata
#'      - upload the DIVGRASS occurences data/Gauzere_2023/20200602_Divgrass.rds
#'      - upload the BiogeoRegions2016 shapefile (data/eea_v_3035_1_mio_biogeo-regions_p_2016_v01_r00/BiogeoRegions2016.shp)
#'      - upload the DIVGRASS trait data data/Gauzere_2023/divgrass_filled_species_traits.Rdata
#'
#'      Gaüzère P. et al. The functional trait distinctiveness of plant species is scale dependent.
#'      Ecography 2023, e06504 (2023).
#'
#'Produces:
#'
#'      - results/Div_grass_french_alps.Rdata: the subset of divgrass for the french alps
#'      - results/alpine_within_radius.Rdata: neighbor points within each buffer of the Div_grass_french_alps
#'      - results/alpine_buffer.Rdata: coordinates of neighbors points
#'      - results/disttrait_divgrass.Rdata: functional distance matrix
#'      - results/pcoas_divgrass.Rdata: pcoa on the functional distance matrix
#'      - results/divgrass_insurance.RData: insurance for all the plots in the alps
#'      - outputs/S2_Fig_1.tiff
#'      - outputs/S2_Fig_2.tiff
#'      - outputs/Fig_4b_left.tiff
#'      - outputs/Fig_4c.tiff
#'      - outputs/Fig_4d.tiff
#'
#'      Note : Fig_4a is produced with tmap (lines 535-548) and pasted into the final figure
#'
#' @author Nicolas Mouquet, \email{nicolas.mouquet@@cnrs.fr},
#' @date 2024/05/10 first created, major update 2025/03/01
##################################################################################################

rm(list = ls())

library(ggplot2)
library(dplyr)
crs_lambert93 <- 2154

####Function Insu----

Insu <- function(insured, insurer, occ_list, dist_mat, D_insu, D_thr) {
  # insured <- buffer_cell_used$neighbor[j]
  # insurer <-focal_cell
  # occ_list <- occ_list
  # dist_mat <- dist_mat
  # D_insu <- D_insu
  # D_thr <- D_thr

  spe_insured <- toupper(occ_list[[as.character(insured)]])
  spe_insurer <- toupper(occ_list[[as.character(insurer)]])

  #Distance and insurance
  dist_insured <- dist_mat[
    colnames(dist_mat) %in% spe_insured,
    rownames(dist_mat) %in% spe_insured
  ]
  insured_Di = funrar::distinctiveness_global(dist_insured)
  colnames(insured_Di) <- c("Species", "Di")

  QD <- as.numeric(quantile(insured_Di$Di, probs = seq(0, 1, 0.01))[D_thr + 1])

  if (D_thr == 0) {
    insured_D = insured_Di$Species
  } else {
    insured_D = insured_Di$Species[insured_Di$Di >= QD]
  }

  Insu <- do.call(
    rbind,
    lapply(insured_D, function(sp_insured) {
      #sp_insured <- insured_D[1]

      #is the distinct species from the focal community is insured in the neighboring com ?
      #threshold for insurance radius

      if (sp_insured %in% spe_insurer) {
        dist_insurer <- dist_mat[
          colnames(dist_mat) %in% spe_insurer,
          rownames(dist_mat) %in% spe_insurer
        ]
        Insurance <- sum(dist_insurer[sp_insured, ] < D_insu)
      } else {
        dist_insurer <- dist_mat[
          colnames(dist_mat) %in% c(spe_insurer, sp_insured),
          rownames(dist_mat) %in% c(spe_insurer, sp_insured)
        ]
        Insurance <- sum(dist_insurer[sp_insured, ] < D_insu) - 1
      }

      cbind.data.frame(Species = sp_insured, Insurance = Insurance)
    })
  )

  Inward <- round(100 * sum(Insu$Insurance > 0) / nrow(Insu), 2)

  spnotinsu = paste0(Insu$Species[Insu$Insurance == 0], collapse = ",")
  spinsu = paste0(Insu$Species[Insu$Insurance > 0], collapse = ",")

  rm(dist_insured)
  rm(Insu)

  cbind.data.frame(
    D_insu = D_insu,
    D_thr = D_thr,
    insured = insured,
    insurer = insurer,
    Div_insured = length(spe_insured),
    Div_insurer = length(spe_insurer),
    Inward = Inward,
    Div_distinct = length(insured_D),
    spnotinsu = spnotinsu,
    spinsu = spinsu
  )
}

Insu_who <- function(insured, insurer, occ_list, dist_mat, D_insu, D_thr) {
  # insured <- focal_cell
  # insurer <-buffer_cell_used$neighbor[2]
  # occ_list <- occ_list
  # dist_mat <- disttrait
  # D_insu <- D_insu
  # D_thr <- D_thr

  spe_insured <- toupper(occ_list[[as.character(insured)]])
  spe_insurer <- toupper(occ_list[[as.character(insurer)]])

  #Distance and insurance
  dist_insured <- dist_mat[
    colnames(dist_mat) %in% spe_insured,
    rownames(dist_mat) %in% spe_insured
  ]
  insured_Di = funrar::distinctiveness_global(dist_insured)
  colnames(insured_Di) <- c("Species", "Di")

  QD <- as.numeric(quantile(insured_Di$Di, probs = seq(0, 1, 0.01))[D_thr + 1])

  if (D_thr == 0) {
    insured_D = insured_Di$Species
  } else {
    insured_D = insured_Di$Species[insured_Di$Di >= QD]
  }

  Insu <- do.call(
    rbind,
    lapply(insured_D, function(sp_insured) {
      #sp_insured <- insured_spe_Di80[1]

      #is the insured sp in the insurerboring com ?
      #threshold for insurance radius

      if (sp_insured %in% spe_insurer) {
        dist_insurer <- dist_mat[
          colnames(dist_mat) %in% spe_insurer,
          rownames(dist_mat) %in% spe_insurer
        ]
        Insurance <- sum(dist_insurer[sp_insured, ] < D_insu)
      } else {
        dist_insurer <- dist_mat[
          colnames(dist_mat) %in% c(spe_insurer, sp_insured),
          rownames(dist_mat) %in% c(spe_insurer, sp_insured)
        ]
        Insurance <- sum(dist_insurer[sp_insured, ] < D_insu) - 1
      }

      cbind.data.frame(Species = sp_insured, Insurance = Insurance)
    })
  )

  colnames(Insu) <- c("Species", insurer)

  Insu
}

#----

####LOAD POINTS DATA & SHAPE_FILE & SUBSET FRENCH ALPS----

# Load shapefile and convert to Lambert-93 immediately
shp_file <- sf::st_read(here::here(
  "data",
  "eea_v_3035_1_mio_biogeo-regions_p_2016_v01_r00",
  "BiogeoRegions2016.shp"
))
shp_file <- sf::st_transform(shp_file, crs_lambert93)

# Load points data
load(here::here("data/Gauzere_2023/divgrass_coms.Rdata"))

# Remove NAs in coordinates and convert to sf object (already in Lambert-93)
points_sf <- divgrass_com[
  !is.na(divgrass_com$X_L93) & !is.na(divgrass_com$Y_L93),
]
points_sf <- sf::st_as_sf(
  points_sf,
  coords = c("X_L93", "Y_L93"),
  crs = crs_lambert93
)

# Remove duplicates based on X_WGS84 and Y_WGS84 and keep only the first occurence for duplicates
points_sf_unique <- points_sf %>%
  dplyr::group_by(X_WGS84, Y_WGS84) %>%
  dplyr::slice(1) %>% # Keep only the first occurrence
  dplyr::ungroup() %>%
  sf::st_as_sf() # Ensure it remains an sf object

# Define bounding box for the French Alps (approximate) and transform to Lambert-93
bbox_french_alps <- sf::st_bbox(
  c(xmin = 4.95, xmax = 7.85, ymin = 43.85, ymax = 46.5),
  crs = 4326
) %>%
  sf::st_as_sfc() %>%
  sf::st_transform(crs_lambert93)

# Extract only the alpine biome from the shapefile
alpine_regions <- shp_file[shp_file$short_name == "alpine", ]

# Ensure French Alps geometry is clean and merged
french_alps <- sf::st_intersection(alpine_regions, bbox_french_alps) %>%
  sf::st_make_valid() %>% # Fix invalid geometries
  sf::st_union() %>% # Merge overlapping parts
  sf::st_cast("MULTIPOLYGON") # Ensure proper geometry type

# Filter points inside the French Alps Box
Div_grass_french_alps <- points_sf_unique[
  sf::st_within(points_sf_unique, french_alps, sparse = FALSE),
]

# Remove points that are in Italy
sites_italy <- c(
  'SOP_T_1937.3.1',
  'SOP_T_1937.4.1',
  'SOP_T_1937.5.1',
  'SOP_T_2303.5.1',
  'SOP_T_862.2.7',
  'SOP_T_491.2.6',
  'SOP_T_2303.4.1',
  'SOP_T_491.8.5',
  'SOP_T_1891.24.4',
  'SOP_T_1194.3.20',
  'SOP_T_194.7.19',
  'SOP_T_194.7.20',
  'SOP_T_1194.4.9',
  'CME1_C16041898',
  'CME1_C16046914'
)

Div_grass_french_alps <- Div_grass_french_alps[
  !Div_grass_french_alps$com %in% sites_italy,
]

# Keep only points with gps precision >1137 meters
Div_grass_french_alps <- Div_grass_french_alps[
  Div_grass_french_alps$PREC_LOC %in% c('704', '1137'),
]

#Map
# tmap::tmap_mode("view")  # Enable interactive mode
#
# tmap::tm_shape(french_alps) +
#   tmap::tm_polygons(col = "lightblue", border.col = "black") +
#   tmap::tm_shape(Div_grass_french_alps) +
#   tmap::tm_dots(fill = "red", size = 0.2)
#
#Save & Load
save(
  Div_grass_french_alps,
  file = here::here("results", "Div_grass_french_alps.Rdata")
)
load(here::here("results", "Div_grass_french_alps.Rdata"))

#----

###POINTS WITHIN BUFFER----
#Long, run once (and load)

# Define radius in meters
radius <- 3000

# Find neighbors within radius for each point in Div_grass_french_alps (this is long)
alpine_within_radius <- sf::st_is_within_distance(
  Div_grass_french_alps,
  Div_grass_french_alps,
  dist = radius
)
save(
  alpine_within_radius,
  file = here::here("results", "alpine_within_radius.Rdata")
)

# Create a data frame showing which points are neighbors
alpine_buffer <- data.frame(
  focal_cell = rep(
    Div_grass_french_alps$com,
    sapply(alpine_within_radius, length)
  ), # Use 'com' as ID
  neighbor_cell = Div_grass_french_alps$com[unlist(alpine_within_radius)] # Get 'com' from indices
)

# Join back to get coordinates of neighbors (this can be long too)
alpine_buffer <- merge(
  alpine_buffer,
  cbind(neighbor_cell = Div_grass_french_alps$com, Div_grass_french_alps),
  by = "neighbor_cell"
)

#Save
save(alpine_buffer, file = here::here("results", "alpine_buffer.Rdata"))

#Load
load(file = here::here("results", "alpine_buffer.Rdata"))
load(file = here::here("results", "alpine_within_radius.Rdata"))

#Map example
# Select a point
focal_com <- "CAL1_848810"
id <- which(Div_grass_french_alps$com %in% focal_com)
reference_point <- Div_grass_french_alps[id, ]

# Select its neighbors
neighbors <- Div_grass_french_alps[alpine_within_radius[[id]], ]

# Create the buffer
buffer_Xkm <- sf::st_buffer(reference_point, dist = 5000)
buffer_Xkm_wgs84 <- sf::st_transform(buffer_Xkm, crs = 4326)

#map
tmap::tmap_mode("view") # Enable interactive mode

tmap::tm_shape(french_alps) +
  tmap::tm_polygons(col = "lightblue", border.col = "black") +
  tmap::tm_shape(Div_grass_french_alps) +
  tmap::tm_shape(buffer_Xkm_wgs84) + # Add the buffer as a polygon
  tmap::tm_polygons(col = "green", fill_alpha = 0.3, border.col = "darkgreen") + # Customize buffer appearance
  tmap::tm_shape(Div_grass_french_alps) +
  tmap::tm_dots(fill = "gray", size = 0.5) +
  tmap::tm_shape(reference_point) +
  tmap::tm_dots(fill = "red", size = 0.5) + # Plot the selected reference point
  tmap::tm_shape(neighbors) +
  tmap::tm_dots(fill = "blue", size = 0.4) # Plot its neighbors within 30 km
#----

####TRAITS #1----
#Select the traits
load(here::here('data', 'Gauzere_2023', 'divgrass_filled_species_traits.Rdata'))

#Keep only SLA", "PH", "SM", "LA", "LDMC (see Gauzere et al. 2023)
filled_species_traits_wide <- filled_species_traits_wide[, c(
  "TAXREF_SHORT",
  "SLA",
  "PH",
  "SM",
  "LA",
  "LDMC"
)]

transformed_species_traits <- na.omit(data.frame(
  filled_species_traits_wide$TAXREF_SHORT,
  log10(filled_species_traits_wide[, c("SLA", "PH", "SM", "LA", "LDMC")])
))
colnames(transformed_species_traits)[1] <- "taxa.name"
rownames(transformed_species_traits) <- transformed_species_traits$taxa.name

#remove POTAMOGETON PECTINATUS which is an aquatic plant
transformed_species_traits <- transformed_species_traits[
  !transformed_species_traits$taxa.name %in% "POTAMOGETON PECTINATUS",
]

#----

#OCCURENCES----
#need to compute rather than save (too long)
divgrass <- readRDS(here::here("data/Gauzere_2023/20200602_Divgrass.rds"))

coms <- Div_grass_french_alps$com
#keep only the alps
divgrass_alps <- divgrass[divgrass$com %in% coms, ]
#rempove species names NA
divgrass_alps <- divgrass_alps[!is.na(divgrass_alps$TAXREF_SHORT), ]
#remove species without traits
divgrass_alps <- divgrass_alps[
  divgrass_alps$TAXREF_SHORT %in% transformed_species_traits$taxa.name,
]

#compute a list of sites with species composition

all_com_divgrass <- pbmcapply::pbmclapply(
  coms,
  function(id_site) {
    as.character(divgrass_alps$TAXREF_SHORT[divgrass_alps$com %in% id_site])
  },
  mc.cores = parallel::detectCores() - 2
)
names(all_com_divgrass) <- coms

#----

####TRAITS #2 S2_Fig_1----
#Compute the trait distance matrix and distinctiveness on the alp plants

#Keep only species within alps
alps_transformed_species_traits <- transformed_species_traits[
  transformed_species_traits$taxa.name %in% unique(divgrass_alps$TAXREF_SHORT),
]

#Distance matrix and distinctiveness
disttrait_divgrass = funrar::compute_dist_matrix(
  alps_transformed_species_traits[, -1],
  metric = "euclidean"
)
save(
  disttrait_divgrass,
  file = here::here("results", "disttrait_divgrass.Rdata")
)

global_Di = funrar::distinctiveness_global(disttrait_divgrass)
colnames(global_Di) <- c("Species", "global_di")

pcoas_divgrass <- ape::pcoa(disttrait_divgrass) #long to run
save(pcoas_divgrass, file = here::here("results", "pcoas_divgrass.Rdata"))

#load
load(here::here("results", "pcoas_divgrass.Rdata"))
load(here::here("results", "disttrait_divgrass.Rdata"))

# Eigenvalues of the RFS

df <- data.frame(
  "Eigenvalues" = paste0('E', 1:8),
  "Value" = pcoas_divgrass$"values"$"Relative_eig"[1:8]
)

a <- ggplot(data = df, aes(x = Eigenvalues, y = Value)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_bw()

# Visualization of the RFS

pc_axes <- data.frame(pcoas_divgrass$vectors)
colnames(pc_axes) <- gsub("\\.", "", colnames(pc_axes))

b <- ggplot(pc_axes, aes(Axis1, Axis2)) +
  geom_point(colour = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.5) +
  xlab(paste0("PC1 (", round(df$Value[1] * 100, 1), " %)")) +
  ylab(paste0("PC2 (", round(df$Value[2] * 100, 1), " %)")) +
  theme_bw()

S2_Fig_1 <- gridExtra::arrangeGrob(a, b, ncol = 2)
ggsave(
  file = here::here("outputs", "S2_Fig_1.tiff"),
  S2_Fig_1,
  width = 20,
  height = 10,
  dpi = 300,
  units = "cm",
  device = 'tiff'
)

#----

####LOOP OVER MANY CELLS WITH BUFFERS, Fig_4b,d----

sp_cell <- names(all_com_divgrass)
per_buff_all <- 1
occ_list <- all_com_divgrass
dist_mat <- disttrait_divgrass
D_insu <- median(disttrait_divgrass) * 0.3
D_thr <- 90
set.seed(1)

load(file = here::here("results", "alpine_buffer.Rdata"))

#Compute insurance
#long to run
divgrass_insurance <- do.call(
  rbind,
  pbmcapply::pbmclapply(
    sp_cell,
    function(focal_cell) {
      #focal_cell=sp_cell[600]

      buffer_cells <- cbind.data.frame(
        reference = focal_cell,
        neighbor = alpine_buffer$neighbor_cell[
          alpine_buffer$focal_cell %in% focal_cell
        ]
      )

      #remove the focal_cell from the buffers

      buffer_cells <- buffer_cells[
        !buffer_cells$reference == buffer_cells$neighbor,
      ]

      sp_focal <- unique(all_com_divgrass[[focal_cell]])
      sp_focal <- sp_focal[!is.na(sp_focal)]

      set.seed(1)
      if (per_buff_all < 1)
        buffer_cells <- buffer_cells[
          which(
            buffer_cells$neighbor %in%
              c(
                focal_cell,
                sample(
                  buffer_cells$neighbor[-1],
                  round(length(buffer_cells$neighbor) * per_buff_all)
                )
              )
          ),
        ]
      buffer_cell_used <- buffer_cells[
        unlist(lapply(as.character(buffer_cells$neighbor), function(id) {
          length(unique(all_com_divgrass[[id]])) >= 5
        })),
      ]

      #Compute insurance when there is more than 5 species in the focal cell
      if ((length(sp_focal) >= 5) & (nrow(buffer_cell_used) > 2)) {
        I_buffer <- do.call(
          rbind,
          (pbmcapply::pbmclapply(
            2:dim(buffer_cell_used)[1],
            function(j) {
              #j=2

              a <- Insu(
                insured = focal_cell,
                insurer = buffer_cell_used$neighbor[j],
                occ_list = all_com_divgrass,
                dist_mat = dist_mat,
                D_insu,
                D_thr
              )
              b <- Insu(
                insured = buffer_cell_used$neighbor[j],
                insurer = focal_cell,
                occ_list = all_com_divgrass,
                dist_mat = dist_mat,
                D_insu,
                D_thr
              )
              a$Outward <- b$Inward
              a

              #beta_funk <- Betafunk(insured=focal_cell,insurer=buffer_cell_used$neighbor[j],occ_list,traits=traits)
            },
            mc.cores = parallel::detectCores()
          ))
        )

        #allspecies
        all_species <- unique(unlist(lapply(
          buffer_cell_used$neighbor,
          function(id) {
            occ_list[[id]]
          }
        )))

        #who is never insured
        sp_insu = table(unlist(strsplit(
          I_buffer$spinsu[!I_buffer$spinsu == ""],
          ","
        )))
        sp_insu <- names(sp_insu[sp_insu > 0])
        sp_notinsu = table(unlist(strsplit(
          I_buffer$spnotinsu[!I_buffer$spnotinsu == ""],
          ","
        )))
        sp_notinsu <- names(sp_notinsu[sp_notinsu > 0])
        sp_notinsu <- paste0(
          sp_notinsu[!sp_notinsu %in% sp_insu],
          collapse = ","
        )

        #means

        mean_one_buffer <- cbind.data.frame(
          D_insu = I_buffer$D_insu[1],
          D_thr = I_buffer$D_thr[1],
          cell_focal = focal_cell,
          Inward = mean(I_buffer$Inward[-1], na.rm = TRUE),
          Outward = mean(I_buffer$Outward[-1], na.rm = TRUE),
          Nsp_loc = I_buffer$Div_insured[1],
          Div_distinct = I_buffer$Div_distinct[1],
          Nsp_reg = mean(I_buffer$Div_insurer[-1], na.rm = TRUE),
          gamma = length(all_species),
          nbuffer = nrow(buffer_cell_used),
          sp_notinsu = sp_notinsu
        )

        mean_one_buffer
      } else {
        NULL
      }
    },
    mc.cores = parallel::detectCores() - 1
  )
)

mean(divgrass_insurance$nbuffer)
sd(divgrass_insurance$nbuffer)

#Add env values (not be used)
rownames(divgrass_insurance) <- divgrass_insurance$cell_focal
divgrass_com <- data.frame(divgrass_com)
rownames(divgrass_com) <- divgrass_com$com
divgrass_insurance <- merge(
  divgrass_insurance,
  divgrass_com,
  by = 'row.names',
  all = FALSE
)
divgrass_insurance <- sf::st_as_sf(
  divgrass_insurance,
  coords = c("X_L93", "Y_L93"),
  crs = 2154
)

#Compute quantiles for classification in Source, Sink, and Intermediate
low_threshold <- quantile(
  divgrass_insurance$Inward - divgrass_insurance$Outward,
  0.1,
  na.rm = TRUE
)
high_threshold <- quantile(
  divgrass_insurance$Inward - divgrass_insurance$Outward,
  0.9,
  na.rm = TRUE
)

divgrass_insurance$Category <- ifelse(
  divgrass_insurance$Inward - divgrass_insurance$Outward >= high_threshold,
  "Sink",
  ifelse(
    divgrass_insurance$Inward - divgrass_insurance$Outward <= low_threshold,
    "Source",
    "Intermediate"
  )
)

divgrass_insurance$Category <- as.factor(divgrass_insurance$Category)

#Compute the mean and sd for species richness

mean(divgrass_insurance$Nsp_loc)
sd(divgrass_insurance$Nsp_loc)
mean(divgrass_insurance$Div_distinct)
sd(divgrass_insurance$Div_distinct)
max(divgrass_insurance$Div_distinct)

#Save
save(
  divgrass_insurance,
  file = here::here("results", "divgrass_insurance.RData")
)

#Load
divgrass_insurance <- get(load(here::here(
  'results',
  'divgrass_insurance.RData'
)))

#plot inward vs outward Fig_4b_right
color_s <- rev(RColorBrewer::brewer.pal(n = 8, name = "RdBu"))

library(ggplot2)
library(RColorBrewer)
all_cells <- divgrass_insurance[
  order(divgrass_insurance$Nsp_loc, decreasing = F),
]
all_cells <- divgrass_insurance[
  sample(1:nrow(divgrass_insurance), nrow(divgrass_insurance)),
]
all_cells$log10_Nsp_loc = log10(all_cells$Nsp_loc)

Fig_4b_right <- ggplot(
  all_cells,
  aes(x = Inward, y = Outward, color = log10_Nsp_loc)
) +
  geom_point(alpha = 0.8) +
  ylab("Outward insurance") +
  xlab("Inward insurance") +
  scale_colour_gradientn(colours = color_s) +
  #scale_colour_gradientn(colours = brewer.pal(n = 8, name = "YlOrRd"))+
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 10)
  ) +
  geom_abline(
    slope = 1,
    intercept = as.numeric(quantile(
      all_cells$Outward - all_cells$Inward,
      prob = 0.8,
      na.rm = T
    )),
    linetype = "dashed",
    size = 0.5
  ) +
  geom_abline(
    slope = 1,
    intercept = as.numeric(quantile(
      all_cells$Outward - all_cells$Inward,
      prob = 0.2,
      na.rm = T
    )),
    linetype = "dashed",
    size = 0.5
  ) +
  ylim(0, 100) +
  xlim(0, 100)

#find examples

#plot(divgrass_insurance$Inward, divgrass_insurance$Outward, pch = 19, col = "blue", main = "Click on a point to get its ID")
#clicked_points <- identify(divgrass_insurance$Inward, divgrass_insurance$Outward, labels = divgrass_insurance$cell_focal, plot = TRUE)
#print(paste("You clicked on ID(s):", divgrass_insurance$cell_focal[clicked_points]))

#map source and sinks

# Convert to factor for proper mapping
divgrass_insurance$Category <- factor(
  divgrass_insurance$Category,
  levels = c("Source", "Sink", "Intermediate")
)

# Define custom colors for the 3 categories
category_palette <- c(
  "Sink" = "#FA6C67",
  "Intermediate" = "#FFF7BA",
  "Source" = "#52A153"
)

# look at the categories to check Fig_S2_Fig_2
S2_Fig_2 <- ggplot(
  divgrass_insurance,
  aes(x = Inward, y = Outward, color = Category)
) +
  geom_point(alpha = 0.7) + # Scatter plot of points
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") + # 1:1 line
  scale_color_manual(values = category_palette) + # Define colors
  theme_bw() +
  labs(x = "Inward", y = "Outward", color = "Insurance status") +
  theme(legend.position = "right")

ggsave(
  file = here::here("outputs", "S2_Fig_2.tiff"),
  S2_Fig_2,
  width = 13,
  height = 9,
  dpi = 300,
  units = "cm",
  device = 'tiff'
)

sum(divgrass_insurance$Category %in% "Source")
sum(divgrass_insurance$Category %in% "Sink")
sum(divgrass_insurance$Category %in% "Intermediate")

# Create the map for all categories with an example cell
subper = 1
sub_divgrass_insurance <- divgrass_insurance[
  sample(1:nrow(divgrass_insurance), round(nrow(divgrass_insurance) * subper)),
]
sub_divgrass_insurance <- sub_divgrass_insurance %>%
  dplyr::arrange(factor(Category, levels = c("Intermediate", "Sink", "Source")))

tmap::tmap_mode("view") # Mode interactif (ou "plot" pour un affichage statique)
tmap::tm_shape(sub_divgrass_insurance) +
  tmap::tm_symbols(
    col = "Category",
    palette = category_palette, # Use palette directly, matching categories
    size = 0.5
  ) +
  tmap::tm_layout(title = "Source & Sink Classification")

#Create the map for only sources and sinks only Fig_3a

filtered_divgrass_insurance <- sub_divgrass_insurance[
  sub_divgrass_insurance$Category %in% c("Source", "Sink"),
]
filtered_divgrass_insurance$Category <- factor(
  filtered_divgrass_insurance$Category,
  levels = c("Source", "Sink")
)
category_palette <- c("Sink" = "#FA6C67", "Source" = "#52A153")

tmap::tmap_mode("view") # Mode interactif (ou "plot" pour un affichage statique)
tmap::tm_shape(filtered_divgrass_insurance) +
  tmap::tm_symbols(
    col = "Category",
    palette = category_palette, # Use palette directly, matching categories
    size = 0.5
  ) +
  tmap::tm_layout(legend.show = FALSE)

#Box plot categories Fig4bleft

divgrass_insurance$Category <- factor(
  divgrass_insurance$Category,
  levels = c("Source", "Intermediate", "Sink")
)

Fig_4b_left <- ggplot(
  divgrass_insurance,
  aes(x = Category, y = Nsp_loc, fill = Category)
) +
  geom_boxplot(width = .5, outlier.shape = NA) +
  ylab("Species richness") +
  #xlab("Status")+
  scale_fill_manual(
    values = c(
      "Source" = "#7DBD80",
      "Intermediate" = "#FFF7BA",
      "Sink" = "#FA6C67"
    )
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 14),
    legend.position = "none"
  ) +
  ylim(0, 50) +
  scale_x_discrete(labels = c("Source", "Intermediate", "Sink")) +
  ggplot2::theme(
    axis.line.x = ggplot2::element_line(linetype = "blank"),
    axis.ticks.x = ggplot2::element_blank()
  )

Fig_4b <- gridExtra::arrangeGrob(Fig_4b_left, Fig_4b_right, ncol = 2)
ggsave(
  file = here::here("outputs", "Fig_4b.tiff"),
  Fig_4b,
  width = 15,
  height = 7,
  dpi = 300,
  units = "cm",
  device = 'tiff'
)

#What are the species not insured and in how much plots ?

sp_notinsu <- unique(unlist(strsplit(
  divgrass_insurance$sp_notinsu[!divgrass_insurance$sp_notinsu == ""],
  ","
)))

sp_notinsu_tb <- table(unlist(strsplit(
  divgrass_insurance$sp_notinsu[!divgrass_insurance$sp_notinsu == ""],
  ","
)))
sp_notinsu_tb <- data.frame(sp_notinsu_tb)
colnames(sp_notinsu_tb) <- c("Species", "noinsu")

#for these species create a table of species occurrences across plots

sp_notinsu_occ <- sapply(sp_notinsu, function(sp) {
  sum(sapply(all_com_divgrass, function(community) sp %in% community))
})

sp_notinsu_occ <- data.frame(
  Species = names(sp_notinsu_occ),
  occ = as.integer(sp_notinsu_occ)
)

#merge and compute the %
sp_notinsu_tb <- merge(sp_notinsu_tb, sp_notinsu_occ)
sp_notinsu_tb$pernotinsu = (sp_notinsu_tb$noinsu / sp_notinsu_tb$occ) * 100

#what are the global distinct species ?

Di = funrar::distinctiveness_global(disttrait_divgrass)
colnames(Di) <- c("Species", "Di")

QD <- as.numeric(quantile(Di$Di, probs = seq(0, 1, 0.01))[D_thr + 1])

insured_D = Di$Species[Di$Di >= QD]

#plot the pcoa Fig_4d

pc_axes <- data.frame(pcoas_divgrass$vectors)
colnames(pc_axes) <- gsub("\\.", "", colnames(pc_axes))

#pc_axes_sp_notinsu
pc_axes_sp_notinsu <- pc_axes[rownames(pc_axes) %in% sp_notinsu, ]
pc_axes_sp_notinsu$Species <- rownames(pc_axes_sp_notinsu)
pc_axes_sp_notinsu <- merge(pc_axes_sp_notinsu, sp_notinsu_tb)
pc_axes_sp_notinsu$Distance <- sqrt(
  pc_axes_sp_notinsu$Axis1^2 + pc_axes_sp_notinsu$Axis2^2
)
pc_axes_sp_notinsu_radius <- min(pc_axes_sp_notinsu$Distance)
pc_axes_sp_notinsu <- pc_axes_sp_notinsu[
  order(pc_axes_sp_notinsu$pernotinsu, decreasing = F),
]

#pc_axes_distinct
pc_axes_distinct <- pc_axes[rownames(pc_axes) %in% insured_D, ]
pc_axes_distinct$Distance <- sqrt(
  pc_axes_distinct$Axis1^2 + pc_axes_distinct$Axis2^2
)
pc_axes_distinct_radius <- min(pc_axes_distinct$Distance)

#target some species for illustration
focal_species <- c(
  "BUPLEURUM BALDENSE",
  "LATHYRUS LATIFOLIUS",
  "PETASITES HYBRIDUS",
  "DIPHASIASTRUM ALPINUM"
)
focal_coords <- pc_axes[rownames(pc_axes) %in% focal_species, ]

#plot the pcoa

ptsize = 2

Fig_4d <- ggplot(pc_axes, aes(x = Axis1, y = Axis2)) +
  geom_point(col = "black", alpha = 0.1) +
  geom_vline(xintercept = 0, linetype = "dashed", col = "#9E9E9E") +
  geom_hline(yintercept = 0, linetype = "dashed", col = "#9E9E9E") +
  ylim(min(pcoas_divgrass$vectors[, 2]), max(pcoas_divgrass$vectors[, 2])) +
  xlim(
    min(pcoas_divgrass$vectors[, 1]) - 1.2,
    max(pcoas_divgrass$vectors[, 1]) + 1.7
  ) +
  theme_bw() +
  geom_point(
    data = pc_axes_sp_notinsu,
    aes(color = pernotinsu),
    size = ptsize - 0.2
  ) +
  scale_color_gradient(low = "lightblue", high = "darkblue") +
  geom_point(
    data = pc_axes_distinct,
    shape = 1,
    col = "#FA6C67",
    size = ptsize
  ) +
  geom_point(
    data = focal_coords,
    shape = 1,
    color = "#FA6C67",
    size = ptsize + 2,
    stroke = 1.15
  ) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank()
  )

#BUPLEURUM BALDENSE up left
#LATHYRUS LATIFOLIUS up right
#DIPHASIASTRUM ALPINUM low left
#PETASITES HYBRIDUS low right

ggsave(
  file = here::here("outputs", "Fig_4d.tiff"),
  Fig_4d,
  width = 14,
  height = 12,
  dpi = 300,
  units = "cm",
  device = 'tiff'
)

#find examples
rownames(pc_axes_sp_notinsu) <- pc_axes_sp_notinsu$Species
plot(
  pc_axes_sp_notinsu$Axis1,
  pc_axes_sp_notinsu$Axis2,
  pch = 19,
  col = "blue",
  main = "Click on a point to get its ID"
)
clicked_points <- identify(
  pc_axes_sp_notinsu$Axis1,
  pc_axes_sp_notinsu$Axis2,
  labels = rownames(pc_axes_sp_notinsu),
  plot = TRUE
)
print(paste(
  "You clicked on ID(s):",
  rownames(pc_axes_sp_notinsu)[clicked_points]
))
pc_axes_sp_notinsu[
  pc_axes_sp_notinsu$Species %in% rownames(pc_axes_sp_notinsu)[clicked_points],
]
#----

####ONE CELL WITHIN A BUFFER----

focal_cell = "CAL1_658505"

per_buff_all <- 1

buffer_cells <- cbind.data.frame(
  reference = focal_cell,
  neighbor = alpine_buffer$neighbor_cell[
    alpine_buffer$focal_cell %in% focal_cell
  ]
)

sp_focal <- unique(all_com_divgrass[[focal_cell]])
sp_focal <- sp_focal[!is.na(sp_focal)]

set.seed(1)
if (per_buff_all < 1)
  buffer_cells <- buffer_cells[
    which(
      buffer_cells$neighbor %in%
        c(
          focal_cell,
          sample(
            buffer_cells$neighbor[-1],
            round(length(buffer_cells$neighbor) * per_buff_all)
          )
        )
    ),
  ]
buffer_cell_used <- buffer_cells[
  unlist(lapply(as.character(buffer_cells$neighbor), function(id) {
    length(unique(all_com_divgrass[[id]])) >= 5
  })),
]

occ_list <- all_com_divgrass
dist_mat <- disttrait_divgrass
D_insu <- median(disttrait_divgrass) * 0.3 #10% of the median value
D_thr <- 90

I_buffer <- do.call(
  rbind,
  (pbmcapply::pbmclapply(
    1:dim(buffer_cell_used)[1],
    function(j) {
      if (focal_cell != buffer_cell_used$neighbor[j]) {
        a <- Insu(
          insured = focal_cell,
          insurer = buffer_cell_used$neighbor[j],
          occ_list = all_com_divgrass,
          dist_mat = dist_mat,
          D_insu,
          D_thr
        )
        b <- Insu(
          insured = buffer_cell_used$neighbor[j],
          insurer = focal_cell,
          occ_list = all_com_divgrass,
          dist_mat = dist_mat,
          D_insu,
          D_thr
        )
        a$Outward <- b$Inward
        a
      }
    },
    mc.cores = parallel::detectCores()
  ))
)

#who are insured

insuwho <- pbmcapply::pbmclapply(
  2:dim(buffer_cell_used)[1],
  function(j) {
    #j=2

    Insu_who(
      insured = focal_cell,
      insurer = buffer_cell_used$neighbor[j],
      occ_list = all_com_divgrass,
      dist_mat = dist_mat,
      D_insu,
      D_thr
    )
  },
  mc.cores = parallel::detectCores()
)

insu_who_all <- insuwho[[1]]

for (i in 2:length(insuwho)) {
  insu_who_all <- merge(insu_who_all, insuwho[[i]])
}

rownames(insu_who_all) <- insu_who_all$Species
insu_who_all <- insu_who_all[, -1]

insu_who_all <- rowSums(insu_who_all)

#means

mean_one_buffer <- cbind.data.frame(
  cell_focal = I_buffer$insured[1],
  Inward = mean(I_buffer$Inward[-1], na.rm = TRUE),
  Outward = mean(I_buffer$Outward[-1], na.rm = TRUE),
  Nsp_loc = I_buffer$Div_insured[1],
  Div_distinct = I_buffer$Div_distinct[1],
  Nsp_reg = mean(I_buffer$Div_insurer[-1], na.rm = TRUE)
)


mean_one_buffer

#Map metrics

id <- which(Div_grass_french_alps$com %in% focal_cell)
reference_point <- Div_grass_french_alps[id, ]

# Select its neighbors and add the Insurance values
neighbors <- Div_grass_french_alps[alpine_within_radius[[id]], ]

neighbors <- neighbors[which(neighbors$com %in% I_buffer$insurer), ]

neighbors$inward = NA
for (i in 1:nrow(neighbors)) {
  neighbors$inward[i] = I_buffer$Inward[I_buffer$insurer %in% neighbors$com[i]]
}

neighbors$outward = NA
for (i in 1:nrow(neighbors)) {
  neighbors$outward[i] = I_buffer$Outward[
    I_buffer$insurer %in% neighbors$com[i]
  ]
}

# Create the buffer
buffer <- sf::st_buffer(reference_point, dist = 3000)
buffer_wgs84 <- sf::st_transform(buffer, crs = 4326)

#map
palette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "YlGn"))(255)

tmap::tmap_mode("view") # Enable interactive mode

tmap::tm_shape(Div_grass_french_alps) +
  # tmap::tm_polygons(col = "lightblue", border.col = "black") +
  tmap::tm_shape(Div_grass_french_alps) +
  tmap::tm_shape(buffer_wgs84) + # Add the buffer as a polygon
  tmap::tm_polygons(col = "white", fill_alpha = 0.3, border.col = "darkgreen") + # Customize buffer appearance
  tmap::tm_shape(Div_grass_french_alps) +
  tmap::tm_dots(fill = "gray", size = 0.5) +
  tmap::tm_shape(neighbors) +
  tmap::tm_symbols(
    col = "inward",
    palette = palette,
    size = 0.4,
    fill.scale = "cont"
  ) +
  tmap::tm_shape(reference_point) +
  tmap::tm_dots(fill = "red", size = 0.5) # Plot the selected reference point

#end----

####TWO CELLS Fig 4c----
library(gridExtra)
library(ggforce)

D_insu <- median(disttrait_divgrass) * 0.3
D_thr <- 90
dist_mat <- disttrait_divgrass

minus = 2
ptsize = 2.5
addymin = 1 - minus
remymax = 1.5 - minus
addxmin = 2 - minus
remxmax = 1.2 - minus

sink <- "CAL1_848810" #ex of fig3
source <- "SOP_T_1491.19.4" #ex of fig3

spe_insured <- all_com_divgrass[[sink]]
spe_insurer <- all_com_divgrass[[source]]

pc_insured <- data.frame(pcoas_divgrass$vectors[
  rownames(pcoas_divgrass$vectors) %in% spe_insured,
])
pc_insurer <- data.frame(pcoas_divgrass$vectors[
  rownames(pcoas_divgrass$vectors) %in% spe_insurer,
])

colnames(pc_insured) <- gsub("\\.", "", colnames(pc_insured))
colnames(pc_insurer) <- gsub("\\.", "", colnames(pc_insurer))

#Compute the localdisctinctivness

dist_insured <- dist_mat[
  colnames(dist_mat) %in% spe_insured,
  rownames(dist_mat) %in% spe_insured
]
insured_Di = funrar::distinctiveness_global(dist_insured)
colnames(insured_Di) <- c("Species", "insured_Di")

dist_insurer <- dist_mat[
  colnames(dist_mat) %in% spe_insurer,
  rownames(dist_mat) %in% spe_insurer
]
insurer_Di = funrar::distinctiveness_global(dist_insurer)
colnames(insurer_Di) <- c("Species", "insurer_Di")

insured_Di_Q <- as.numeric(as.numeric(quantile(
  insured_Di$insured_Di,
  probs = seq(0, 1, 0.01)
)[D_thr + 1]))
insured_spe_Di_Q <- insured_Di$Species[insured_Di$insured_Di > insured_Di_Q]
pc_insured_Di_Q <- pc_insured[rownames(pc_insured) %in% insured_spe_Di_Q, ]

insurer_Di_Q <- as.numeric(as.numeric(quantile(
  insurer_Di$insurer_Di,
  probs = seq(0, 1, 0.01)
)[D_thr + 1]))
insurer_spe_Di_Q <- insurer_Di$Species[insurer_Di$insurer_Di > insurer_Di_Q]
pc_insurer_Di_Q <- pc_insurer[rownames(pc_insurer) %in% insurer_spe_Di_Q, ]

#Project the species on the global PCA
col_gray <- '#8F8F8F10'

pc_axes <- data.frame(pcoas_divgrass$vectors)

colnames(pc_axes) <- gsub("\\.", "", colnames(pc_axes))

fig_sink <- ggplot(pc_axes, aes(x = Axis1, y = Axis2)) +
  geom_point(col = col_gray) +
  geom_vline(xintercept = 0, linetype = "dashed", col = "#9E9E9E") +
  geom_hline(yintercept = 0, linetype = "dashed", col = "#9E9E9E") +
  ylim(
    min(pcoas_divgrass$vectors[, 2]) + addymin,
    max(pcoas_divgrass$vectors[, 2]) - remymax
  ) +
  xlim(
    min(pcoas_divgrass$vectors[, 1]) + addxmin,
    max(pcoas_divgrass$vectors[, 1]) - remxmax
  ) +
  theme_bw() +
  geom_point(data = pc_insured, col = "black", size = ptsize) +
  geom_point(data = pc_insured_Di_Q, col = "#FA6C67", size = ptsize) +
  geom_circle(
    aes(x0 = Axis1, y0 = Axis2, r = D_insu),
    data = pc_insurer_Di_Q,
    col = "#52A153",
    linetype = 2
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank()
  )

fig_source <- ggplot(pc_axes, aes(x = Axis1, y = Axis2)) +
  geom_point(col = col_gray) +
  geom_vline(xintercept = 0, linetype = "dashed", col = "#9E9E9E") +
  geom_hline(yintercept = 0, linetype = "dashed", col = "#9E9E9E") +
  ylim(
    min(pcoas_divgrass$vectors[, 2]) + addymin,
    max(pcoas_divgrass$vectors[, 2]) - remymax
  ) +
  xlim(
    min(pcoas_divgrass$vectors[, 1]) + addxmin,
    max(pcoas_divgrass$vectors[, 1]) - remxmax
  ) +
  theme_bw() +
  geom_point(data = pc_insurer, col = "black", size = ptsize) +
  # geom_point(data = pc_insured_Di_Q,
  #            col = "#52A153",size=ptsize)+
  geom_point(data = pc_insurer_Di_Q, col = "#52A153", size = ptsize) +
  geom_circle(
    aes(x0 = Axis1, y0 = Axis2, r = D_insu),
    data = pc_insured_Di_Q,
    col = "#FA6C67",
    linetype = 2
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank()
  )

Insured_Insurer <- gridExtra::arrangeGrob(fig_sink, fig_source, ncol = 1)
ggsave(
  file = here::here("outputs", "Fig_4c.tiff"),
  Insured_Insurer,
  width = 7,
  height = 14,
  dpi = 300,
  units = "cm",
  device = 'tiff'
)

#----

####SPATIAL AUTOCORELATION----

# Extract coordinates
coords <- sf::st_coordinates(divgrass_insurance)

# Define spatial neighbors using k=8 nearest neighbors
nb <- spdep::knn2nb(spdep::knearneigh(coords, k = 8))
lw <- spdep::nb2listw(nb, style = "W")

# Compute Join Count Statistics
join_count_results <- spdep::joincount.multi(divgrass_insurance$Category, lw)

#Convert to data frame for easier manipulation
jc_df <- as.data.frame(join_count_results)
jc_df$`z-value` <- as.numeric(jc_df$`z-value`)

# Compute two-tailed p-values from z-values
jc_df$p_value <- 2 * (1 - pnorm(abs(jc_df$`z-value`)))

# Add a significance label
jc_df$significance <- cut(
  jc_df$p_value,
  breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
  labels = c("***", "**", "*", "ns")
)

# Print results with significance
print(jc_df)

#end----
