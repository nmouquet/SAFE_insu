###################################################################################################
#' Spatial Insurance, analysis for birds worldwide
#'
#'Uses: 
#'
#'    -  occurrences from IUCN used in Loiseau et al. 2020 data/Loiseau2021/birds_occurrences.RData
#'    -  birds traits from the AVOINET database (Tobias et al. 2022) data/AVONET/ELEData/TraitData/AVONET1_BirdLife.csv
#'    -  WDPA shapefiles data/BIG_FILES/WDPA_Mar2025_Public_shp/WDPA_Mar2025_Public_shp_1, WDPA_Mar2025_Public_shp_2, 
#'       WDPA_Mar2025_Public_shp_3
#'    -  Human impact from Mu et al. 2022 data/BIG_FILES/hfp2018.tiff
#'  
#'    Loiseau N., Mouquet N. et al., Global distribution and conservation status of ecologically rare mammal and bird species. 
#'    Nature Communications 11,  (2020).
#'    
#'    H. Mu, X. Li, Y. Wen, J. Huang, P. Du, W. Su, S. Miao, M. Geng, A global record of annual terrestrial Human Footprint 
#'    dataset from 2000 to 2018. Scientific Data 9, 176 (2022).
#'    
#'    Tobias J.A. et al., AVONET: morphological, ecological and geographical data for all birds. Ecology Letters 25, 581-597 (2022).
#'  
#'Produces:
#'
#'    - results/birds_insurance.RData : insurance for all cells 
#'    - results/birds_insurance_pa.RData : insurance for all cells + % of protected area
#'    - results/wdpa_I_II_union.RData
#'    - outputs/Fig_5a_1.tiff
#'    - outputs/Fig_5a_2.tiff
#'    - outputs/Fig_5a_3.tiff
#'    - outputs/Fig_5b.tiff
#'    - outputs/Fig_5c.tiff
#'    - outputs/Fig_5d.tiff
#'    - outputs/Fig_6.tiff
#'    - outputs/S3_Fig_1
#'    - outputs/S3_Fig_2
#'    - outputs/S3_Fig_3
#'    - outputs/S3_Fig_4
#'    - outputs/S3_Fig_5
#'        
#' @author Nicolas Mouquet, \email{nicolas.mouquet@@cnrs.fr},
#' @date 2023/05/22 first created, major update 2024/05/10 & 2025/03/01
##################################################################################################
rm(list=ls(all=TRUE))
library(ggplot2)
library(dplyr)

source(here::here("R","world_map.R"))
source(here::here("R","buffer.R"))
source(here::here("R","get_cell_coords.R"))
source(here::here("R","map_variable.R"))

color_s <- rev(RColorBrewer::brewer.pal(n = 8, name = "RdBu"))

#### Functions----
Insu <- function(insured,insurer,occ_list,dist_mat,D_insu,D_thr){
  
  # insured <- focal_cell
  # insurer <-buffer_cells_used$neighbor[2]
  # occ_list <- occ_list
  # dist_mat <- disttrait
  # D_insu <- 0.05
  # D_thr <- 90

  spe_insured <- occ_list[[as.character(insured)]]
  spe_insurer <- occ_list[[as.character(insurer)]]
  
  #Distance and insurance 
  dist_insured <-dist_mat[colnames(dist_mat)%in%spe_insured,rownames(dist_mat)%in%spe_insured]
  insured_Di = funrar::distinctiveness_global(dist_insured)
  colnames(insured_Di) <- c("Species","Di")
  
  QD <- as.numeric(quantile(insured_Di$Di, probs = seq(0, 1, 0.01))[D_thr+1])
  
  if (D_thr==0) {
    insured_D=insured_Di$Species
  } else {
    insured_D=insured_Di$Species[insured_Di$Di>=QD]
  }
  
  Insu <- do.call(rbind,lapply(insured_D, function(sp_insured){
    #sp_insured <- insured_spe_Di80[1]
    
    #is the insured sp in the insurerboring com ? 
    #threshold for insurance radius 
    
    if (sp_insured%in%spe_insurer){
      
      dist_insurer <-dist_mat[colnames(dist_mat)%in%spe_insurer,rownames(dist_mat)%in%spe_insurer]
      Insurance <- sum(dist_insurer[sp_insured,]<D_insu)
      
    } else {
      dist_insurer <-dist_mat[colnames(dist_mat)%in%c(spe_insurer,sp_insured),rownames(dist_mat)%in%c(spe_insurer,sp_insured)]
      Insurance <- sum(dist_insurer[sp_insured,]<D_insu)-1
    }
    
    cbind.data.frame(Species=sp_insured,Insurance=Insurance)
    
  }))
  
  Inward <- round(100*sum(Insu$Insurance>0)/nrow(Insu),2)
  
  rm(dist_insured)
  rm(Insu)
  
  cbind.data.frame(D_insu=D_insu,D_thr=D_thr,insured=insured,insurer=insurer,Div_insured=length(spe_insured),
                   Div_insurer=length(spe_insurer),Inward=Inward,Div_distinct=length(insured_D))
  
}
#----

#### Set param, Load Data & Compute functional space, Fig_5b ----

    options(taxa = "birds", dispersion_max = 350000)

  ## Load data from Loiseau_et_al
    occ_list <- get(load(file = here::here("data","Loiseau2021","birds_occurrences.RData")))
    Loiseau_et_al    <- get(load(file = here::here("data","Loiseau2021","Loiseau_et_al.RData")))
    Loiseau_et_al <- Loiseau_et_al[Loiseau_et_al$class%in%"birds",]
    Loiseau_et_al <- Loiseau_et_al[,c("species_code","scientific_name")]
    Loiseau_et_al <- Loiseau_et_al[Loiseau_et_al$species_code%in%unique(unlist(occ_list)),]
    
 ## Load bird traits from AVOINET and compute the disttrait
    
    traits <- read.csv(file=here::here("data","AVONET","ELEData","TraitData","AVONET1_BirdLife.csv"))
    traits <- traits[c("Species1", "Beak.Length_Culmen", "Beak.Length_Nares", 
                       "Beak.Width", "Beak.Depth", "Tarsus.Length", "Wing.Length", 
                       "Kipps.Distance", "Secondary1", "Hand.Wing.Index", "Tail.Length", "Mass")]
    colnames(traits)[1] <- "scientific_name" 
    rownames(traits) <- traits$scientific_name
    
    traits <- merge(Loiseau_et_al,traits)
    rownames(traits) <- traits$species_code

    traits_for_hypervolumes <- na.omit(traits[-c(1:2)])
    rownames(traits_for_hypervolumes) <- na.omit(traits)$species_code
    traits_for_hypervolumes[, 1:11] <- scale(log10(traits_for_hypervolumes[, 1:11]), center = F)
    # flightless birds have wing length == 0
    traits_for_hypervolumes[traits_for_hypervolumes$Wing.Length < 0.1, c("Wing.Length", "Secondary1", "Tail.Length")] <- 0
  
    disttrait = funrar::compute_dist_matrix(traits_for_hypervolumes,metric = "euclidean")
    
    #pcoas<-ape::pcoa(disttrait) #long to run
    #save(pcoas,file=here::here("results","birds_pcoa_avoinet.RData"))
    pcoas <- get(load(file = here::here("results","birds_pcoa_avoinet.RData")))
    
    global_Di = funrar::distinctiveness_global(disttrait)
    colnames(global_Di) <- c("Species","global_di")
    
    hist(global_Di$global_di)
    
  ##plot global map 
    
    tiff(filename = "outputs/Fig_5b.tiff", width = 2000, height = 1300, res = 300)
    world_map_robin(id=NA,color_s = color_s)
    if (dev.cur() != 1) dev.off()
    
#----

#### VISU FUNCTIONAL SPACE---- 

  library(ggplot2)
    
  load(file = here::here("results","birds_pcoa_avoinet.RData"))
    

  ## Eigenvalues of the RFS
  
  df <- data.frame("Eigenvalues" = paste0('E', 1:8),
                   "Value"       = pcoas$"values"$"Relative_eig"[1:8])
  
  a <- ggplot(data = df, aes(x = Eigenvalues, y = Value)) + 
    geom_bar(stat = "identity", fill = "steelblue") +
    geom_text(aes(label = round(Value, digits = 2)), vjust = 1.6, color = "white",
              size = 3.5) +
    theme_bw()
  
  ## Visualization of the RFS
  
    pc_axes <- data.frame(pcoas$vectors)
    colnames(pc_axes) <- gsub("\\.","",colnames(pc_axes))
  
  b <- ggplot(pc_axes, aes(Axis1, Axis2)) + 
    geom_point(colour = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 0.5) + 
    xlab(paste0("PC1 (",round(df$Value[1]*100,1)," %)"))+
    ylab(paste0("PC2 (",round(df$Value[2]*100,1)," %)"))+
    theme_bw()
  
  #S2_Fig_X <- gridExtra::arrangeGrob(a,b, ncol=2)
  #ggsave(file=here::here("outputs","S2_Fig_X.tiff"), S2_Fig_X1,width = 20, height = 10, dpi = 300, units = "cm", device='tiff') 
  
  #check some species (if needed)
    #plot(pc_axes$Axis1, pc_axes$Axis2, pch = 19, col = "blue", main = "Click on a point to get its ID")
    #clicked_points <- identify(pc_axes$Axis1, pc_axes$Axis2, labels = rownames(pc_axes), plot = TRUE)
    #print(paste("You clicked on ID(s):", traits$scientific_name[traits$species_code%in%rownames(pc_axes)[clicked_points]]))

    
#----
        
####TWO CELLS----
  library(gridExtra)
  library(ggforce)
        
  D_thr <- 90
  D_insu <- 0.1
  dist_mat <- disttrait
  
  #ex inward >> Outward
  # insured <- "180652"
  # insurer <- "179475"
  
  #ex inward >> Outward
  insured <- "61256"
  insurer <- "60673"
        
  spe_insured <- occ_list[[insured]]
  spe_insurer <- occ_list[[insurer]]
        
  pc_insured <- data.frame(pcoas$vectors[rownames(pcoas$vectors)%in%spe_insured,])
  pc_insurer <- data.frame(pcoas$vectors[rownames(pcoas$vectors)%in%spe_insurer,])
  
  colnames(pc_insured) <- gsub("\\.","",colnames(pc_insured))
  colnames(pc_insurer) <- gsub("\\.","",colnames(pc_insurer))
  
        
  #Compute the localdisctinctivness  
        
    dist_insured <-dist_mat[colnames(dist_mat)%in%spe_insured,rownames(dist_mat)%in%spe_insured]
    insured_Di = funrar::distinctiveness_global(dist_insured)
    colnames(insured_Di) <- c("Species","insured_Di")
        
    dist_insurer <-dist_mat[colnames(dist_mat)%in%spe_insurer,rownames(dist_mat)%in%spe_insurer]
    insurer_Di = funrar::distinctiveness_global(dist_insurer)
    colnames(insurer_Di) <- c("Species","insurer_Di")
        
    insured_Di_Q <- as.numeric(as.numeric(quantile(insured_Di$insured_Di, probs = seq(0, 1, 0.01))[D_thr+1]))
    insured_spe_Di_Q <- insured_Di$Species[insured_Di$insured_Di>insured_Di_Q]
    pc_insured_Di_Q <- pc_insured[rownames(pc_insured)%in%insured_spe_Di_Q,]
        
    insurer_Di_Q <- as.numeric(as.numeric(quantile(insurer_Di$insurer_Di, probs = seq(0, 1, 0.01))[D_thr+1]))
    insurer_spe_Di_Q <- insurer_Di$Species[insurer_Di$insurer_Di>insurer_Di_Q]
    pc_insurer_Di_Q <- pc_insurer[rownames(pc_insurer)%in%insurer_spe_Di_Q,]
        
        #Project the species on the global PCA 
        col_gray <- '#a6a6a510'
        
        pc_axes <- data.frame(pcoas$vectors)
        
        colnames(pc_axes) <- gsub("\\.","",colnames(pc_axes))
        
        fig_pc_insured_1 <- ggplot(pc_axes, aes(x = Axis1, y = Axis2))+
          geom_point(col = col_gray)+
          ylim(min(pcoas$vectors[,2]),max(pcoas$vectors[,2]))+
          xlim(min(pcoas$vectors[,1]),max(pcoas$vectors[,1]))+
          theme_bw()+
          geom_point(data = pc_insured, 
                     col = "black")+
          geom_point(data = pc_insured_Di_Q, 
                     col = "red")
        
        fig_pc_insurer_1 <- ggplot(pc_axes, aes(x = Axis1, y = Axis2))+
          geom_point(col = col_gray)+
          ylim(min(pcoas$vectors[,2]),max(pcoas$vectors[,2]))+
          xlim(min(pcoas$vectors[,1]),max(pcoas$vectors[,1]))+
          theme_bw()+
          geom_point(data = pc_insurer, 
                     col = "black")+
          geom_point(data = pc_insured_Di_Q,shape=1,
                     col = "red")+
          geom_circle(aes(x0 = Axis1, y0 = Axis2, r = D_insu), data = pc_insured_Di_Q,
                      col = "#F78C8C",linetype=2)
        
        fig_pc_insured_2 <- ggplot(pc_axes, aes(x = Axis1, y = Axis2))+
          geom_point(col = col_gray)+
          ylim(min(pcoas$vectors[,2]),max(pcoas$vectors[,2]))+
          xlim(min(pcoas$vectors[,1]),max(pcoas$vectors[,1]))+
          theme_bw()+
          geom_point(data = pc_insurer, 
                     col = "black")+
          geom_point(data = pc_insurer_Di_Q, 
                     col = "red")
        
        fig_pc_insurer_2 <- ggplot(pc_axes, aes(x = Axis1, y = Axis2))+
          geom_point(col = col_gray)+
          ylim(min(pcoas$vectors[,2]),max(pcoas$vectors[,2]))+
          xlim(min(pcoas$vectors[,1]),max(pcoas$vectors[,1]))+
          theme_bw()+
          geom_point(data = pc_insured, 
                     col = "black")+
          geom_point(data = pc_insurer_Di_Q,shape=1,
                     col = "red")+
          geom_circle(aes(x0 = Axis1, y0 = Axis2, r = D_insu), data = pc_insurer_Di_Q,
                      col = "#F78C8C",linetype=2)
        
        Insured_Insurer <- gridExtra::arrangeGrob(fig_pc_insured_1,fig_pc_insurer_1,fig_pc_insurer_2,fig_pc_insured_2, ncol=2)
        ggsave(file=here::here("outputs",paste0(insured,"_vs_",insurer,".tiff")), Insured_Insurer,width = 30, height = 30, dpi = 300, units = "cm", device='tiff') 
#----
        
####COMPUTE ALL BUFFERS only once----
        
  all_cell_occ <- as.integer(unique(names(occ_list)))
  max_disp=350000 #distance based on M. W. McKnight, P. S. White, R. I. McDonald, J. F. Lamoreux, W. Sechrest, R. S. Ridgely, S. N. Stuart, Putting beta-diversity on the map: Broad-scale congruence and coincidence in the extremes. Plos Biology 5, 2424-2432 (2007).
  
  buffer(cell_id=all_cell_occ,plot=FALSE,buff_rad=max_disp,distance=TRUE)

  get_buffer <- function(id_cell){
          an.error.occured <- FALSE
          tryCatch( {buffer_cells <<- buffer(cell_id=id_cell,plot=FALSE,buff_rad=max_disp,distance=TRUE) }
                    , error = function(e) {an.error.occured <<- TRUE})
          
          save(buffer_cells,file=here::here("results",'BIG_FILES',"all_buffer",paste0("cell_",id_cell,".RData")))
        }
        
  pbmcapply::pbmcmapply(get_buffer, all_cell_occ, mc.cores = paralell::parallel::detectCores()-1)
        
#----

####LOOP OVER MANY CELLS WITH BUFFERS, Fig_5a_1,Fig_5c,Fig_5d----
  #veryyy long 

  sp_cell <- names(occ_list)
  
  set.seed(27)
  disp_max=350000
  per_buff_all <- 1 #used to sample the neighboring cells if needed (to speed up computation)
  occ_list <- occ_list
  dist_mat <- disttrait
  D_insu <- median(disttrait)*0.1
  D_thr <- 90
  
  #sp_cell <- sample(sp_cell,10000) #used to sample the cells if needed (to speed up computation)

  error_log_file <- here::here("error_log.txt")
  file.create(error_log_file)
  birds_insurance <- do.call(rbind,pbmcapply::pbmclapply(sp_cell, function(focal_cell){
    
    #focal_cell=sp_cell[1]
    
    tryCatch({
      
      an.error.occured <- FALSE
      tryCatch( {
        buffer_cells <- get(load(file = here::here("data",'BIG_FILES', "all_buffer",paste0("cell_",focal_cell,".RData"))))
        buffer_cells <- buffer_cells[buffer_cells$dist<=disp_max,]
      }, error = function(e) {an.error.occured <<- TRUE})
      
      buffer_cells <- buffer_cells[!buffer_cells$reference==buffer_cells$neighbor,]
      set.seed(1)
      if (per_buff_all<1) buffer_cells_used <- buffer_cells[which(buffer_cells$neighbor %in% c(focal_cell,sample(buffer_cells$neighbor[-1],round(length(buffer_cells$neighbor)*per_buff_all)))),]
      
      if (per_buff_all==1) buffer_cells_used <- buffer_cells
      
      buffer_cells_used <- buffer_cells_used[unlist(lapply(as.character(buffer_cells_used$neighbor),function(id){length(unique(occ_list[[id]]))>=5})),]
      
      sp_focal <- unique(occ_list[[as.character(focal_cell)]])
      
      #compute insurance when there are more than 10 species only in the focal cell
      if  ((length(sp_focal)>=10) & (an.error.occured==FALSE) & (nrow(buffer_cells_used) > 2)){
        
        sp_focal <- unique(occ_list[[as.character(focal_cell)]])
        
        I_buffer <- do.call(rbind,( parallel::mclapply(1:dim(buffer_cells_used)[1],function(j){
          
          #j=1
          
          a <- Insu(insured=focal_cell,insurer=buffer_cells_used$neighbor[j],occ_list,dist_mat=disttrait,D_insu,D_thr)
          b <- Insu(insured=buffer_cells_used$neighbor[j],insurer=focal_cell,occ_list,dist_mat=disttrait,D_insu,D_thr)
          a$Outward <- b$Inward
          a$dist_cells= dist_cells=buffer_cells_used$dist[j]
          a
          
        },mc.cores = 1)))
        
        #Compute betadiv (subset neighbor cell with more than 5 species)
        
        birds_insurance_buffer <- as.character(buffer_cells_used$neighbor)
        
        all_species <- unique(unlist(lapply(birds_insurance_buffer, function(id){occ_list[[id]]})))
        all_species <- all_species[order(all_species)]
        
        buffer_com <- as.data.frame(matrix(0, length(birds_insurance_buffer),length(all_species)))
        rownames(buffer_com) <- birds_insurance_buffer
        colnames(buffer_com) <- all_species
        for (i in 1:length(birds_insurance_buffer)){
          buffer_com[i,colnames(buffer_com) %in% occ_list[[rownames(buffer_com)[i]]]]=1
        }

        #means 
        
        mean_one_buffer <- cbind.data.frame(D_insu=I_buffer$D_insu[1],
                                            D_thr=I_buffer$D_thr[1],
                                            cell_focal=I_buffer$insured[1],
                                            Inward=mean(I_buffer$Inward[-1],na.rm= TRUE),
                                            Outward=mean(I_buffer$Outward[-1],na.rm= TRUE),
                                            Nsp_loc=I_buffer$Div_insured[1],
                                            Div_distinct=I_buffer$Div_distinct[1],
                                            Nsp_reg=mean(I_buffer$Div_insurer[-1],na.rm = TRUE),
                                            gamma=length(all_species),
                                            cells=nrow(buffer_cells),
                                            cells_used=nrow(buffer_cells_used))
        # Free up memory
        rm(I_buffer, all_species)
        
        #return mean_one_buffer
        mean_one_buffer
        
      } else {NULL}
      
      
    }, error = function(e) {
      
      # Append error information to the log file using cat()
      cat(paste(Sys.time(), "- Error at cell =", focal_cell," Nsploc =",length(sp_focal),"\n"), 
          file = error_log_file, append = TRUE)
    })
    
    
    
  },mc.cores = parallel::detectCores()-1))
  message(readLines(error_log_file))
  birds_insurance <- birds_insurance[!is.na(birds_insurance$Inward), ]
  
  #Add coordinates
  
    # we retrieve cell ids from another raster
    grd <- suppressWarnings(
      raster::raster(here::here("data", "Loiseau2021","grid_area_with_cell_ids.tif")))
  
    # we use a raster with species richness as base line 
    ras <- suppressWarnings(raster::raster(here::here("data","Loiseau2021","birds_richness.tif")))
    
    coords <- raster::xyFromCell(ras, 1:raster::ncell(ras))
    sp <- sp::SpatialPoints(coords, proj4string = raster::crs(ras))
    sp_transformed <- sp::spTransform(sp, sp::CRS("+proj=longlat +datum=WGS84"))
    
    longlat <- cbind.data.frame("real_cell_id" = 1:raster::ncell(grd),
                                longitude=sp::coordinates(sp_transformed)[, 1],
                                latitude=sp::coordinates(sp_transformed)[, 2])
    
    ids <- data.frame("real_cell_id" = 1:raster::ncell(grd), 
                      "old_cell_id"  = grd[])
    
    data <- merge(ids, longlat)
  
    birds_insurance <- merge(birds_insurance,data, by.x = "cell_focal", by.y = "old_cell_id", 
                           all = FALSE)
  
  # Classify points into Source, Sink, and Intermediate
  
    low_threshold <- quantile(birds_insurance$Inward - birds_insurance$Outward, 0.1, na.rm = TRUE)
    high_threshold <- quantile(birds_insurance$Inward - birds_insurance$Outward, 0.9, na.rm = TRUE)
    birds_insurance$SS_status <- ifelse(birds_insurance$Inward - birds_insurance$Outward >= high_threshold, "Sink",
                                        ifelse(birds_insurance$Inward - birds_insurance$Outward <= low_threshold, "Source", "Intermediate"))
    birds_insurance$SS_status <- factor(birds_insurance$SS_status, levels = c("Source", "Intermediate", "Sink"))
    
  #Save  
    save(birds_insurance,file=here::here("results","birds_insurance.RData"))
    
  #load
    load(here::here('results','birds_insurance.RData'))
    
  # look at the categories to check S3_Fig_1
    S3_Fig_1 <- ggplot(birds_insurance, aes(x = Inward, y = Outward, color = SS_status)) +
      geom_point(alpha = 0.7) +  # Scatter plot of points
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  # 1:1 line
      scale_color_manual(values = c("Source" = "#7DBD80", "Sink" = "#FA6C67", "Intermediate" = "#FFF7BA")) +  # Define colors
      theme_bw() +
      labs(x = "Inward", y = "Outward", color = "Insurance status") +
      theme(legend.position = "right")
    
    ggsave(file=here::here("outputs","S3_Fig_1.tiff"), S3_Fig_1,width = 13, height = 9, dpi = 300, units = "cm", device='tiff') 
    
    
  # look at species richness within each categories S3_Fig_2
    birds_insurance$SS_status <- factor(birds_insurance$SS_status, levels = c("Source", "Intermediate", "Sink"))
    
    S3_Fig_2 <- ggplot(birds_insurance, aes(x = SS_status, y = Nsp_loc,fill = SS_status))+
      geom_boxplot(width = .5, outlier.shape = NA) + 
      ylab("Species richness") +
      #xlab("Status")+
      scale_fill_manual(values = c("Source" = "#7DBD80", "Intermediate" = "#FFF7BA", "Sink" = "#FA6C67")) +
      theme_bw()+
      theme(axis.text.x = element_text(size = 12),
            axis.title.x = element_blank(),
            axis.text.y = element_text(size = 8),
            axis.title.y = element_text(size = 12),
            legend.position = "none")+
      ylim(0,500)+
      scale_x_discrete(labels = c("Source", "Intermediate", "Sink"))+
      ggplot2::theme(axis.line.x  = ggplot2::element_line(linetype = "blank"),
                     axis.ticks.x = ggplot2::element_blank())
    
    ggsave(file=here::here("outputs","S3_Fig_2.tiff"), S3_Fig_2,width = 14, height = 10, dpi = 300, units = "cm", device='tiff') 
    
    kruskal.test(Nsp_loc ~ SS_status, data = birds_insurance)
    FSA::dunnTest(birds_insurance$Nsp_loc ~ birds_insurance$SS_status, method="bonferroni")    
    
    #Look at the number of cells within buffers for each Category (if needed)
      birds_insurance$SS_status <- factor(birds_insurance$SS_status, levels = c("Source", "Intermediate", "Sink"))
      
      ggplot(birds_insurance, aes(x = SS_status, y = cells,fill = SS_status))+
        geom_boxplot(width = .5, outlier.shape = NA) + 
        ylab("Number of cells within buffers") +
        #xlab("Status")+
        scale_fill_manual(values = c("Source" = "#7DBD80", "Intermediate" = "#FFF7BA", "Sink" = "#FA6C67")) +
        theme_bw()+
        theme(axis.text.x = element_text(size = 12),
              axis.title.x = element_blank(),
              axis.text.y = element_text(size = 8),
              axis.title.y = element_text(size = 12),
              legend.position = "none")+
        ylim(0,200)+
        scale_x_discrete(labels = c("Source", "Intermediate", "Sink"))+
        ggplot2::theme(axis.line.x  = ggplot2::element_line(linetype = "blank"),
                       axis.ticks.x = ggplot2::element_blank())
    
  #plot with quantiles Fig_5c

    library(ggplot2)
    library(RColorBrewer)
    birds_insurance <- birds_insurance[order(birds_insurance$Nsp_loc,decreasing = F),]
    Fig_5c <- ggplot(birds_insurance, aes(x=Inward, y=Outward, color=Nsp_loc)) +
      geom_point(alpha=0.8,size=0.9) +
      ylab("Outward insurance")+xlab("Inward insurance")+
      scale_colour_gradientn(colours = color_s)+
      theme_bw()+
      theme(legend.position="none")+
      geom_abline (slope=1, intercept=as.numeric(quantile(birds_insurance$Outward-birds_insurance$Inward,prob=0.9,na.rm = T)),linetype = "dashed", size=0.5)+
      geom_abline (slope=1, intercept=as.numeric(quantile(birds_insurance$Outward-birds_insurance$Inward,prob=0.1,na.rm = T)),linetype = "dashed", size=0.5)+
      ylim(0,100)+xlim(0,100)
    ggsave(file=here::here("outputs","Fig_5c.tiff"), Fig_5c,width = 6, height = 6, dpi = 300, units = "cm", device='tiff') 
  
  #Select the functional sources and sinks and plot the global map Fig_5a_1 
  
    birds_insurance$SS=((birds_insurance$Outward-birds_insurance$Inward)+abs(min(birds_insurance$Outward-birds_insurance$Inward,na.rm = T)))/200
    birds_insurance$SS[(birds_insurance$Outward-birds_insurance$Inward>= quantile(birds_insurance$Outward-birds_insurance$Inward,prob=0.9,na.rm = TRUE))]=1
    birds_insurance$SS[(birds_insurance$Outward-birds_insurance$Inward<= quantile(birds_insurance$Outward-birds_insurance$Inward,prob=0.1,na.rm = TRUE))]=0
  
    if (dev.cur() != 1) dev.off()
      tiff(filename = "outputs/Fig_4a_1.tiff", width = 2000, height = 1300, res = 300)
      color_ss <- c("#FA6964","#FFF7BA","#57BD5C")
      world_map2_robin(birds_insurance[ , c("cell_focal", "SS")],color_s=color_ss) 
    if (dev.cur() != 1) dev.off()
  
  #Find examples
  
    #sinks 
      ex <- birds_insurance[((birds_insurance$Inward>60) & (birds_insurance$Outward<40)),]
      ex <- ex[!is.na(ex$D_insu),]
      ex[ex$Nsp_loc>150,]
    
    #sources
      ex <- birds_insurance[((birds_insurance$Inward<40) & (birds_insurance$Outward>60)),]
      ex <- ex[!is.na(ex$D_insu),]
      ex[ex$Nsp_loc>120,]
    
     
#end----
    
####Fig_5a_2 Source within a buffer----
    
    focal_cell=61549
    disp_max=350000
    std <- TRUE
    per_buff_all <- 1
    occ_list <- occ_list
    dist_mat <- disttrait
    D_insu <- median(disttrait)*0.1
    D_thr <- 90
    
    if (dev.cur() != 1) dev.off()
    world_map_robin(id=focal_cell,color_s = color_s)
    
    an.error.occured <- FALSE
    tryCatch( {
      buffer_cells <- get(load(file = here::here("data", "BIG_FILES","all_buffer",paste0("cell_",focal_cell,".RData"))))
      buffer_cells <- buffer_cells[buffer_cells$dist<=disp_max,]
    }, error = function(e) {an.error.occured <<- TRUE})
    
    sp_focal <- unique(occ_list[[as.character(focal_cell)]])
    
    set.seed(1)
    if (per_buff_all<1) buffer_cells <- buffer_cells[which(buffer_cells$neighbor %in% c(focal_cell,sample(buffer_cells$neighbor[-1],round(length(buffer_cells$neighbor)*per_buff_all)))),]
    buffer_cells_used <- buffer_cells[unlist(lapply(as.character(buffer_cells$neighbor),function(id){length(unique(occ_list[[id]]))>=5})),]
    
    I_buffer <- do.call(rbind,(pbmcapply::pbmclapply(1:dim(buffer_cells_used)[1],function(j){
      
      #j=1
      
      a <- Insu(insured=focal_cell,insurer=buffer_cells_used$neighbor[j],occ_list,dist_mat=disttrait,D_insu,D_thr)
      b <- Insu(insured=buffer_cells_used$neighbor[j],insurer=focal_cell,occ_list,dist_mat=disttrait,D_insu,D_thr)
      a$Outward <- b$Inward
      a$dist_cells= dist_cells=buffer_cells_used$dist[j]
      a
      
      
    },mc.cores = parallel::detectCores())))
    
    all_cells_buffer <- as.character(buffer_cells_used$neighbor)
    all_species <- unique(unlist(lapply(all_cells_buffer, function(id){occ_list[[id]]})))
    all_species <- all_species[order(all_species)]
    
    mean_one_buffer <- cbind.data.frame(cell_focal=I_buffer$insured[1],
                                        Inward=mean(I_buffer$Inward[-1],na.rm= TRUE),
                                        Outward=mean(I_buffer$Outward[-1],na.rm= TRUE),
                                        Nsp_loc=I_buffer$Div_insured[1],
                                        Div_distinct=I_buffer$Div_distinct[1],
                                        Nsp_reg=mean(I_buffer$Div_insurer[-1],na.rm = TRUE),
                                        gamma=length(all_species),
                                        cells=nrow(I_buffer))
    
    mean_one_buffer
    
    #Map metrics 
    
    I_buffer$insurer <- as.numeric(I_buffer$insurer)
    I_buffer$insured <- I_buffer$insured
    
    if (dev.cur() != 1) dev.off()
    tiff(filename = "outputs/Fig_5a_2.tiff", width = 2000, height = 1300, res = 300)
    par(mfrow = c(1, 3))
      map_variable(I_buffer[,-c(1:3)], "Div_insurer",disp_max=disp_max,reference=TRUE, grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(name = "RdBu", 9)))(255)) #YlOrRd
      map_variable(I_buffer[,-c(1:3)], "Inward",disp_max=disp_max,reference=TRUE, breaks = seq(0, 100, length.out = 256),col=grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "YlGn"))(255))
      map_variable(I_buffer[,-c(1:3)], "Outward",disp_max=disp_max,reference=TRUE, breaks = seq(0, 100, length.out = 256),col=grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "YlGn"))(255))
    if (dev.cur() != 1) dev.off()
    
    # if (dev.cur() != 1) dev.off()
    # tiff(filename = "outputs/Fig_4a_2.tiff", width = 2000, height = 1300, res = 300)
    # par(mfrow = c(1, 3))
    # map_variable(I_buffer[,-c(1:3)], "Div_insurer",disp_max=disp_max,reference=TRUE, rev(viridis::viridis(255))) 
    # map_variable(I_buffer[,-c(1:3)], "Inward",disp_max=disp_max,reference=TRUE, breaks = seq(0, 100, length.out = 256),col=rev(viridis::viridis(255)))
    # map_variable(I_buffer[,-c(1:3)], "Outward",disp_max=disp_max,reference=TRUE, breaks = seq(0, 100, length.out = 256),col=rev(viridis::viridis(255)))
    # if (dev.cur() != 1) dev.off()
    
#end----
    
####Fig_5a_3 Sink within a buffer----
  
    focal_cell=180652
    disp_max=350000
    std <- TRUE
    per_buff_all <- 1
    occ_list <- occ_list
    dist_mat <- disttrait
    D_insu <- median(disttrait)*0.1
    D_thr <- 90
    
    if (dev.cur() != 1) dev.off()
    world_map_robin(id=focal_cell,color_s = color_s)
    
    an.error.occured <- FALSE
    tryCatch( {
      buffer_cells <- get(load(file = here::here("data", "BIG_FILES","all_buffer",paste0("cell_",focal_cell,".RData"))))
      buffer_cells <- buffer_cells[buffer_cells$dist<=disp_max,]
    }, error = function(e) {an.error.occured <<- TRUE})
    
    sp_focal <- unique(occ_list[[as.character(focal_cell)]])
    
    set.seed(1)
    if (per_buff_all<1) buffer_cells <- buffer_cells[which(buffer_cells$neighbor %in% c(focal_cell,sample(buffer_cells$neighbor[-1],round(length(buffer_cells$neighbor)*per_buff_all)))),]
    buffer_cells_used <- buffer_cells[unlist(lapply(as.character(buffer_cells$neighbor),function(id){length(unique(occ_list[[id]]))>=5})),]
    
    I_buffer <- do.call(rbind,(pbmcapply::pbmclapply(1:dim(buffer_cells_used)[1],function(j){
      
      #j=1
      
      a <- Insu(insured=focal_cell,insurer=buffer_cells_used$neighbor[j],occ_list,dist_mat=disttrait,D_insu,D_thr)
      b <- Insu(insured=buffer_cells_used$neighbor[j],insurer=focal_cell,occ_list,dist_mat=disttrait,D_insu,D_thr)
      a$Outward <- b$Inward
      a$dist_cells= dist_cells=buffer_cells_used$dist[j]
      a
      
      
    },mc.cores = parallel::detectCores())))
    
    all_cells_buffer <- as.character(buffer_cells_used$neighbor)
    all_species <- unique(unlist(lapply(all_cells_buffer, function(id){occ_list[[id]]})))
    all_species <- all_species[order(all_species)]
    
    mean_one_buffer <- cbind.data.frame(cell_focal=I_buffer$insured[1],
                                        Inward=mean(I_buffer$Inward[-1],na.rm= TRUE),
                                        Outward=mean(I_buffer$Outward[-1],na.rm= TRUE),
                                        Nsp_loc=I_buffer$Div_insured[1],
                                        Div_distinct=I_buffer$Div_distinct[1],
                                        Nsp_reg=mean(I_buffer$Div_insurer[-1],na.rm = TRUE),
                                        gamma=length(all_species),
                                        cells=nrow(I_buffer))
    
    mean_one_buffer
    
    #Map metrics 
    
    I_buffer$insurer <- as.numeric(I_buffer$insurer)
    I_buffer$insured <- I_buffer$insured
    
    if (dev.cur() != 1) dev.off()
    tiff(filename = "outputs/Fig_5a_3.tiff", width = 2000, height = 1300, res = 300)
    par(mfrow = c(1, 3))
    map_variable(I_buffer[,-c(1:3)], "Div_insurer",reference=TRUE,disp_max=disp_max, grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(name = "RdBu", 9)))(255)) #YlOrRd
    map_variable(I_buffer[,-c(1:3)], "Inward",reference=TRUE,disp_max=disp_max, breaks = seq(0, 100, length.out = 256),col=grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "YlGn"))(255))
    map_variable(I_buffer[,-c(1:3)], "Outward",reference=TRUE,disp_max=disp_max, breaks = seq(0, 100, length.out = 256),col=grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "YlGn"))(255))
    if (dev.cur() != 1) dev.off()
    
    # if (dev.cur() != 1) dev.off()
    # tiff(filename = "outputs/Fig_4a_3.tiff", width = 2000, height = 1300, res = 300)
    # par(mfrow = c(1, 3))
    # map_variable(I_buffer[,-c(1:3)], "Div_insurer",reference=TRUE,disp_max=disp_max, rev(viridis::viridis(255))) 
    # map_variable(I_buffer[,-c(1:3)], "Inward",reference=TRUE,disp_max=disp_max, breaks = seq(0, 100, length.out = 256),col=rev(viridis::viridis(255)))
    # map_variable(I_buffer[,-c(1:3)], "Outward",reference=TRUE,disp_max=disp_max, breaks = seq(0, 100, length.out = 256),col=rev(viridis::viridis(255)))
    # if (dev.cur() != 1) dev.off()
    
#----
    
####SPATIAL AUTOCOR TEST---- 

  # Remove NA values in Inward
    data_clean <- birds_insurance[!is.na(birds_insurance$Inward), ]
  
  # Convert to sf object
    data_sf <- sf::st_as_sf(data_clean, coords = c("longitude", "latitude"), crs = 4326)  # WGS 84
  
  # Extract coordinates
    coords <- sf::st_coordinates(data_sf)
    
  # Define spatial neighbors using k=8 nearest neighbors
    nb <- spdep::knn2nb(spdep::knearneigh(coords, k = 8))
    listw <- spdep::nb2listw(nb, style = "W")
      
  # Compute Join Count Statistics
    join_count_results <- spdep::joincount.multi(data_clean$SS_status, listw)
  
  #Convert to data frame for easier manipulation
    jc_df <- as.data.frame(join_count_results)
    jc_df$`z-value` <- as.numeric(jc_df$`z-value`)
    
  # Compute two-tailed p-values from z-values
    jc_df$p_value <- 2 * (1 - pnorm(abs(jc_df$`z-value`)))
    
  # Add a significance label 
    jc_df$significance <- cut(jc_df$p_value,
                              breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                              labels = c("***", "**", "*", "ns"))
    
    # Joincount Expected  Variance    z-value p_value significance
    #   Source:Source              1588.6250   303.76  30.27303  233.52296       0          ***
    #   Intermediate:Intermediate 21980.5625 19442.64 118.96654  232.68369       0          ***
    #   Sink:Sink                  1675.1250   303.76  30.27303  249.24425       0          ***
    #   Intermediate:Source        2577.0625  4860.86 165.25414 -177.65663       0          ***
    #   Sink:Source                 338.3125   607.62  61.24944  -34.41101       0          ***
    #   Sink:Intermediate          2219.8125  4860.86 165.25414 -205.44711       0          ***
    #   Jtot                       5135.1875 10329.34 269.23414 -316.55533       0          ***
    # 
    
#----
      
####Latitudinal gradient S3_Fig_3----

    library(RColorBrewer)
    birds_insurance$nsplog10 <- log10(birds_insurance$Nsp_loc)
      
      a <- ggplot(birds_insurance, aes(x=nsplog10, y=latitude, colour = Inward)) +
        geom_point(alpha=0.5,size=0.3)+theme_minimal()+
        scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 8, name = "YlGnBu")))+
        geom_smooth(aes(x = nsplog10),orientation = "y",method="gam",linetype="solid",linewidth=1.5,colour="#F78B8B",alpha=0.8)+
        theme_bw()+
        theme(legend.position = "none",panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank())+
        ylim(-75,85)+
        scale_y_continuous(breaks=seq(-80,90,20))+
        xlab("log10(richness)")+
        ylab("Lattitude")+
        geom_hline(yintercept = 0,
                   color = "#757575", linetype = "dashed",linewidth = 0.5)
      
      b <- ggplot(birds_insurance, aes(x=Inward, y=latitude, colour = nsplog10)) +
        geom_point(alpha=0.5,size=0.3)+theme_minimal()+
        scale_colour_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 8, name = "YlGnBu")))+
        geom_smooth(aes(x = Inward),orientation = "y",method="gam",linetype="solid",linewidth=1.5,colour="#F78B8B",alpha=0.8)+
        theme_bw()+
        theme(legend.position = "none",panel.grid.minor = element_blank())+
        ylim(-75,85)+
        scale_y_continuous(breaks=seq(-80,90,20))+
        xlab("Inward Insurance)")+
        ylab("Lattitude")+
        geom_hline(yintercept = 0,
                   color = "#757575", linetype = "dashed",linewidth = 0.5)
      
      c <- ggplot(birds_insurance, aes(x=Outward, y=latitude, colour = nsplog10)) +
        geom_point(alpha=0.5,size=1)+theme_minimal()+
        scale_colour_gradientn(colours = rev(brewer.pal(n = 8, name = "YlGnBu")))+
        geom_smooth(aes(x = Outward),orientation = "y",method="gam",linetype="solid",linewidth=1.5,colour="#F78B8B",alpha=0.8)+
        theme_bw()+
        theme(legend.position = "none")+
        ylim(-75,85)+
        scale_y_continuous(breaks=seq(-80,90,20))+
        xlab("Outward Insurance)")+
        ylab("Lattitude")+
        geom_hline(yintercept = 0,
                   color = "#757575", linetype = "dashed",linewidth = 0.5)
      
      S3_Fig_3 <- gridExtra::arrangeGrob(b,a,c, ncol=3)
      ggsave(file=here::here("outputs","S3_Fig_3.tiff"), S3_Fig_3,width = 20, height = 7, dpi = 300, units = "cm", device='tiff') 
      
#----

###HI Fig_6_a, S3_Fig_4----
  #Data from Mu et al. 2018 https://doi.org/10.1038/s41597-022-01284-8
  library(ggplot2)    
  library(terra)  
      
  load(here::here('results','birds_insurance.RData'))
      
  # Load the raster with cell IDs
      grd <- suppressWarnings(
        terra::rast(here::here("data", "Loiseau2021", "grid_area_with_cell_ids.tif"))
      )
      
  # Load the raster HI (from https://www.nature.com/articles/s41597-022-01284-8)
      hi <- suppressWarnings(
        terra::rast(here::here("data", "BIG_FILES", "hfp2018.tif"))
      )
      
  #have a look 
      terra::plot(hi, main = "Human Footprint 2018")
      
  # re-project hi to match grd #LONG
      hi_reprojected <- terra::project(hi, terra::crs(grd)) 
      terra::writeRaster(hi_reprojected, filename = here::here("data", "BIG_FILES", "hi_reprojected.tif"), filetype = "GTiff", overwrite = TRUE)
      
  #load hi_reprojected
      hi_reprojected <- suppressWarnings(
        terra::rast(here::here("data", "BIG_FILES", "hi_reprojected.tif"))
      )
      
  # Calculate the mean value of hi_resampled within each cell of grd
      hi_resampled <- terra::resample(hi_reprojected, grd, method = "bilinear")
      mean_values <- zonal(hi_resampled, grd, fun = "mean", na.rm = TRUE)
      
  #S3_Fig_4 
      mean_raster <- terra::rast(grd)
      mean_raster[] <- NA_real_  # Ensure it's a numeric raster
      
      idx <- match(grd[], mean_values$grid_area_with_cell_ids)
      mean_raster[] <- mean_values$hfp2018[idx]
      
      terra::plot(mean_raster, main = "Mean Human Footprint per 50x50 km Cell")
      
      # Reproject the raster to Mollweide and plot 
      mean_raster_moll <- terra::project(mean_raster, "ESRI:54009")
      
      tiff(filename = "outputs/S3_Fig_4.tiff", width = 2000, height = 1300, res = 300)
      terra::plot(mean_raster_moll)
      if (dev.cur() != 1) dev.off()
      
  #merge with bird_insurance and plot 
      
      colnames(mean_values) <- c("cell_focal","hi")
      birds_insurance_hi <- merge(birds_insurance,mean_values,all = TRUE)
      birds_insurance_hi <- birds_insurance_hi[!is.na(birds_insurance_hi$SS_status),]
      
      Fig_6a <- ggplot(birds_insurance_hi, aes(x = SS_status, y = hi, fill = SS_status)) +
        geom_boxplot(outlier.shape = NA) +
        # Add count labels
        # geom_text(data = box_labels, 
        #           aes(x = SS_status, y = med, label = paste0("n=", n)), 
        #           inherit.aes = FALSE, size = 4, color = "black")+
        ylab("Human footprint") +
        scale_fill_manual(values = c("Source" = "#7DBD80", 
                                     "Intermediate" = "#FFF7BA", 
                                     "Sink" = "#FA6C67")) +
        ylim(0,30)+
        theme_bw() +
        theme(
          axis.text.x = element_text(size = 12),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 8),
          axis.title.y = element_text(size = 12),
          legend.position = "none"
        )
      
      kruskal.test(hi ~ SS_status, data = birds_insurance_hi)
      FSA::dunnTest(birds_insurance_hi$hi ~ birds_insurance_hi$SS_status, method="bonferroni")    
      
      save(birds_insurance_hi,file=here::here("results","birds_insurance_hi.RData"))
#---- 
          
###PA Fig_6_b and Fig_S3_4----
library(sf)  
library(terra)  
library(exactextractr)  
library(dplyr)
library(ggplot2)
      
  load(here::here('results','birds_insurance_hi.RData'))
      
  # wdpa shapefile dirs
      shapefile_dirs <- c(here::here("data","BIG_FILES","WDPA_Mar2025_Public_shp","WDPA_Mar2025_Public_shp_0"), 
                          here::here("data","BIG_FILES","WDPA_Mar2025_Public_shp","WDPA_Mar2025_Public_shp_1"), 
                          here::here("data","BIG_FILES","WDPA_Mar2025_Public_shp","WDPA_Mar2025_Public_shp_2"))
      
  # load the polygons
      polygons_list <- lapply(shapefile_dirs, function(dir) {
        sf::st_read(paste0(dir, "/WDPA_Mar2025_Public_shp-polygons.shp"))  # Nom générique, adaptez si nécessaire
      })
      
  # merge within a single sf object and select only terrestrial PAs
      wdpa_polygons <- do.call(rbind, polygons_list)
      wdpa_polygons_ter <- wdpa_polygons[wdpa_polygons$MARINE%in%c(0,1),]
      
  # subset WPA (very long, wdpa_I_union, wdpa_I_II_union will be saved and must be loaded)
      
      wdpa_I <- wdpa_polygons_ter[wdpa_polygons_ter$IUCN_CAT%in%c("Ia",'Ib'),]
      sf::sf_use_s2(FALSE)
      wdpa_I_union <- sf::st_union(wdpa_I)
      save(wdpa_I_union,file=here::here("results","wdpa_I_union.RData"))
      
      wdpa_I_II <- wdpa_polygons_ter[wdpa_polygons_ter$IUCN_CAT%in%c("Ia",'Ib','II'),]
      sf::sf_use_s2(FALSE)
      wdpa_I_II_union <- sf::st_union(wdpa_I_II)
      save(wdpa_I_II_union,file=here::here("results","wdpa_I_II_union.RData"))
      
      wdpa_I_IV <- wdpa_polygons_ter[wdpa_polygons_ter$IUCN_CAT%in%c("Ia",'Ib','II','III','IV'),]
      sf::sf_use_s2(FALSE)
      wdpa_I_IV_union <- sf::st_union(wdpa_I_IV)
      save(wdpa_I_IV_union,file=here::here("results","wdpa_I_IV_union.RData"))
      
  # Load wdpa_I_II_union and wdpa_I_IV_union
      load(here::here("results","wdpa_I_union.RData"))
      load(here::here("results","wdpa_I_II_union.RData"))
      load(here::here("results","wdpa_I_IV_union.RData"))
      
  # Load the raster with cell IDs
      grd <- suppressWarnings(
        terra::rast(here::here("data", "Loiseau2021", "grid_area_with_cell_ids.tif"))
      )
      
  #Merge birds_insurance with wdpa_I_union
      wdpa_extract_I <- exactextractr::exact_extract(grd, wdpa_I_union, include_area = TRUE)
      
      wdpa_extract_iucn_I <- do.call(rbind, wdpa_extract_I)
      
      wdpa_extract_iucn_I <- wdpa_extract_iucn_I[,c(1,3)]
      colnames(wdpa_extract_iucn_I) <- c("cell_focal","PA_I")
      
      birds_insurance_hi_pa <- merge(birds_insurance_hi,wdpa_extract_iucn_I,all = TRUE)
      
      birds_insurance_hi_pa$PA_I[is.na(birds_insurance_hi_pa$PA_I)]=0
      birds_insurance_hi_pa$PA_I <- birds_insurance_hi_pa$PA_I*100
      
  #Merge birds_insurance with wdpa_I_II_union
      wdpa_extract_I_II <- exactextractr::exact_extract(grd, wdpa_I_II_union, include_area = TRUE)
      
      wdpa_extract_iucn_I_II <- do.call(rbind, wdpa_extract_I_II)
      
      wdpa_extract_iucn_I_II <- wdpa_extract_iucn_I_II[,c(1,3)]
      colnames(wdpa_extract_iucn_I_II) <- c("cell_focal","PA_I_II")
      
      birds_insurance_hi_pa <- merge(birds_insurance_hi_pa,wdpa_extract_iucn_I_II,all = TRUE)
      
      birds_insurance_hi_pa$PA_I_II[is.na(birds_insurance_hi_pa$PA_I_II)]=0
      birds_insurance_hi_pa$PA_I_II <- birds_insurance_hi_pa$PA_I_II*100
      
  #Merge birds_insurance with wdpa_I_IV_union
      wdpa_extract_I_IV <- exactextractr::exact_extract(grd, wdpa_I_IV_union, include_area = TRUE)
      
      wdpa_extract_iucn_I_IV <- do.call(rbind, wdpa_extract_I_IV)
      
      wdpa_extract_iucn_I_IV <- wdpa_extract_iucn_I_IV[,c(1,3)]
      colnames(wdpa_extract_iucn_I_IV) <- c("cell_focal","PA_I_IV")
      
      birds_insurance_hi_pa <- merge(birds_insurance_hi_pa,wdpa_extract_iucn_I_IV,all = TRUE)
      
      birds_insurance_hi_pa$PA_I_IV[is.na(birds_insurance_hi_pa$PA_I_IV)]=0
      birds_insurance_hi_pa$PA_I_IV <- birds_insurance_hi_pa$PA_I_IV*100
  
  #Remove the few D_insu NA 
      birds_insurance_hi_pa <- birds_insurance_hi_pa[!is.na(birds_insurance_hi_pa$D_insu),] #there is few D_insu NA (less than 400)
      
  #save
      
      birds_insurance_hi_pa$SS_status<- factor(birds_insurance_hi_pa$SS_status, levels = c("Source", "Intermediate", "Sink"))
      save(birds_insurance_hi_pa, file=here::here("results","birds_insurance_hi_pa.RData"))
    
  #load 
      load(here::here("results","birds_insurance_hi_pa.RData"))
      birds_insurance_hi_pa <- birds_insurance_hi_pa[!is.na(birds_insurance_hi_pa$SS_status),]
      
  #Fig_6_b PA I-IV
    #plot the % of cells which have at least 1% of their area within a PA I-IV Fig_6a
      
      ThtPA <- 1
      # Create empty results list
        results <- list()
      
      # Define  conditions
        conditions <- c("n_zero", "plus")
      
      # Loop through each level of SS_status
        for (status in levels(birds_insurance_hi_pa$SS_status)) {
        
        #status <- "Source"
        
        # Subset data for the current status
        subset_data <- birds_insurance_hi_pa[birds_insurance_hi_pa$SS_status == status, ]
        total_n <- nrow(subset_data)
        
        # Count how many values fall into each condition
        n_zero       <- sum(subset_data$PA_I_IV < ThtPA, na.rm = TRUE)
        n_plus   <- sum(subset_data$PA_I_IV >= ThtPA, na.rm = TRUE)
        
        results[[status]] <- data.frame(
          SS_status = status,
          condition = conditions,
          percentage = round(100 * c(n_zero, n_plus) / total_n, 1)
        )
      }
      
      # Combine all results into a single data.frame, subset and plot
        summary_table <- do.call(rbind, results)
        rownames(summary_table) <- NULL
      
        summary_table$SS_status<- factor(summary_table$SS_status, 
                                       levels = c("Source", "Intermediate", "Sink"))
        summary_plus <- subset(summary_table, condition == "plus")
      
        Fig_S3_4e <- ggplot(summary_plus, aes(x = SS_status, y = percentage, fill = SS_status)) +
                    geom_bar(stat = "identity", width = 0.6,color = "black") +
                    # Force labels at fixed y = 10
                    geom_text(aes(label = paste0(percentage, "%")), 
                              y = 10, color = "black", size = 4) +
                    ylab("cells % with at least 1% of PA (I-IV)") +
                    xlab("") +
                    scale_fill_manual(values = c("Source" = "#7DBD80", 
                                                 "Intermediate" = "#FFF7BA", 
                                                 "Sink" = "#FA6C67")) +
                    theme_bw() +
                    theme(
                      axis.text.x = element_text(size = 12),
                      axis.title.x = element_blank(),
                      axis.text.y = element_text(size = 8),
                      axis.title.y = element_text(size = 12),
                      legend.position = "none"
                    )
      
    #for the cell with PA (>1%) plot the % of protection, Fig_5b
      
      sub_birds_insurance_pa <- birds_insurance_hi_pa[birds_insurance_hi_pa$PA_I_IV>=ThtPA,]

      box_labels <- aggregate(PA_I_IV ~ SS_status, data = sub_birds_insurance_pa, 
                              FUN = function(x) c(n = length(x), med = median(x, na.rm = TRUE)))
      box_labels$n <- box_labels$PA_I_IV[, "n"]
      box_labels$med <- box_labels$PA_I_IV[, "med"]+5
      box_labels$PA_I_IV <- NULL  # Remove the matrix column
      
      Fig_6b <- ggplot(sub_birds_insurance_pa, aes(x = SS_status, y = PA_I_IV, fill = SS_status)) +
        geom_boxplot(outlier.shape = NA) +
        # Add count labels
        # geom_text(data = box_labels, 
        #           aes(x = SS_status, y = med, label = paste0("n=", n)), 
        #           inherit.aes = FALSE, size = 4, color = "black")+
        ylab("% of PA (I-IV) within cells") +
        scale_fill_manual(values = c("Source" = "#7DBD80", 
                                     "Intermediate" = "#FFF7BA", 
                                     "Sink" = "#FA6C67")) +
        theme_bw() +
        theme(
          axis.text.x = element_text(size = 12),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 8),
          axis.title.y = element_text(size = 12),
          legend.position = "none"
        )
  
    #save the figure 6
      
      Fig_6 <- gridExtra::arrangeGrob(Fig_6a,Fig_6b, ncol=2)
      ggsave(file=here::here("outputs","Fig_6.tiff"), Fig_6,width = 20, height = 8, dpi = 300, units = "cm", device='tiff') 
      
    #mean comparison
      kruskal.test(PA_I_IV ~ SS_status, data = birds_insurance_hi_pa[birds_insurance_hi_pa$PA_I_IV>=1,])
      FSA::dunnTest(birds_insurance_hi_pa[birds_insurance_hi_pa$PA_I_IV>=1,]$PA_I_IV ~ birds_insurance_hi_pa[birds_insurance_hi_pa$PA_I_IV>=1,]$SS_status, method="bonferroni")    
      
  #Fig_S3_5a,b PA I
    #plot the % of cells which have at least 1% of their area within a PA I
      
      ThtPA <- 1
      # Create empty results list
      results <- list()
      
      # Define  conditions
      conditions <- c("n_zero", "plus")
      
      # Loop through each level of SS_status
      for (status in levels(birds_insurance_hi_pa$SS_status)) {
        
        #status <- "Source"
        
        # Subset data for the current status
        subset_data <- birds_insurance_hi_pa[birds_insurance_hi_pa$SS_status == status, ]
        total_n <- nrow(subset_data)
        
        # Count how many values fall into each condition
        n_zero       <- sum(subset_data$PA_I < ThtPA, na.rm = TRUE)
        n_plus   <- sum(subset_data$PA_I >= ThtPA, na.rm = TRUE)
        
        results[[status]] <- data.frame(
          SS_status = status,
          condition = conditions,
          percentage = round(100 * c(n_zero, n_plus) / total_n, 1)
        )
      }
      
      # Combine all results into a single data.frame, subset and plot
      summary_table <- do.call(rbind, results)
      rownames(summary_table) <- NULL
      
      summary_table$SS_status<- factor(summary_table$SS_status, 
                                       levels = c("Source", "Intermediate", "Sink"))
      summary_plus <- subset(summary_table, condition == "plus")
      
      Fig_S3_5a <- ggplot(summary_plus, aes(x = SS_status, y = percentage, fill = SS_status)) +
        geom_bar(stat = "identity", width = 0.6,color = "black") +
        # Force labels at fixed y = 10
        geom_text(aes(label = paste0(percentage, "%")), 
                  y = 5, color = "black", size = 4) +
        ylab("cells % with at least 1% of PA (Ia-b)") +
        xlab("") +
        scale_fill_manual(values = c("Source" = "#7DBD80", 
                                     "Intermediate" = "#FFF7BA", 
                                     "Sink" = "#FA6C67")) +
        theme_bw() +
        theme(
          axis.text.x = element_text(size = 12),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 8),
          axis.title.y = element_text(size = 12),
          legend.position = "none"
        )
      
      #for the cell with PA (>1%) plot the % of protection, Fig_S3_4b
      
      sub_birds_insurance_pa <- birds_insurance_hi_pa[birds_insurance_hi_pa$PA_I>=ThtPA,]
      
      box_labels <- aggregate(PA_I ~ SS_status, data = sub_birds_insurance_pa, 
                              FUN = function(x) c(n = length(x), med = median(x, na.rm = TRUE)))
      box_labels$n <- box_labels$PA_I[, "n"]
      box_labels$med <- box_labels$PA_I[, "med"]+5
      box_labels$PA_I <- NULL  # Remove the matrix column
      
      Fig_S3_5b <- ggplot(sub_birds_insurance_pa, aes(x = SS_status, y = PA_I, fill = SS_status)) +
        geom_boxplot(outlier.shape = NA) +
        # Add count labels
        # geom_text(data = box_labels, 
        #           aes(x = SS_status, y = med, label = paste0("n=", n)), 
        #           inherit.aes = FALSE, size = 4, color = "black")+
        ylab("% of PA (Ia-b) within cells") +
        scale_fill_manual(values = c("Source" = "#7DBD80", 
                                     "Intermediate" = "#FFF7BA", 
                                     "Sink" = "#FA6C67")) +
        theme_bw() +
        theme(
          axis.text.x = element_text(size = 12),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 8),
          axis.title.y = element_text(size = 12),
          legend.position = "none"
        )
      
      #mean comparison
      kruskal.test(PA_I ~ SS_status, data = birds_insurance_hi_pa[birds_insurance_hi_pa$PA_I>=1,])
      FSA::dunnTest(birds_insurance_hi_pa[birds_insurance_hi_pa$PA_I>=1,]$PA_I ~ birds_insurance_hi_pa[birds_insurance_hi_pa$PA_I>=1,]$SS_status, method="bonferroni")    
      
      
  #Fig_S3_5c,d PA I-II
      #plot the % of cells which have at least 1% of their area within a PA I-II Fig_5a
      
      ThtPA <- 1
      # Create empty results list
      results <- list()
      
      # Define  conditions
      conditions <- c("n_zero", "plus")
      
      # Loop through each level of SS_status
      for (status in levels(birds_insurance_hi_pa$SS_status)) {
        
        #status <- "Source"
        
        # Subset data for the current status
        subset_data <- birds_insurance_hi_pa[birds_insurance_hi_pa$SS_status == status, ]
        total_n <- nrow(subset_data)
        
        # Count how many values fall into each condition
        n_zero       <- sum(subset_data$PA_I_II < ThtPA, na.rm = TRUE)
        n_plus   <- sum(subset_data$PA_I_II >= ThtPA, na.rm = TRUE)
        
        results[[status]] <- data.frame(
          SS_status = status,
          condition = conditions,
          percentage = round(100 * c(n_zero, n_plus) / total_n, 1)
        )
      }
      
      # Combine all results into a single data.frame, subset and plot
      summary_table <- do.call(rbind, results)
      rownames(summary_table) <- NULL
      
      summary_table$SS_status<- factor(summary_table$SS_status, 
                                       levels = c("Source", "Intermediate", "Sink"))
      summary_plus <- subset(summary_table, condition == "plus")
      
      Fig_S3_5c <- ggplot(summary_plus, aes(x = SS_status, y = percentage, fill = SS_status)) +
        geom_bar(stat = "identity", width = 0.6,color = "black") +
        # Force labels at fixed y = 10
        geom_text(aes(label = paste0(percentage, "%")), 
                  y = 12.5, color = "black", size = 4) +
        ylab("cells % with at least 1% of PA (Ia-b, II)") +
        xlab("") +
        scale_fill_manual(values = c("Source" = "#7DBD80", 
                                     "Intermediate" = "#FFF7BA", 
                                     "Sink" = "#FA6C67")) +
        theme_bw() +
        theme(
          axis.text.x = element_text(size = 12),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 8),
          axis.title.y = element_text(size = 12),
          legend.position = "none"
        )
      
      #for the cell with PA (>1%) plot the % of protection, Fig_5b
      
      sub_birds_insurance_pa <- birds_insurance_hi_pa[birds_insurance_hi_pa$PA_I_II>=ThtPA,]
      
      box_labels <- aggregate(PA_I_II ~ SS_status, data = sub_birds_insurance_pa, 
                              FUN = function(x) c(n = length(x), med = median(x, na.rm = TRUE)))
      box_labels$n <- box_labels$PA_I_II[, "n"]
      box_labels$med <- box_labels$PA_I_II[, "med"]+5
      box_labels$PA_I_II <- NULL  # Remove the matrix column
      
      Fig_S3_5d <- ggplot(sub_birds_insurance_pa, aes(x = SS_status, y = PA_I_II, fill = SS_status)) +
        geom_boxplot(outlier.shape = NA) +
        # Add count labels
        # geom_text(data = box_labels, 
        #           aes(x = SS_status, y = med, label = paste0("n=", n)), 
        #           inherit.aes = FALSE, size = 4, color = "black")+
        ylab("% of PA (Ia-b, II) within cells") +
        scale_fill_manual(values = c("Source" = "#7DBD80", 
                                     "Intermediate" = "#FFF7BA", 
                                     "Sink" = "#FA6C67")) +
        theme_bw() +
        theme(
          axis.text.x = element_text(size = 12),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 8),
          axis.title.y = element_text(size = 12),
          legend.position = "none"
        )
      
      #mean comparison
      kruskal.test(PA_I_II ~ SS_status, data = birds_insurance_hi_pa[birds_insurance_hi_pa$PA_I_II>=1,])
      FSA::dunnTest(birds_insurance_hi_pa[birds_insurance_hi_pa$PA_I_II>=1,]$PA_I_II ~ birds_insurance_hi_pa[birds_insurance_hi_pa$PA_I_II>=1,]$SS_status, method="bonferroni")    
      
      
    #save the figure 
      
      S3_Fig_5 <- gridExtra::arrangeGrob(Fig_S3_5a,Fig_S3_5b,Fig_S3_5c,Fig_S3_5d,Fig_S3_5e,Fig_6b, ncol=2)
      ggsave(file=here::here("outputs","S3_Fig_4.tiff"), S3_Fig_4,width = 20, height = 22, dpi = 300, units = "cm", device='tiff') 
      
      
#----
      

  