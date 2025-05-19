###################################################################################################
#' Spatial Insurance, simulated data 
#'
#'Produces:
#'        
#'        - results/sim.RData dataframe of the simulation runs 
#'        - outputs/FIG_3.tiff (& legend FIG_3_legend.tiff) 
#'        - outputs/S1_Fig_1.tiff
#'        
#' @author Nicolas Mouquet, \email{nicolas.mouquet@@cnrs.fr},
#' @date 2024/04/22 first created, major update 2024/05/10
##################################################################################################

rm(list = ls())
library(ggplot2)

#INIT ----
#set species traits and focal and neighbor community pairs 

#set species traits 
  #We consider 5000 species with three traits T1,T2,T3 draw into a normal distribution  
  
set.seed(1)

  traits <- cbind.data.frame(Species=paste0(rep("SP",5000),1:5000),T1=rnorm(5000,sd = 10),T2=rnorm(5000,sd = 10),
                              T3=rnorm(5000,sd = 10))
  traits <- cbind.data.frame(Species=traits$Species,T1=scale(traits$T1),T2=scale(traits$T2),T3=scale(traits$T3))
  rownames(traits) <- traits$Species
  
#Compute distinctiviness
  
  dist_mat = funrar::compute_dist_matrix(traits[,!colnames(traits) %in% "Species"],metric = "euclidean")
  global_Di = funrar::distinctiveness_global(dist_mat)
  colnames(global_Di) <- c("Species","global_di")
  
  median(dist_mat)
  mean(dist_mat)
  sd(dist_mat)
  boot::boot(dist_mat, function(d, i) median(d[i]), R = 1000)
  
#Compute the pca
  pca_traits <- ade4::dudi.pca(traits[, !colnames(traits) %in% "Species"], scannf = FALSE, nf = 3)
  
# Combine 'Species' column with PCA loadings and global_Di
  pc_axes <- cbind(Species = traits$Species, pca_traits$li)
  pc_axes <- merge(pc_axes, global_Di, by = "Species", all.x = TRUE)
  
#set focal and neighbor community pairs
  #we consider 10000 pairs of communities
  #each community contains 10 to 300 species chosen randomly within the 
  #5000 species pool
  
  set.seed(1)
  occ_list <- pbmcapply::pbmclapply(1:10000,function(id){
    n=sample(10:300,1)
    paste0("SP",sample(1:5000,n))
  },mc.cores = parallel::detectCores())
  names(occ_list) <- 1:10000

#----
  
#FUNCTION INSU----
  
Insu <- function(insured,insurer,occ_list,dist_mat,D_insu,D_thr){
    
    # insured <- Ex[1]
    # insurer <- Ex[2]
    # occ_list <- occ_list
    # dist_mat <- dist_mat
    # D_insu=D_insu
    # D_thr=D_thr
    # 
    spe_insured <- toupper(occ_list[[as.character(insured)]])
    spe_insurer <- toupper(occ_list[[as.character(insurer)]])
    
    #Distance and insurance 
    dist_insured <-dist_mat[colnames(dist_mat)%in%spe_insured,rownames(dist_mat)%in%spe_insured]
    insured_Di = funrar::distinctiveness_global(dist_insured)
    colnames(insured_Di) <- c("Species","Di")
    
    QD <- as.numeric(quantile(insured_Di$Di, probs = seq(0, 1, 0.01))[D_thr+1])
    
    if (D_thr==0) {
      insured_D=insured_Di$Species
      insured_C=insured_Di$Species
    } else {
      insured_D=insured_Di$Species[insured_Di$Di>=QD]
      insured_C=insured_Di$Species[insured_Di$Di<QD]
    }
    
    Insu_D <- do.call(rbind,lapply(insured_D, function(sp_insured){
      #sp_insured <- insured_D[1]
      
      #is the distinct species from the focal community is insured in the neighboring com ? 
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
    
    Insu_C <- do.call(rbind,lapply(insured_C, function(sp_insured){
      #sp_insured <- insured_D[1]
      
      #is the distinct species from the focal community is insured in the neighboring com ? 
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
    
    Inward_D <- round(100*sum(Insu_D$Insurance>0)/nrow(Insu_D),2)
    
    Inward_C <- round(100*sum(Insu_C$Insurance>0)/nrow(Insu_C),2)
    
    rm(dist_insured)
    rm(Insu_D,Insu_C)
    
    cbind.data.frame(D_insu=D_insu,D_thr=D_thr,insured=insured,insurer=insurer,Div_insured=length(spe_insured),
                     Div_insurer=length(spe_insurer),
                     Inward_D=Inward_D,
                     Inward_C=Inward_C,
                     Div_distinct=length(insured_D))
    
}
  
Betafunk <- function(insured,insurer,occ_list,traits){
    
  # insured <- Ex[1]
  # insurer <- Ex[2]
  # occ_list <- occ_list
  # traits <- traits
    
    spe_insured <- occ_list[[insured]]
    spe_insurer <- occ_list[[insurer]]
    
    spelist <- unique(c(spe_insured,spe_insurer))
    spelist <- spelist[order(spelist)]
    
    df <- matrix(ncol = length(spelist), nrow = 2)
    colnames(df) <- spelist
    df[1,] <- ifelse(spelist%in%spe_insured,1,0)
    df[2,] <- ifelse(spelist%in%spe_insurer,1,0)
    rownames(df) <- c("A","B")
    
    traits.sub <-traits[rownames(traits)%in%colnames(df),]
    #traits.sub <- traits.sub[order(rownames(traits.sub)),-1]
    
    traits.sub <- traits.sub[spelist,-1]
    
    colnames(df)
    row.names(traits.sub)
    
    test.pair_funk<-betapart::functional.beta.pair(x=df, traits=traits.sub, index.family = "jaccard")

    rm(traits.sub)
    rm(df)
    
    cbind.data.frame(insured=insured,insurer=insurer,
                     funct.beta.jac=as.numeric(test.pair_funk$funct.beta.jac),
                     funct.beta.jne=as.numeric(test.pair_funk$funct.beta.jne),
                     funct.beta.jtu=as.numeric(test.pair_funk$funct.beta.jtu)
                     )
    
  }
  
#----

####EX WITH TWO CELLS----
library(ggforce)

D_thr <- 90
D_insu <- 0.4

all_div <- do.call(rbind,pbmcapply::pbmclapply(occ_list, function(occ){
  cbind.data.frame(div=length(occ)) 
},mc.cores = parallel::detectCores()))
all_div$id <- rownames(all_div)

insured <- "2"
insurer <- "7"

spe_insured <- occ_list[[insured]]
spe_insurer <- occ_list[[insurer]]

pc_insured <- pc_axes[pc_axes$Species%in%spe_insured,]
pc_insurer <- pc_axes[pc_axes$Species%in%spe_insurer,]

#Compute the localdisctinctivness  

  dist_insured <-dist_mat[colnames(dist_mat)%in%spe_insured,rownames(dist_mat)%in%spe_insured]
  insured_Di = funrar::distinctiveness_global(dist_insured)
  colnames(insured_Di) <- c("Species","insured_Di")
  
  dist_insurer <-dist_mat[colnames(dist_mat)%in%spe_insurer,rownames(dist_mat)%in%spe_insurer]
  insurer_Di = funrar::distinctiveness_global(dist_insurer)
  colnames(insurer_Di) <- c("Species","insurer_Di")
  
  insured_Di_Q <- as.numeric(as.numeric(quantile(insured_Di$insured_Di, probs = seq(0, 1, 0.01))[D_thr+1]))
  insured_spe_Di_Q <- insured_Di$Species[insured_Di$insured_Di>insured_Di_Q]
  pc_insured_Di_Q <- pc_insured[pc_insured$Species%in%insured_spe_Di_Q,]

  insurer_Di_Q <- as.numeric(as.numeric(quantile(insurer_Di$insurer_Di, probs = seq(0, 1, 0.01))[D_thr+1]))
  insurer_spe_Di_Q <- insurer_Di$Species[insurer_Di$insurer_Di>insurer_Di_Q]
  pc_insurer_Di_Q <- pc_insurer[pc_insurer$Species%in%insurer_spe_Di_Q,]

#Project the species on the global PCA 
col_gray <- '#a6a6a510'

fig_pc_insured_1 <- ggplot(pc_axes, aes(x = Axis1, y = Axis2))+
  geom_point(col = col_gray)+
  ylim(min(min(pc_axes$Axis1),min(pc_axes$Axis2)),max(max(pc_axes$Axis1),max(pc_axes$Axis2)))+
  xlim(min(min(pc_axes$Axis1),min(pc_axes$Axis2)),max(max(pc_axes$Axis1),max(pc_axes$Axis2)))+
  theme_bw()+
  geom_point(data = pc_insured, 
             col = "black")+
  geom_point(data = pc_insured_Di_Q, 
             col = "red")

fig_pc_insurer_1 <- ggplot(pc_axes, aes(x = Axis1, y = Axis2))+
  geom_point(col = col_gray)+
  ylim(min(min(pc_axes$Axis1),min(pc_axes$Axis2)),max(max(pc_axes$Axis1),max(pc_axes$Axis2)))+
  xlim(min(min(pc_axes$Axis1),min(pc_axes$Axis2)),max(max(pc_axes$Axis1),max(pc_axes$Axis2)))+
  theme_bw()+
  geom_point(data = pc_insurer, 
             col = "black")+
  geom_point(data = pc_insured_Di_Q,shape=1,
             col = "red")+
  geom_circle(aes(x0 = Axis1, y0 = Axis2, r = D_insu), data = pc_insured_Di_Q,
              col = "#F78C8C",linetype=2)

fig_pc_insured_2 <- ggplot(pc_axes, aes(x = Axis1, y = Axis2))+
  geom_point(col = col_gray)+
  ylim(min(min(pc_axes$Axis1),min(pc_axes$Axis2)),max(max(pc_axes$Axis1),max(pc_axes$Axis2)))+
  xlim(min(min(pc_axes$Axis1),min(pc_axes$Axis2)),max(max(pc_axes$Axis1),max(pc_axes$Axis2)))+
  theme_bw()+
  geom_point(data = pc_insurer, 
             col = "black")+
  geom_point(data = pc_insurer_Di_Q, 
             col = "red")

fig_pc_insurer_2 <- ggplot(pc_axes, aes(x = Axis1, y = Axis2))+
  geom_point(col = col_gray)+
  ylim(min(min(pc_axes$Axis1),min(pc_axes$Axis2)),max(max(pc_axes$Axis1),max(pc_axes$Axis2)))+
  xlim(min(min(pc_axes$Axis1),min(pc_axes$Axis2)),max(max(pc_axes$Axis1),max(pc_axes$Axis2)))+
  theme_bw()+
  geom_point(data = pc_insured, 
             col = "black")+
  geom_point(data = pc_insurer_Di_Q,shape=1,
             col = "red")+
  geom_circle(aes(x0 = Axis1, y0 = Axis2, r = D_insu), data = pc_insurer_Di_Q,
              col = "#F78C8C",linetype=2)

gridExtra::grid.arrange(fig_pc_insured_1,fig_pc_insurer_1,fig_pc_insurer_2,fig_pc_insured_2, ncol=2)
#----

#INWARD OUTWARD PAIRS ALL-----
  
  D_insu_grad <- c(0.4,0.6,0.8)
  Dthr_grad <- c(90)

  set.seed(1)
  dist_mat <- dist_mat
  traits <- traits
  
  sim <- do.call(rbind,pbmcapply::pbmclapply(1:10000,function(i){
    #i=1
    Ex <- sample(names(occ_list),2)
    
    spe_insured <- occ_list[[names(occ_list)[names(occ_list)%in%Ex[1]]]]
    spe_insurer <- occ_list[[names(occ_list)[names(occ_list)%in%Ex[2]]]]
    
    if ((length(spe_insured)>10) & (length(spe_insurer)>10)){
      out <- do.call(rbind,lapply(Dthr_grad,function(Dthr){
        
        do.call(rbind,lapply(D_insu_grad,function(D_insu){
            #D_insu <- D_insu_grad[1]
            a <- Insu(Ex[1],Ex[2],occ_list,dist_mat,D_insu=D_insu,D_thr=Dthr)
            b <- Insu(Ex[2],Ex[1],occ_list,dist_mat,D_insu=D_insu,D_thr=Dthr)
            a$Outward_D <- b$Inward_D
            a$Outward_C <- b$Inward_C
            a
          }))
       
      }))
      beta <- Betafunk(Ex[1],Ex[2],occ_list,traits)
      merge(out,beta)
    }
    
  },mc.cores = parallel::detectCores()-2))
  
  save(sim,file=(here::here('results','sim.RData')))
  
  #Figure final
  library(ggplot2)
  library(RColorBrewer)
  library(gridExtra)
  load(here::here('results','sim.RData'))

  sim_90 <- sim[sim$D_thr%in%90,]
  sim_90$D_insu <- factor(sim_90$D_insu, levels=c('0.8', '0.6', '0.4'))
  
  sim_90_long <- tidyr::pivot_longer(
    data = sim_90,
    cols = c("Inward_D", "Inward_C"),
    names_to = "Inward_D_C",
    values_to = "Inward"
  )

  a <- ggplot(sim_90_long,aes(D_insu,Inward,colors=Inward_D_C))+
    geom_boxplot(aes(fill=Inward_D_C),outlier.shape = NA)+
    ylab("Inward insurance (%)") +
    xlab("Insurance radius") +
    scale_fill_manual(values=c("#D6D6D6","#90DB9B"))+
    theme_bw()+
    theme(legend.position="none")
  
  D_insu=0.6
  D_thr=90
  
  distinct_sim_90 <- sim_90[sim_90$D_insu%in%D_insu,]
  distinct_sim_90 <- distinct_sim_90[distinct_sim_90$D_thr%in%D_thr,]

  distinct_sim_90$diff_div <- distinct_sim_90$Div_insured-distinct_sim_90$Div_insurer
  
  distinct_sim_90 <- distinct_sim_90[order(distinct_sim_90$diff_div,decreasing = F),]
  distinct_sim_90 <- distinct_sim_90[sample(1:nrow(distinct_sim_90),nrow(distinct_sim_90)),]
  
  
  b <- ggplot(distinct_sim_90, aes(x=funct.beta.jac, y=Inward_D, color=Div_insured)) +
    geom_point(alpha=0.7,size = 0.9) +
    geom_density_2d(color = "black", alpha = 0.2, size = 0.3) +  
    ylim(0,100)+
    ylab("Distinct Inward (%)")+xlab("Functional beta diversity")+
    scale_colour_gradientn(colours = rev(brewer.pal(n = 8, name = "RdBu")))+
    theme_bw()+
    theme(legend.position="none")+
    geom_smooth(method =lm,color="#999999")+
    geom_text(x=0.2, 
              y=97, 
              label=   paste0(" r= ",round(cor(distinct_sim_90$funct.beta.jac,distinct_sim_90$Inward_D,use="pairwise.complete.obs"),2)),
              color='#999999',hjust = 0,size=3.5)
  
  c <- ggplot(distinct_sim_90, aes(x=funct.beta.jac, y=Inward_C, color=Div_insured)) +
    geom_point(alpha=0.7,size = 0.9) +
    geom_density_2d(color = "black", alpha = 0.2, size = 0.3) +  
    ylim(0,100)+
    ylab("Common Inward (%)")+xlab("Functional beta diversity")+
    scale_colour_gradientn(colours = rev(brewer.pal(n = 8, name = "RdBu")))+
    theme_bw()+
    theme(legend.position="none")+
    geom_smooth(method =lm,color="#999999")+
    geom_text(x=0.2, 
              y=3, 
              label=   paste0(" r= ",round(cor(distinct_sim_90$funct.beta.jac,distinct_sim_90$Inward_C,use="pairwise.complete.obs"),2)),
              color='#999999',hjust = 0,size=3.5)

  d <- ggplot(distinct_sim_90, aes(x=Inward_D, y=Outward_D, color=Div_insured)) +
    geom_point(alpha=0.7,size = 0.9) +
    ylim(0,100)+xlim(0,100)+
    ylab("Outward (%)")+xlab("Inward (%)")+
    scale_colour_gradientn(colours = rev(brewer.pal(n = 8, name = "RdBu")))+
    theme_bw()+
    theme(legend.position="none")+
    geom_abline (slope=1, linetype = "dashed", color="#999999", size=0.5)
  
  FIG_3 <- gridExtra::arrangeGrob(a,b,c,d, ncol=2)
  ggsave(file=here::here("outputs","FIG_3.tiff"), FIG_3,width = 15, height = 12, dpi = 300, units = "cm", device='tiff') 

  #plot the scale for panel b,c & d to show on the figure 3
  
  p <- ggplot(distinct_sim_90, aes(x = 1, y = 1, color = Div_insured)) +
    geom_point() +
    scale_colour_gradientn(
      colours = rev(brewer.pal(n = 8, name = "RdBu")),
      name = "Nsp"  # Greek Delta
    ) +
    theme_void() +
    theme(legend.position = "right")
  legend <- cowplot::get_legend(p)
  legend <- cowplot::plot_grid(NULL, legend, NULL, 
                                        ncol = 1, 
                                        rel_heights = c(0.3, 0.4, 0.3))
  ggsave(file=here::here("outputs","FIG_3_legend.tiff"), legend,width = 15, height = 15, dpi = 300, units = "cm", device='tiff') 
  
  #S1_Fig_1
  
  a <- ggplot(distinct_sim_90, aes(x=funct.beta.jac, y=Inward_D, color=Div_insured)) +
    geom_point(alpha=0.7) +
    geom_density_2d(color = "black", alpha = 0.2, size = 0.3) +  
    ylim(0,100)+
    ylab("Distinct Inward (%)")+xlab("funct.beta.jac")+
    scale_colour_gradientn(colours = rev(brewer.pal(n = 8, name = "RdBu")))+
    theme_bw()+
    theme(legend.position="none")+
    geom_smooth(method =lm,color="#999999")+
    geom_text(x=0.2, 
              y=97, 
              label=   paste0(" r= ",round(cor(distinct_sim_90$funct.beta.jac,distinct_sim_90$Inward_D,use="pairwise.complete.obs"),2)),
              color='#999999',hjust = 0,size=3.5)
  
  b <- ggplot(distinct_sim_90, aes(x=funct.beta.jne, y=Inward_D, color=Div_insured)) +
    geom_point(alpha=0.7) +
    geom_density_2d(color = "black", alpha = 0.2, size = 0.3) +  
    ylim(0,100)+
    ylab("Distinct Inward (%)")+xlab("funct.beta.jne")+
    scale_colour_gradientn(colours = rev(brewer.pal(n = 8, name = "RdBu")))+
    theme_bw()+
    theme(legend.position="none")+
    geom_smooth(method =lm,color="#999999")+
    geom_text(x=0.2, 
              y=97, 
              label=   paste0(" r= ",round(cor(distinct_sim_90$funct.beta.jne,distinct_sim_90$Inward_D,use="pairwise.complete.obs"),2)),
              color='#999999',hjust = 0,size=3.5)
  
  c <- ggplot(distinct_sim_90, aes(x=funct.beta.jtu, y=Inward_D, color=Div_insured)) +
    geom_point(alpha=0.7) +
    geom_density_2d(color = "black", alpha = 0.2, size = 0.3) +  
    ylim(0,100)+
    ylab("Distinct Inward (%)")+xlab("funct.beta.jtu")+
    scale_colour_gradientn(colours = rev(brewer.pal(n = 8, name = "RdBu")))+
    theme_bw()+
    theme(legend.position="none")+
    geom_smooth(method =lm,color="#999999")+
    geom_text(x=0.2, 
              y=97, 
              label=   paste0(" r= ",round(cor(distinct_sim_90$funct.beta.jtu,distinct_sim_90$Inward_D,use="pairwise.complete.obs"),2)),
              color='#999999',hjust = 0,size=3.5)
  
  
  
  d <- ggplot(distinct_sim_90, aes(x=funct.beta.jac, y=Inward_C, color=Div_insured)) +
    geom_point(alpha=0.7) +
    geom_density_2d(color = "black", alpha = 0.2, size = 0.3) +  
    ylim(0,100)+
    ylab("Common Inward (%)")+xlab("funct.beta.jac")+
    scale_colour_gradientn(colours = rev(brewer.pal(n = 8, name = "RdBu")))+
    theme_bw()+
    theme(legend.position="none")+
    geom_smooth(method =lm,color="#999999")+
    geom_text(x=0.2, 
              y=3, 
              label=   paste0(" r= ",round(cor(distinct_sim_90$funct.beta.jac,distinct_sim_90$Inward_C,use="pairwise.complete.obs"),2)),
              color='#999999',hjust = 0,size=3.5)
  
  e <- ggplot(distinct_sim_90, aes(x=funct.beta.jne, y=Inward_C, color=Div_insured)) +
    geom_point(alpha=0.7) +
    geom_density_2d(color = "black", alpha = 0.2, size = 0.3) +  
    ylim(0,100)+
    ylab("Common Inward (%)")+xlab("funct.beta.jne")+
    scale_colour_gradientn(colours = rev(brewer.pal(n = 8, name = "RdBu")))+
    theme_bw()+
    theme(legend.position="none")+
    geom_smooth(method =lm,color="#999999")+
    geom_text(x=0.2, 
              y=3, 
              label=   paste0(" r= ",round(cor(distinct_sim_90$funct.beta.jne,distinct_sim_90$Inward_C,use="pairwise.complete.obs"),2)),
              color='#999999',hjust = 0,size=3.5)
  
  f <- ggplot(distinct_sim_90, aes(x=funct.beta.jtu, y=Inward_C, color=Div_insured)) +
    geom_point(alpha=0.7) +
    geom_density_2d(color = "black", alpha = 0.2, size = 0.3) +  
    ylim(0,100)+
    ylab("Common Inward (%)")+xlab("funct.beta.jtu")+
    scale_colour_gradientn(colours = rev(brewer.pal(n = 8, name = "RdBu")))+
    theme_bw()+
    theme(legend.position="none")+
    geom_smooth(method =lm,color="#999999")+
    geom_text(x=0.2, 
              y=3, 
              label=   paste0(" r= ",round(cor(distinct_sim_90$funct.beta.jtu,distinct_sim_90$Inward_C,use="pairwise.complete.obs"),2)),
              color='#999999',hjust = 0,size=3.5)
  
  S1_Fig_1 <- gridExtra::arrangeGrob(a,b,c,d,e,f, ncol=3)
  ggsave(file=here::here("outputs","S1_Fig_1.tiff"), S1_Fig_1,width = 20, height = 12, dpi = 300, units = "cm", device='tiff') 
  
  
#----
