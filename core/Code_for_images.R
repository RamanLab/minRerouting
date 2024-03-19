library(readr)
library(dplyr)
library(ggplot2)
library(readxl)
library(ggbreak)
library(ggrepel)
library(gtable)
library(cowplot)
library(grid)
library(ggsignif)
library(gtools)

# Functions needed for generating figures
shift_legend <- function(p){
  
  # check if p is a valid object
  if(!"gtable" %in% class(p)){
    if("ggplot" %in% class(p)){
      gp <- ggplotGrob(p) # convert to grob
    } else {
      message("This is neither a ggplot object nor a grob generated from ggplotGrob. Returning original plot.")
      return(p)
    }
  } else {
    gp <- p
  }
  
  # check for unfilled facet panels
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  if(length(empty.facet.panels) == 0){
    message("There are no unfilled facet panels to shift legend into. Returning original plot.")
    return(p)
  }
  
  # establish extent of unfilled facet panels (including any axis cells in between)
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  empty.facet.panels <- list(min(empty.facet.panels[["t"]]), min(empty.facet.panels[["l"]]),
                             max(empty.facet.panels[["b"]]), max(empty.facet.panels[["r"]]))
  names(empty.facet.panels) <- c("t", "l", "b", "r")
  
  # extract legend & copy over to location of unfilled facet panels
  guide.grob <- which(gp[["layout"]][["name"]] == "guide-box")
  if(length(guide.grob) == 0){
    message("There is no legend present. Returning original plot.")
    return(p)
  }
  gp <- gtable_add_grob(x = gp,
                        grobs = gp[["grobs"]][[guide.grob]],
                        t = empty.facet.panels[["t"]],
                        l = empty.facet.panels[["l"]],
                        b = empty.facet.panels[["b"]],
                        r = empty.facet.panels[["r"]],
                        name = "new-guide-box")
  
  # squash the original guide box's row / column (whichever applicable)
  # & empty its cell
  guide.grob <- gp[["layout"]][guide.grob, ]
  if(guide.grob[["l"]] == guide.grob[["r"]]){
    gp <- gtable_squash_cols(gp, cols = guide.grob[["l"]])
  }
  if(guide.grob[["t"]] == guide.grob[["b"]]){
    gp <- gtable_squash_rows(gp, rows = guide.grob[["t"]])
  }
  gp <- gtable_remove_grobs(gp, "guide-box")
  
  return(gp)
}



### Set working directory to where output of minRerouting was saved
models = c("e_coli_core","iIT341","iML1515", "iPC815" ,"iSSON_1240","iYL1228", "STM_v1_0", "iEK1008")
dir = "/home/kramanlab/Tanisha/examplesss/"

df = data.frame()
one = list()

for (k in 1:length(models)){
  
  setwd(paste0(dir,models[k]))
  
  # Read lists of double lethals and submodules of reactions
  SL_dl <- read_excel("one_SyntheticLethal_list.xlsx", sheet = 2)[c(-1,-2),]
  subclass = read_csv("subclasses.csv", col_names = F)
  
  setwd(paste0(dir,"/",  models[k]))
  one_minre = read_xlsx("one_minReroutingSets.xlsx")
  one = as.data.frame(na.omit(one_minre[-1,1:7]))
  one[,8] = as.numeric(one$data35) + as.numeric(one$data36)
  
  one = merge(one, subclass, by.x = "data31", by.y = "X1")
  one = merge(one, subclass, by.x = "data32", by.y = "X1")
  
  colnames(one) = c("Rxn_1", "Rxn_2", "size", "flux", "short", "long", "common", "syna", "Rxn1_submodule", "Rxn2_submodule")
  
  # Reading the files for pFBA reaction type annotation and classification of SLs into RSL/PSL
  setwd(paste0(dir,"results/",  models[k]))
  rxn_classes = read.csv(list.files(pattern = "pFBA"))
  
  rxn_classes$Rxn_1 = gsub("'", "", rxn_classes$Rxn_1)
  rxn_classes$Rxn_2 = gsub("'", "", rxn_classes$Rxn_2)
  rxn_classes$Rxn_1_Class = gsub("'", "", rxn_classes$Rxn_1_Class)
  rxn_classes$Rxn_2_Class = gsub("'", "", rxn_classes$Rxn_2_Class)
  
  type2 = read.csv(list.files(pattern = "norm"))
  type1 = read.csv(list.files(pattern = "ambi"))
  type3 = merge(type1, rxn_classes,by.x = c("Rxn_1", "Rxn_2"), by.y=c("Rxn_1", "Rxn_2"))
  
  # Using our approach for classification
  RSL_zero = type3[which(type3$Rxn1_Prod > 0 & type3$Rxn2_Prod > 0),]
  RSL_zero = type2[which(type2$Rxn1_Prod > 0 & type2$Rxn2_Prod > 0),]
  type1 = type3[which(type3$Rxn1_Prod <= 0 | type3$Rxn2_Prod <= 0),]
  
  nonNA = na.omit(type1)
  RSL_one = nonNA[which((sign(nonNA$g_min_Fluxes * nonNA$g_max_Fluxes) > 0) & (sign(nonNA$l_min_Fluxes*nonNA$l_max_Fluxes) > 0)) ,]
  NAexist = type1[!(type1$Rxn_1 %in% nonNA$Rxn_1),]
  
  RSL_two = NAexist[which((sign(NAexist$g_min_Fluxes * NAexist$g_max_Fluxes) > 0) | (sign(NAexist$l_min_Fluxes*NAexist$l_max_Fluxes) > 0)) ,]
  PSL_two = type1[which((sign(type1$g_min_Fluxes * type1$g_max_Fluxes) <= 0) | (sign(type1$l_min_Fluxes*type1$l_max_Fluxes) <= 0)) ,]
  
  RSL = rbind(RSL_zero, RSL_one, RSL_two)
  
  PSL = PSL_two
  
  PSL$class = "PSL"
  RSL$class = "RSL"
  class = rbind(RSL,PSL)
  
  compiled = merge(class, one,by.x = c("Rxn_1", "Rxn_2"), by.y=c("Rxn_2", "Rxn_1"))
  # Write the compiled file
  # write.csv(compiled, paste0("Compiled_", models[k], ".csv"))
  
  # Code for finding centrality of each reaction
  DLs = list()
  
  setwd(paste0(dir,models[k]))
  Smat = list.files(pattern = "Smat")
  file = read.csv(Smat)
  rxns = list.files(pattern = "Rxns")
  cnames = read.csv(rxns, header = F)
  names(file) = cnames[,1]
    
  rxns = list()
  for (i in 1:nrow(file)){
    rxns[[i]] = which(file[i,] != 0)
  }
  
  an = list()
  for (j in 1:ncol(file)){
    ant = data.frame();
    for (i in 1:nrow(file)){
      if (j %in% rxns[[i]])
        ant = c(ant , rxns[[i]])
    }
    ant = t(as.data.frame(ant))
    if (nrow(ant) == 0){
      an[j] = 0
    }
    else {
      an[j] = length(unique(ant[,1]))
    }
  }
  
  finnn = do.call(rbind, an)
  row.names(finnn) = cnames[,1] 
  
  minre = read_excel("one_minReroutingSets.xlsx")
    DLs = list()
  for (i in 1:((nrow(minre) - 1)/2)){
    rxns_sl = t(minre[(2*i),c(-1,-2,-3,-4,-5,-6,-7)])
    if (length(finnn[which(row.names(finnn) %in% rxns_sl)]) == 0){
      DLs[[i]] = as.data.frame(NA)
    } else {
      DLs[[i]] = as.data.frame(finnn[which(row.names(finnn) %in% rxns_sl)]/(ncol(file)-1))
      DLs[[i]]$rownames = rxns_sl[1:nrow(DLs[[i]]),]
    }
  }
    rxns_sl = list()
  for (i in 1:((nrow(minre) - 1)/2)){
    rxns_sl[[i]] = t(minre[(2*i),c(1,2)])
  }
  
  
  # For iPC815, delete 176th element from list since the changes in flux are below our cutoff)
  #DLs[[176]]= NULL
  #rxns_sl[[176] = NULL
  
  list_rxn = do.call(cbind, rxns_sl)
  list_rxn = unique(list_rxn)
  
  list_rxn1 = do.call(rbind, rxns_sl)
  list_rxn1= unique(list_rxn1)
  
  RI = list()
  for (i in 1:nrow(list_rxn1)){
    sum = 0
    for (j in 1:ncol(list_rxn)){
      if(length(DLs[[j]])>1){
        a = sum(as.numeric(list_rxn1[i] %in% as.array(list_rxn[1:2,j])))
        sum= sum + a
        sum
      }
    }
    RI[[i]] = sum
  }
  
  RCI = list()
  for (i in 1:nrow(list_rxn1)){
    sum = 0
    for (j in 1:length(DLs)){
      if(length(DLs[[j]])>1){
        a = sum(as.numeric(list_rxn1[i] %in% as.array(DLs[[j]][,2])))
        sum= sum + a
        sum
      }
    }
    RCI[[i]] = sum
  }
  
  RI_ar = as.data.frame(do.call(rbind, RI))/length(DLs)
  row.names(RI_ar) = list_rxn1[,1]
  
  RCI_ar = as.data.frame(do.call(rbind, RCI))/length(DLs)
  row.names(RCI_ar) = list_rxn1[,1]
  RCI_ar$rownames = row.names(RCI_ar)
  
  RI_ar$rownames = row.names(RI_ar)
  r = merge(RCI_ar, RI_ar, by = "rownames")
  
  dlcent = (do.call(rbind, DLs))
  dlcent = dlcent[!duplicated(dlcent[,"rownames"]),]
  row.names(dlcent) = dlcent$rownames
  
 
  setwd(paste0(dir, "/results/",models[k], "/pfba"))
  Ess = read.csv(list.files(pattern = "Ess"), header = F)
  if(k != 1){
    block = read.csv(list.files(pattern = "Blocked"), header = F)
  }
  Opt = read.csv(list.files(pattern = "Opt"),header = F)
  Noflux = read.csv(list.files(pattern = "No"),header = F)
  ELE = read.csv(list.files(pattern = "ELE"),header = F)
  MLE = read.csv(list.files(pattern = "MLE"),header = F)
  
  for (i in 1:nrow(dlcent)){
    dlcent[which(row.names(dlcent) %in% Opt[,1]),3] = "Opt"
    dlcent[which(row.names(dlcent) %in% MLE[,1]),3] = "MLE"
    if(k !=1){
      dlcent[which(row.names(dlcent) %in% block[,1]),3] = "Block"  
    }
    dlcent[which(row.names(dlcent) %in% ELE[,1]),3] = "ELE"
    dlcent[which(row.names(dlcent) %in% Ess[,1]),3] = "Ess"
    dlcent[which(row.names(dlcent) %in% Noflux[,1]),3] = "Noflux"
    
  }
  
  finnnnn = merge(dlcent, r, by.x = "rownames")
  row.names(finnnnn) = finnnnn$rownames
  finnnnn = finnnnn[,c(-1)]
  colnames(finnnnn) = c("Centrality", "Class", "RCI", "RI")
  # Write the csv for centrality and ri,rci values
  # write.csv(finnnnn, "RI_Centrality.csv")
  
}

# Image for ri vs reaction ty[es]
all_df = data_frame()
i = 8
for (i in 1:7){
  setwd(paste0(dir, "results/", models[i], "/pfba/"))
  df = read.csv(list.files(pattern = "_Centrality"))
  
  colnames(df) = c("Rxn", "Centrality", "Class", "RCI", "RI")
  
  setwd(paste0(dir, models[i]))
  trial = read_excel(list.files(pattern = "one_minRerout"))
  ndiv = (nrow(trial)-1)/2
  df$RCI = df$RCI/ndiv
  df$RI = df$RI/ndiv

  df$species = models[i]
  all_df =rbind( all_df, df)
  
}

Centrality_compiled = all_df
p = ggplot(Centrality_compiled,aes(x = Class, y = RI, fill = Class))+
  geom_boxplot(aes(fill = Class))+
  facet_wrap(~species, scales = "free")+
  xlab("Type of Reaction")+
  guides(fill = guide_legend(title = "Type of Reaction"))+
  theme_minimal()+
  ylab("Redundancy Index")+
  ggtitle("Redundancy Index vs Type of Reaction across species")+
  # Add a border around each facet
  theme(panel.border=element_rect(fill=NA, colour="grey40"),text =  element_text(size = 20))+
  theme(text = element_text(size = 24, face = "bold"), legend.text = element_text(size = 22), legend.title = element_text(size = 24),title = element_text(size = 26), legend.key.size = unit(1, 'cm'))

grid.draw(shift_legend(p))

# Image for ri and rci
Centrality = Centrality_compiled %>% group_by(species) %>% mutate(med_ri = median(RI), mean_ri = mean(RI), mean_rci = mean(RCI), med_rci = median(RCI))

df1 = unique(Centrality[,c('species', 'mean_ri')])
df2 = unique(Centrality[,c('species', 'mean_rci')])
df1$type = "Redundancy Index"
df2$type = "Reaction Compensation Index"

colnames(df1)[2] = colnames(df2)[2] = "Mean"
df = rbind(df1,df2)
colnames(df)[3] = "Characteristic"


ggplot(df, aes(x = species, y = Mean, shape = Characteristic, colour = Characteristic))+
  geom_point(size = 10)+
  theme_bw()+
  scale_y_break(c(0.25,0.6), ticklabels = c(0.6,0.65))+
  ylim(0,0.66)+
  xlab("")+
  ylab("Mean Value of Characteristic")+
  ggtitle("Mean Redundany Index and Reaction Compensation Index of Species")+
  theme(text = element_text(size = 24, face = "bold"), title = element_text(size = 28))+
  theme(legend.spacing.y = unit(1, 'cm')) +
  guides(shape = guide_legend(byrow = TRUE))

# Image for synthetic accessibility, flux diff, and size
all_df = data.frame()
for (i in 1:length(models)){
  
  setwd(paste0(dir, models[i]))
  setwd(paste0("/home/kramanlab/Tanisha/examplesss/","results/",  models[i]))
  every = read.csv(list.files(pattern = "Compiled"))
  df = every
  #colnames(df) = c("Rxn_1", "Rxn_2", "Size", "Flux_diff", "common", "SynA", "Rxn_1_class", "Rxn_2_class", "Class","slsize")
  df$species = models[i]
  df$syna = as.numeric(df$syna)/as.numeric(df$size)
  all_df =rbind( all_df, df)
  
}

# Outlier identification:

for (i in 1:length(models)){

  setwd(paste0(dir,"/results/",models[i]))
  Combined_everything <- read.csv(list.files(pattern = "Compiled"))

  out <- boxplot.stats(Combined_everything$flux)$out
  out_ind <- which(Combined_everything$flux %in% c(out))
  
  dat = Combined_everything[out_ind,]
  setwd(paste0(dir,"/results"))
  write.csv(dat,paste0(models[i],"_outliers.csv"))
}


# Finding out the relevant significant interactions
anno <- all_df %>% group_by(species)
anno = anno %>% mutate(pvalue = t.test(
  all_df[class == "PSL", "flux"],
  all_df[class == "RSL", "flux"]
)$p.value)

anno = unique(anno[,c("species","pvalue")])
# Representing the numbers as corresponding asterisks
anno1 = c("NS", "NS","***","***","***","***","***","***")

# Fixing the position and ranges of the significance levels
annotation_df = as.data.frame(models)
colnames(annotation_df)[1] = "species"
annotation_df$start = "PSL"
annotation_df$end = "RSL"
annotation_df$y = c(75,400,500,500,500,500,550,550)
blank_data <- data.frame(species = c("e_coli_core", "e_coli_core","iIT341", "iIT341","iML1515","iML1515","iPC815","iPC815",
                                     "iSSON_1240","iSSON_1240","iYL1228","iYL1228","STM_v1_0","STM_v1_0", "iEK1008", "iEK1008" ), x = c("PSL","RSL"), y = c(-10, 
                                                                                                                                                            100, -10,550,-10,700,-10,600,-10,700,-10,600,-10,700,-10,700))


p = ggplot(all_df, aes(y=as.numeric(size),x = class)) + 
  geom_boxplot(aes(fill=class))+
  geom_signif(data = annotation_df,
              aes(xmin = start, xmax = end, annotations = anno1, y_position = y),
              textsize = 10, vjust = 0.0001,
              manual = TRUE)+
  guides(fill = guide_legend(title = "Class"))+
  xlab("")+
  ylab("SL Cluster Size")+
  theme(strip.text.x = element_text( margin = margin( b = 0, t = 0) ) )+
  geom_blank(data = blank_data, aes(x = x, y = y))+
  expand_limits(y = 0) + scale_y_continuous(expand = c(0, 0))+
  facet_wrap(~species, scales = "free")+
  theme_minimal(base_size = 15) + 
  ggtitle("SL Cluster Size of PSL and RSL Clusters for all species")+
  # Add a border around each facet
  theme(panel.border=element_rect(fill=NA, colour="grey40"),text =  element_text(size = 20))+
  theme(text = element_text(size = 24, face = "bold"), legend.text = element_text(size = 24), legend.title = element_text(size = 26),title = element_text(size = 28), legend.key.size = unit(1.5, 'cm'))

grid.draw(shift_legend(p))

# Code for Fig 7
p = ggplot(all_df,aes(x=syna, y=size, size=flux, color=class)) +
  geom_point(alpha=0.5) +
  scale_size(range = c(4, 12), name="Net Flux Difference")+
  xlab("Synthetic Accessibility")+
  ylab("SL Cluster Size")+
  theme(strip.text.x = element_text( margin = margin( b = 0, t = 0) ) )+
  #geom_blank(data = blank_data, aes(x = x, y = y))+
  #expand_limits(y = 0) + scale_y_continuous(expand = c(0, 0))+
  facet_wrap(~species, scales = "free")+
  theme_minimal(base_size = 15) + 
  ggtitle("SL Cluster Size, Synthetic Accessibility, and Net Flux Difference of Clusters for all species")+
  # Add a border around each facet
  theme(panel.border=element_rect(fill=NA, colour="grey40"),text =  element_text(size = 20))+
  theme(text = element_text(size = 20, face = "bold"), legend.text = element_text(size = 20), legend.title = element_text(size = 22),title = element_text(size = 24), legend.key.size = unit(1, 'cm'))+
  guides(color = guide_legend(title = "Class",override.aes = list(size = 10)))+
  theme(legend.box = "horizontal")

grid.draw(shift_legend(p))

# Read files for outliers and produce gif 8:
all_df = data.frame()
for (i in 2:length(models)){
  setwd(paste0(dir, "results/", models[i]))
  df = read.csv(list.files(pattern = "outliers"))
  #colnames(df) = c("Rxn_1", "Rxn_2", "Size", "Flux_diff", "common", "SynA", "Rxn_1_class", "Rxn_2_class", "Class","slsize")
  df$species = models[i]
  df$syna = as.numeric(df$syna)/as.numeric(df$size)
  all_df =rbind( all_df, df)
  
}

p = ggplot(all_df,aes(x=flux, y=size, size=syna, color=class)) +
  geom_point(alpha=0.5) +
  scale_size(range = c(2, 12), name="Synthetic Accessibility")+
  xlab("Net Flux Difference")+
  ylab("SL Cluster Size")+
  theme(strip.text.x = element_text( margin = margin( b = 0, t = 0) ) )+
  #geom_blank(data = blank_data, aes(x = x, y = y))+
  #expand_limits(y = 0) + scale_y_continuous(expand = c(0, 0))+
  facet_wrap(~species, scales = "free")+
  theme_minimal(base_size = 15) + 
  ggtitle("SL Cluster Size, Synthetic Accessibility, and Net Flux Difference of Clusters for all species")+
  # Add a border around each facet
  theme(panel.border=element_rect(fill=NA, colour="grey40"),text =  element_text(size = 20))+
  theme(text = element_text(size = 20, face = "bold"), legend.text = element_text(size = 20), legend.title = element_text(size = 22),title = element_text(size = 24), legend.key.size = unit(1, 'cm'))+
  guides(color = guide_legend(title = "Class",override.aes = list(size = 10)))+
  theme(legend.box = "horizontal")

grid.draw(shift_legend(p))

