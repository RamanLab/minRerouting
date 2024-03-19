library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)

models = c("e_coli_core","iIT341","iML1515", "iPC815" ,"iSSON_1240","iYL1228", "STM_v1_0", "iEK1008")
dir = "your resluts directory"
######### PART1: rsl vs psl percentage bar plot ###########

all_df = data.frame()
for (i in 1:length(models)){
  
  setwd(paste0(dir, "results/", models[i]))
  df = read.csv(list.files(pattern = "Compiled"))
  #colnames(df) = c("Rxn_1", "Rxn_2", "Size", "Flux_diff", "common", "SynA", "Rxn_1_class", "Rxn_2_class", "Class","slsize")
  df$species = models[i]
  df$syna = as.numeric(df$syna)/as.numeric(df$size)
  all_df =rbind( all_df, df)
}
all_df = all_df %>% group_by(class, species) %>% summarise(n_count = n())
all_df = all_df %>% group_by(species) %>% mutate(perc = percent(n_count/sum(n_count), 0.1))

ggplot(all_df, aes(fill=class, y=n_count, x=species)) + 
  geom_bar(position="fill", stat="identity")+
  xlab("Species")+
  ylab("Percentage (%)")+
  #coord_flip()+
  theme_bw()+
  theme(legend.text=element_text(size=20), legend.title = element_text(size = 22), axis.title = element_text(size = 20), axis.text = element_text(size = 20 ,face = "bold"))+
  guides(fill = guide_legend(title = "Class"))




######### PART2: rsl psl pie charts with opt or etc ############

all_df = data.frame()
for (i in 1:length(models)){

  setwd(paste0(dir, "results/", models[i]))
  df = read.csv(list.files(pattern = "Compiled"))
  #colnames(df) = c("Rxn_1", "Rxn_2", "Size", "Flux_diff", "common", "SynA", "Rxn_1_class", "Rxn_2_class", "Class","slsize")
  df$species = models[i]
  df$syna = as.numeric(df$syna)/as.numeric(df$size)
  #all_df =rbind( all_df, df)
  all_df = df
  RSL = all_df[(which(all_df$class == "RSL")),]
  PSL = all_df[(which(all_df$class == "PSL")),]
  rsl = RSL %>% group_by(Rxn_1_Class, Rxn_2_Class)
  
  rsl2 = rsl %>% summarise(n_count = n())
  psl = PSL %>% group_by(Rxn_1_Class, Rxn_2_Class)
  psl2 = psl %>% summarise(n_count = n())
  rsl2$Class = paste0(rsl2$Rxn_1_Class, ', ', rsl2$Rxn_2_Class)
  psl2$Class = paste0(psl2$Rxn_1_Class, ', ', psl2$Rxn_2_Class)
  library(dplyr)
  rsl2 = rsl2 %>% mutate(perc = percent(n_count/sum(rsl2$n_count), 0.1))
  psl2 = psl2 %>% mutate(perc = percent(n_count/sum(psl2$n_count), 0.1))
  library(scales)
value = rsl2$n_count/sum(rsl2$n_count)
rsl2$value = (value*100)

df2 <- rsl2 %>%
  mutate(csum = rev(cumsum(rev(value))),
         pos = value/2 + lead(csum, 1),
         pos = if_else(is.na(pos), value/2, pos))

jpeg(file = paste0("RSL_", models[i], ".jpg"), width=800, height=450)
ggplot(rsl2, aes(x = "" , y = value, fill = (Class))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Set3") +
  geom_label_repel(data = df2,
                   aes(y = pos, label = perc),
                   size = 12, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Class")) +
  ggtitle("Distribution of Reaction Classes in RSLs")+
  theme_void()+
  theme(legend.text=element_text(size=20), legend.title = element_text(size = 22), title = element_text(size = 22, face = "bold"))
dev.off()

value = psl2$n_count/sum(psl2$n_count)
psl2$value = (value*100)

df2 <- psl2 %>%
  mutate(csum = rev(cumsum(rev(value))),
         pos = value/2 + lead(csum, 1),
         pos = if_else(is.na(pos), value/2, pos))

jpeg(file=paste0("PSL_", models[i], ".jpg"), width=800, height=450)
ggplot(psl2, aes(x = "" , y = value, fill = (Class))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_brewer(palette = "Dark2") +
  geom_label_repel(data = df2,
                   aes(y = pos, label = perc),
                   size = 12, show.legend = FALSE) +
  guides(fill = guide_legend(title = "Class")) +
  ggtitle("Distribution of Reaction Classes in PSLs")+
  theme_void()+
  theme(legend.text=element_text(size=20), legend.title = element_text(size = 22),title = element_text(size = 22, face = "bold"))
dev.off()
#geom_label_repel()
}

ggplot(rsl2, aes(x="", y=n_count, fill=(Class)))+
  geom_bar(width = 1, stat = "identity")+
  coord_polar("y", start = 0)+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x=element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold"))+
  scale_colour_brewer(palette = "Set2")+ 
  geom_label_repel(aes(label = perc),
                   position = position_stack(vjust = 0.6))

######### PART3: submodule composition bar plot?  ############

df = data.frame()
path_df = data.frame()
comp = data.frame()

for (k in 1:length(models)){

  setwd(paste0(dir,  models[k]))
  SL_dl <- read_excel("one_SyntheticLethal_list.xlsx", sheet = 2)[c(-1,-2),]
  subclass = read_csv("subclasses.csv", col_names = F)
  
  setwd(paste0(dir,"results/",models[k]))
  Combined_everything <- read_csv(list.files(pattern = "Compiled"))

  every = read.csv(list.files(pattern = "Compiled"))
  df_needed = Combined_everything[,c(2,3,25,26)] %>% group_by(Rxn1_submodule,Rxn2_submodule) %>% summarise(count = sum(n()))
  df_needed = na.omit(df_needed)
  df_needed$new_count = 0
  for(i in 1:nrow(df_needed)){
    for (j in 1:i){
      if(df_needed$Rxn1_submodule[i] == df_needed$Rxn2_submodule[j] & i != j &df_needed$Rxn2_submodule[i] == df_needed$Rxn1_submodule[j]){
        #df_needed$Rxn1_submodule[i] = 
        print(j)
        print(i)
        df_needed$new_count[i]= (df_needed$count[i] + df_needed$count[j])
        df_needed = df_needed[-j,]
      }
    }
  }
  
  
  for (i in 1:nrow(df_needed)){
    if(df_needed$new_count[i] != 0){
      df_needed$count[i] = df_needed$new_count[i]
    }
  }

  df_needed$rxns = paste0(df_needed$Rxn1_submodule, ' and ', df_needed$Rxn2_submodule)
  
  df2 = df_needed[which(df_needed$count>mean(df_needed$count)),]
  df2$perc = (df2$count/nrow(every))*100
  
  df2$Type = ifelse(df2$Rxn1_submodule==df2$Rxn2_submodule, "Intra-pathway", "Inter-pathway")

  df2$label =  str_wrap(df2$rxns, width = 40)
  ggplot(df2, aes(x=reorder(label,(count)), y=perc)) +
    geom_segment( aes(x=label, xend=label, y=0, yend=perc, color=Type)) +
    geom_point(aes(color = Type), size=10, alpha=0.6 )+
    theme_light() +
    coord_flip() +
    scale_y_continuous(breaks =seq(0,20,5))+
    theme(
      panel.grid.major.y = element_blank(),
      panel.border = element_blank(),
      axis.ticks.y = element_blank()
    )+
    theme_bw()+
    theme(axis.text  = element_text(size = 20, face = "bold"))+
    xlab("")+
    ylab("Percentage (%)")+
    theme(legend.position = c(.7, .3))+
    theme(legend.text = element_text(size=24), legend.title = element_text(size = 26), legend.box.background = element_rect(colour = "darkgrey"))+
    ggtitle(paste0("Submodule distribution of reaction pairs for model ", models[k]))+
    theme(title = element_text(size = 24, face = "bold"))
  
