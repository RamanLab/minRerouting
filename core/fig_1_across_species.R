library(dplyr)
library(ggplot2)
# setwd to results dir
dir = "your working dir"

setwd(paste0(dir, "models" ))
models = c("e_coli_core","iIT341","iML1515", "iPC815" ,"iSSON_1240","iYL1228", "STM_v1_0", "iEK1008")
DL = list()
for (k in 1:length(models)){
  setwd(paste0(dir, "/",models[k]))
  DL[[k]] = as.data.frame(read_xlsx(list.files(pattern = "one_Syn"), sheet = 2)[c(-1,-2),])
}

list = do.call(rbind, DL)
list = list %>% distinct(Var1, Var2)

Abs = list()
for (i in 1:nrow(list)){
  a = 0
  # fin = 0
  for (k in 1:length(check)){
    a = a + sum(DL[[k]][which(DL[[k]]$Var1 == list[i,1]),2] == list[i,2])
  }
  Abs[[i]] = a
}

RI_species = do.call(rbind,Abs)
list$RI_species = RI_species

ggplot(list, aes(x=RI_species[,1])) +
  ggtitle("Frequency Distribution of Synthetic Lethals across Species")+
  geom_histogram(binwidth=0.5)+
  scale_x_continuous(breaks=seq(0,10,1))+
  scale_y_continuous(breaks=seq(0,600,50))+
  xlab("Number of Species")+
  ylab("Number of Synthetic Lethal Pairs")+
  theme_bw()+
  theme(axis.text.x = element_text(face="bold",
                                   size=14),
        axis.text.y = element_text(face="bold",
                                   size=14),
        plot.title = element_text(face = "bold", size = (24)),
        axis.title = element_text(face = "bold", size = 20))

write.csv(list, "Across_species_final.csv")

