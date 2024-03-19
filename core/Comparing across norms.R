library(readxl)
library(ggplot2)
library(dplyr)
library(ggplot2)

#### Set working dir with where your compiled csv is saved
models = c("e_coli_core","iIT341","iML1515", "iPC815" ,"iSSON_1240","iYL1228", "STM_v1_0", "iEK1008")
dir = "your working directory"

two = list()
flux_two = list()
one = list()
flux_one =list()
zero = list()

for (i in 1:length(models)){
  
  setwd(paste0(dir,"/",  models[i]))
  two_minre = read_xlsx(list.files(pattern = "two_minRe"))
  two[[i]] = as.data.frame(na.omit(two_minre[-1,c("data33","data34","data37")]))
  two[[i]][,4] = models[i]
  two[[i]][,5] = "two"

  one_minre = read_xlsx(list.files(pattern = "one_minRe"))
  one[[i]] = as.data.frame(na.omit(one_minre[-1,c("data33","data34","data37")]))
  one[[i]][,4] = models[i]
  one[[i]][,5] = "one"
  
  zero_minre = read_xlsx(list.files(pattern = "zero_minRe"))
  zero[[i]] = as.data.frame(na.omit(zero_minre[-1,c("data33","data34","data37")]))
  zero[[i]][,4] = models[i]
  zero[[i]][,5] = "zero"

}

two = do.call(rbind, two)
one = do.call(rbind, one)
zero = do.call(rbind, zero)
final = rbind(two, one, zero)

colnames(final) = c("size", "flux", "common", "model", "Norm")
final <- final %>% mutate_at(c("size", "flux", "common"), as.numeric)

df = final %>% group_by(model, Norm) %>% mutate(mean_size = mean(size), mean_flux = mean(flux), mean_common = mean(common))

#### Replace variable for y with variable of study (mean_common, mean_size, mean_flux)
ggplot(df, aes(x = model, y = mean_common))+
  geom_point(size = 10, position = position_dodge(width = 0.2), aes(shape = Norm, color= Norm))+
  theme_bw()+
  #guides(shape = guide_legend(override.aes = list(size = 5), title = "Norm"))+
  ggtitle("Mean of common reactions cluster size for synthetic lethal pairs across species")+
  scale_colour_brewer(palette = "Set2")+ 
  xlab("")+
  ylab("Common SL Cluster Size")+
  theme(text = element_text(size=26, face= "bold"), axis.text=element_text(size=22, face= "bold"),
        axis.title=element_text(size=24,face="bold"))+
  scale_shape_manual(values = c(15,17,19))+
  theme(legend.spacing.y = unit(1, 'cm')) +
  guides(color = guide_legend(byrow = TRUE))








  
