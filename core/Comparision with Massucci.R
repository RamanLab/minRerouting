library(readxl)
library(ggplot2)
library(ggsignif)
library(dplyr)
# set wd to where results are saved
#Read in the file from Massucci et al. supplementary data
pcbi_1005949_s001_1_ = read_excel(list.files(pattern = "pcbi"))
pcbi_1005949_s001_1_ = pcbi_1005949_s002
# Read in minRerouting results file
minre = read_excel("one_minReroutingSets.xlsx")
minre = one_minReroutingSets
Rxn1 = which(pcbi_1005949_s001_1_$`React. type` == "Switch")
SL = list()

for (i in 1:(length(Rxn1)-1)){
  SL[[i]] = Rxn1[i+1]- Rxn1[i]
}
SL[[length(Rxn1)]] = nrow(pcbi_1005949_s001_1_) - Rxn1[length(Rxn1)] +1
SL_l = as.data.frame(do.call(rbind, SL))
row.names(SL_l) = t(as.data.frame(pcbi_1005949_s001_1_$Switch[Rxn1]))

flux_diff = list()

for (j in 1:length(Rxn1)){
  zero = 0
  for (i in Rxn1[j]:(SL[[j]]+Rxn1[j])){
    zero = zero + abs(pcbi_1005949_s001_1_[i,9]- pcbi_1005949_s001_1_[i,10])
  }
  flux_diff[[j]] = as.numeric(zero)
}

for ( i in 1:length(SL)){
  SL_l$flux[i] = as.numeric(flux_diff[i])
}

DLs = list()

rxns_sl = list()
for (i in 1:((nrow(minre) - 1)/2)){
  #rxns_sl[[i]] = t(minre[(2*i),c(-1,-2,-3,-4,-5,-6,-7)])
  rxns_sl[[i]] = (minre[(2*i),c(1,2,3,4)])

}

DL_l = do.call(rbind, rxns_sl)
DL_l = as.data.frame(DL_l)
DL_l$flux = as.numeric(DL_l$data34)
DL_l$size = as.numeric(DL_l$data33)
DL_l = group_by(DL_l,data31) %>% summarise(n_count = mean(size), flux = mean(flux))

SL_l$rownames = row.names(SL_l)
SL_l$rownames = gsub("R_", "", SL_l$rownames)

compare = merge(DL_l,SL_l, by.x = "data31", by.y = "rownames")

df1 = compare[,c(1,3)]
df1$type = "Minrerouting"
df2 = compare[,c(1,5)]
df2$type = "Massucci"
colnames(df1) = colnames(df2) = c("Rxn", "Value","type")
d1 = rbind(df1,df2)
d1$Property = "Flux Difference"
  
dff1 = compare[,c(1,2)]
dff1$type = "Minrerouting"
dff2 = compare[,c(1,4)]
dff2$type = "Massucci"
colnames(dff1) = colnames(dff2) = c("Rxn", "Value", "type")
d2 = rbind(dff1,dff2)
d2$Property = "Cluster Size"

df = rbind(d1,d2)

anno = as.data.frame(c( "Cluster Size","Flux Difference"))
for(i in 1:nrow(anno)){

  sub = df[which(df$Property == anno[i,1]),]
  anno[i,2] = t.test(as.numeric(Value) ~as.factor(type), data= sub)$p.value
}

anno$stars = stars.pval(anno$V2)
anno$stars = gsub(" ", "NS", anno$stars)

annotation_df = as.data.frame(c("Cluster Size", "Flux Difference"))
colnames(annotation_df)[1] = "Property"
annotation_df$start = "Minrerouting"
annotation_df$end = "Massucci"
annotation_df$y = c(500,17000)
blank_data <- data.frame(type = c("Cluster Size", "Cluster Size","Flux Difference", "Flux Difference"), x = c("Minrerouting","Massucci"), y = c(0, 
                                                                                                                                                600, 0,19000))


ggplot(df, aes(x = type,y = Value, fill = Property)) + 
  geom_boxplot(aes(fill=Property))+
  geom_signif(data = annotation_df,
              aes(xmin = start, xmax = end, annotations = anno$stars, y_position = y),
              textsize = 10, vjust = 0.0001,
              manual = TRUE)+
  ggtitle(bquote("Synthetic Lethal Comparison for "~italic("S. sonnei")))+
  facet_wrap(~Property,scales="free")+
  theme(axis.title = element_text(face="bold"), plot.title = element_text(face = "bold"))+
  xlab("Data")+
  theme(title = element_text(size = 24))+
  theme(strip.text.x = element_text( margin = margin( b = 0, t = 0) ) )+
  theme_minimal(base_size = 22)+
  theme(panel.border=element_rect(fill=NA, colour="grey40"),text =  element_text(size = 20, face = "bold"))+
  theme(legend.key.size = unit(1,"cm"))

