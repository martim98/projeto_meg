# INSTALLING AND IMPORTING genArise PACKAGE

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("genArise")

library(genArise)


# Parte 1 


file_name<- "chip1.txt"
patient1 <- read.spot(file_name, cy3=3, cy5=2, bg.cy3=5, bg.cy5=4, ids=1, header = TRUE, sep="\t")
bg_corrected_patient1 <- bg.correct(patient1)
normalized_bg_patient1 <- global.norm(mySpot = bg_corrected_patient1)
z_normalized_bg_patient1 <- Zscore(normalized_bg_patient1, type="ma")



# e) normalizing the other 2 arrays

patient2 <- read.spot("chip2.txt", cy3=3, cy5=2, bg.cy3=5, bg.cy5=4, ids=1, header = TRUE, sep="\t")
patient3 <- read.spot("chip3.txt", cy3=3, cy5=2, bg.cy3=5, bg.cy5=4, ids=1, header = TRUE, sep="\t")
bg_corrected_patient2 <- bg.correct(patient2)
bg_corrected_patient3 <- bg.correct(patient3)
normalized_bg_patient2 <- global.norm(mySpot = bg_corrected_patient2)
normalized_bg_patient3 <- global.norm(mySpot = bg_corrected_patient3)
par(mfrow=c(1,4))
cys.plot(patient2, "yellow") #dispersion with log
cys.plot(normalized_bg_patient2, "yellow") #dispersion with log
ma.plot(patient2, "yellow") #ma plot
ma.plot(normalized_bg_patient2, "yellow") #ma plot
par(mfrow=c(1,4))
cys.plot(patient3, "yellow") #dispersion with log
cys.plot(normalized_bg_patient3, "yellow") #dispersion with log
ma.plot(patient3, "yellow") #ma plot
ma.plot(normalized_bg_patient3, "yellow") #ma plot


# f) Normalizar entre arrays e boxplots

par(mfrow=c(1,3))

#equivalente a:
#boxplot(log2(patient1@spotData$Cy5/patient1@spotData$Cy3),col=rainbow(6))
#boxplot(log2(patient2@spotData$Cy5/patient2@spotData$Cy3),col=rainbow(6))
#boxplot(log2(patient3@spotData$Cy5/patient3@spotData$Cy3),col=rainbow(6))
M_array <- array(c(log2(patient1@spotData$Cy5/patient1@spotData$Cy3), 
                   log2(patient2@spotData$Cy5/patient2@spotData$Cy3), 
                   log2(patient3@spotData$Cy5/patient3@spotData$Cy3)),
                   dim = c(2994,3))

boxplot(M_array,col=rainbow(6))

M_array_bg_normalized <- array(c(log2(normalized_bg_patient1@spotData$Cy5/normalized_bg_patient1@spotData$Cy3), 
                                 log2(normalized_bg_patient2@spotData$Cy5/normalized_bg_patient2@spotData$Cy3), 
                                 log2(normalized_bg_patient3@spotData$Cy5/normalized_bg_patient3@spotData$Cy3)),
                                 dim = c(2994,3))

boxplot(M_array_bg_normalized,col=rainbow(6))

median_patients=c(0,0,0)
mad_patients=c(0,0,0)
centered_M_array= matrix(0, nrow = 2994, ncol = 3)
for(j in 1:3){
  median_patients[j]= median(M_array[,j])
  mad_patients[j]=mad(M_array[,j])
  for(i in 1:2994){
    centered_M_array[i,j]=(M_array[i,j]-median_patients[j])/mad_patients[j]
  }
}

boxplot(centered_M_array,col=rainbow(6))

# confirmar q apos normalizacao, a media/mediana deve ser 0 e a std/mad 1
median_patients_pos=c(0,0,0)
mad_patients_pos=c(0,0,0)
for(j in 1:3){
  median_patients_pos[j]= median(centered_M_array[,j])
  mad_patients_pos[j]=mad(centered_M_array[,j])
}


# g) Zscore for all patients

z_normalized_bg_patient2 <- Zscore(normalized_bg_patient2, type="ma")
z_normalized_bg_patient3 <- Zscore(normalized_bg_patient3, type="ma")
par(mfrow=c(1,3))
data(z_normalized_bg_patient1)
Zscore.plot(z_normalized_bg_patient1)
data(z_normalized_bg_patient2)
Zscore.plot(z_normalized_bg_patient2)
data(z_normalized_bg_patient3)
Zscore.plot(z_normalized_bg_patient3)

### tabela Zscores
cut_off <- qnorm(0.975)

z_scores <- z_normalized_bg_patient2@dataSets$Zscore
count_2=0
up_count_2=0
for(i in 1:2994){
  if(abs(z_scores[i])>=2) count_2=count_2+1
  if(z_scores[i]>=2) up_count_2=up_count_2+1
}
ids_2 <- z_normalized_bg_patient2@dataSets$Id
tabela_2 <- data.frame("Gene_ID" = ids_2,"Z_score" = z_scores)
tabela_cutoff_2 <- subset(tabela_2, subset = abs(Z_score) > cut_off)
tabela_cutoff_2 <- tabela_cutoff_2[order(-tabela_cutoff_2$Z_score),]

z_scores <- z_normalized_bg_patient3@dataSets$Zscore
count_3=0
up_count_3=0
for(i in 1:2994){
  if(abs(z_scores[i])>=2) count_3=count_3+1
  if(z_scores[i]>=2) up_count_3=up_count_3+1
}
ids_3 <- z_normalized_bg_patient3@dataSets$Id
tabela_3 <- data.frame("Gene_ID" = ids_3,"Z_score" = z_scores)
tabela_cutoff_3 <- subset(tabela_3, subset = abs(Z_score) > cut_off)
tabela_cutoff_3 <- tabela_cutoff_3[order(-tabela_cutoff_3$Z_score),]


# h) Bayesian method for differentially expressed

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")

library(limma)

A_array <- array(c(0.5*log2(patient1@spotData$Cy5*patient1@spotData$Cy3), 
                   0.5*log2(patient2@spotData$Cy5*patient2@spotData$Cy3), 
                   0.5*log2(patient3@spotData$Cy5*patient3@spotData$Cy3)),
                   dim = c(2994,3))

MA_list_patients <- new("MAList", list("M"= M_array, "A"=A_array))
#design <- c(-1,1,-1,1) 
# -1 corresponds to dye-swap
MA_list_fit <- lmFit(MA_list_patients)
MA_list_eBayes <- eBayes(MA_list_fit, proportion = 0.1)
table<-topTable(MA_list_eBayes,number=20,adjust="BH", sort.by="M")
table
write.table(table,"out.txt")
par(mfrow=c(1,1))
volcanoplot(MA_list_eBayes,highlight=20,main="Swirl vs WT")


#eBayes(normalized_bg_patient1, proportion = 0.01, stdev.coef.lim = c(0.1,4), trend = FALSE, robust = FALSE, winsor.tail.p = c(0.05,0.1))
#treat(normalized_bg_patient1, lfc = log2(1.2), trend = FALSE, robust = FALSE, winsor.tail.p = c(0.05,0.1))
