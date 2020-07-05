# INSTALLING AND IMPORTING genArise PACKAGE

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("genArise")

library(genArise)

# READING DATA FILE AND CREATING SPOT OBJECT

file_name<- "data/chip1.txt"

#patient1 <- read.dataset(file.name = 'chip1.txt', cy3=3, cy5=2, ids=1, symdesc=NULL, zscore=NULL, type=5, header = TRUE, sep="\t")
# os dados n?o s?o do tipo RI nem MA e por isso n?o conseguem ser lidos pelo read.dataset usa-se read.table que converte em tabela 
patient1_table <- read.table(file_name, header = TRUE, sep="\t", quote = "", dec = ".")
#mas pode se usar read.spot para criar logo classe Spot

patient1 <- read.spot(file_name, cy3=3, cy5=2, bg.cy3=5, bg.cy5=4, ids=1, header = TRUE, sep="\t")
patient1

# a) representar graficamente as intensidades observadas nos dois canais atrav?s do diagrama de dispers?o e do MA plot

par(mfrow=c(1,2), bg='black')
plot(patient1_table[,3], patient1_table[,2],pch=20, col="yellow", xlab = "Cy3", ylab = "Cy5") # dispersion
cys.plot(patient1, col="yellow") #dispersion with log

par(mfrow=c(1,2), bg='white')
ma.plot(patient1, "yellow")
ri.plot(patient1, "yellow")

# b) background subtraction


par(mfrow=c(1,3), bg='white')
datos <- attr(patient1, "spotData")# Extract spot data 
M <- log(datos$Cy3, 2) - log(datos$Cy5, 2) 
K <- log(datos$BgCy3, 2) - log(datos$BgCy5, 2)
imageLimma(z = M, row = 499, column = 6, meta.row = 1, meta.column = 1, low = NULL, high = NULL) 
imageLimma(z = K, row = 499, column = 6, meta.row = 1, meta.column = 1, low = NULL, high = NULL) 
imageLimma(z = M-K, row = 499, column = 6, meta.row = 1, meta.column = 1, low = NULL, high = NULL) 

bg_corrected_patient1 <- bg.correct(patient1)

par(mfrow=c(1,2))
cys.plot(bg_corrected_patient1, "yellow") #dispersion with log
ma.plot(bg_corrected_patient1, "yellow") #ma plot
#ri.plot(bg_corrected_patient1) #ri plot

# c) normalizing the array:
# non linear (loess) regression of log ratio against average intensity
#1- for each feature construct A = avergae log intensity ; M = log ratio
#2- produce MA plot
#3 - apply loess regression to the data
#4- for each feature calculate NormalizedRatio=FittedLogRatio-RawLogRatio

normalized_bg_patient1 <- global.norm(mySpot = bg_corrected_patient1)

par(mfrow=c(1,2))
cys.plot(bg_corrected_patient1, "yellow") #dispersion with log
cys.plot(normalized_bg_patient1, "yellow") #dispersion with log
ma.plot(bg_corrected_patient1, "yellow") #ma plot
ma.plot(normalized_bg_patient1, "yellow") #ma plot

# d) Zscore

z_normalized_bg_patient1 <- Zscore(normalized_bg_patient1, type="ma")
par(mfrow=c(1,1))
data(z_normalized_bg_patient1)
Zscore.plot(z_normalized_bg_patient1)

### tabela Zscores
### considerar uma regi?o de rejei??o a um alpha 0.05? 
cut_off <- qnorm(0.975)

z_scores <- z_normalized_bg_patient1@dataSets$Zscore
count=0
for(i in 1:2994){
  if(abs(z_scores[i])>=2) count=count+1
}
  
ids <- z_normalized_bg_patient1@dataSets$Id
tabela <- data.frame("Gene_ID" = ids,
                     "Z_score" = z_scores)

tabela_cutoff <- subset(tabela, subset = abs(Z_score) > cut_off)

tabela_cutoff <- tabela_cutoff[order(-tabela_cutoff$Z_score),]

