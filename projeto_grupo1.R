# installing genArise package

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("genArise")

library(genArise)

# READING DATA FILE AND CREATING SPOT OBJECT

file_name<- "chip1.txt"

#patient1 <- read.dataset(file.name = 'chip1.txt', cy3=3, cy5=2, ids=1, symdesc=NULL, zscore=NULL, type=5, header = TRUE, sep="\t")
# os dados não são do tipo RI nem MA e por isso não conseguem ser lidos pelo read.dataset usa-se read.table que converte em tabela 
patient1_table <- read.table(file_name, header = TRUE, sep="\t", quote = "", dec = ".")
#mas pode se usar read.spot para criar logo classe Spot

patient1 <- read.spot(file_name, cy3=3, cy5=2, bg.cy3=5, bg.cy5=4, ids=1, header = TRUE, sep="\t")
patient1

# a) representar graficamente as intensidades observadas nos dois canais através do diagrama de dispersão e do MA plot

par(mfrow=c(1,2), bg='white')
plot(patient1_table[,3], patient1_table[,2],pch=20, col="dodgerblue") # dispersion
box(col="black")
axis(1, col="black", col.ticks = "black", col.axis="black")
axis(2, col="black", col.ticks = "black", col.axis="black")
mtext("Cy3", side=1, line=3, col="black")
mtext("Cy5", side=2, line=3, col="black")
cys.plot(patient1, col="dodgerblue") #dispersion with log
box(col="black")
axis(1, col="black", col.ticks = "black", col.axis="black")
axis(2, col="black", col.ticks = "black", col.axis="black")
mtext("log2(Cy3)", side=1, line=3, col="black")
mtext("log2(Cy5)", side=2, line=3, col="black")

par(mfrow=c(1,2), bg='white')
ma.plot(patient1)
ri.plot(patient1)

# b) background subtraction

bg_corrected_patient1 <- bg.correct(patient1)

par(mfrow=c(1,2))
cys.plot(bg_corrected_patient1) #dispersion with log
ma.plot(bg_corrected_patient1) #ma plot
#ri.plot(bg_corrected_patient1) #ri plot

# c) normalizing the array:
# non linear (loess) regression of log ratio against average intensity
#1- for each feature construct A = avergae log intensity ; M = log ratio
#2- produce MA plot
#3 - apply loess regression to the data
#4- for each feature calculate NormalizedRatio=FittedLogRatio-RawLogRatio

normalized_bg_patient1 <- global.norm(mySpot = bg_corrected_patient1)

par(mfrow=c(1,2))
cys.plot(normalized_bg_patient1) #dispersion with log
ma.plot(normalized_bg_patient1) #ma plot

# d) Zscore

z_normalized_bg_patient1 <- Zscore(normalized_bg_patient1, type="ma")
par(mfrow=c(1,1))
data(z_normalized_bg_patient1)
Zscore.plot(z_normalized_bg_patient1)


### tabela Zscores
### considerar uma região de rejeição a um alpha 0.05? 
cut_off <- qnorm(0.975)

z_scores <- z_normalized_bg_patient1@dataSets$Zscore
ids <- z_normalized_bg_patient1@dataSets$Id
tabela <- data.frame("Gene_ID" = ids,
                           "Z_score" = z_scores)

tabela_cutoff <- subset(tabela, subset = abs(Z_score) > cut_off)

tabela_cutoff <- tabela_cutoff[order(-tabela_cutoff$Z_score),]
