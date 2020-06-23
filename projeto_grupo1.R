# installing genArise package

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("genArise")

library(genArise)

#reading data file and creating class Spot

#patient1 <- read.dataset(file.name = 'chip1.txt', cy3=3, cy5=2, ids=1, symdesc=NULL, zscore=NULL, type=5, header = TRUE, sep="\t")
# os dados não são do tipo RI nem MA e por isso não conseguem ser lidos pelo read.dataset usa-se read.table que converte em tabela 
patient1_table <- read.table("chip1.txt", header = TRUE, sep="\t", quote = "", dec = ".")
#mas pode se usar read.spot para criar logo classe Spot

patient1 <- read.spot("chip1.txt", 3, 2, 5, 4, 1, symdesc=NULL, header = TRUE, sep="\t")
patient1
# representar graficamente as intensidades observadas nos dois canais através do diagrama de dispersão e do MA plot

par(mfrow=c(1,1))
plot(patient1_table[,3], patient1_table[,2], type="p", col="yellow", main ="Diagrama de dispersão", xlab= "Cy3", ylab="Cy5") # dispersion
cys.plot(patient1, col="yellow") #dispersion with log

ma.plot(patient1, col="yellow")

#background subtraction

bg_corrected_patient1 <- bg.correct(patient1)

cys.plot(patient1, col="yellow") #dispersion with log
ri.plot(patient1, col="yellow")
ma.plot(patient1, col="yellow") #ma plot

#normalizing the array:
# non linear (loess) regression of log ratio against average intensity
#1- for each feature construct A = avergae log intensity ; M = log ratio
#2- produce MA plot
#3 - apply loess regression to the data
#4- for each feature calculate NormalizedRatio=FittedLogRatio-RawLogRatio

#m_values <- m.arise(bg_corrected_patient1)
#a_values <- a.arise(bg_corrected_patient1)
#a_values_ <- sort(a_values, decreasing=FALSE)
#MA_dataframe <- data.frame(a_values_, m_values)
#loessMod10 <- loess(m_values ~ a_values_, data=MA_dataframe, span=0.10)
#smoothed10 <- predict(loessMod10)
#lines(x=a_values, y=smoothed10, col="red")

normalized_bg_patient1 <- global.norm(mySpot = bg_corrected_patient1)
data(Simon)
n.spot <- global.norm(Simon)
