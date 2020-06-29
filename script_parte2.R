### parte 2 ### 


library(genArise)
library(limma)
path1 <- "data/chip1.txt"
path2 <- "data/chip2.txt"
path3 <- "data/chip3.txt"

spot1 <- read.spot(path1, cy3 = 3, cy5 = 2,
                  bg.cy3 = 5, bg.cy5 = 4, ids = 1, header = T)
spot2 <- read.spot(path2, cy3 = 3, cy5 = 2,
                   bg.cy3 = 5, bg.cy5 = 4, ids = 1, header = T)
spot3 <- read.spot(path3, cy3 = 3, cy5 = 2,
          bg.cy3 = 5, bg.cy5 = 4, ids = 1, header = T)

## a)
## corrigir o background
bg_corrected_patient1 <- bg.correct(spot1)
bg_corrected_patient2 <- bg.correct(spot2)
bg_corrected_patient3 <- bg.correct(spot3)

## normalizar como no primeiro procedimento
normalized_bg_patient1 <- global.norm(mySpot = bg_corrected_patient1)
normalized_bg_patient2 <- global.norm(mySpot = bg_corrected_patient2)
normalized_bg_patient3 <- global.norm(mySpot = bg_corrected_patient3)

## b) 
## comparar os dados das três amostras
## cy3

boxplot(normalized_bg_patient1@spotData$Cy3, normalized_bg_patient2@spotData$Cy3, 
        normalized_bg_patient3@spotData$Cy3,
        ylim = c(0, 5000),
        xlab = 'Chips',
        ylab = 'Normalized intensities',
        main = 'Cy3 Normalized intensities')
axis(1, at=1:3, labels=c('Chip 1', 'Chip 2', 'Chip 3'))

## cy5

boxplot(normalized_bg_patient1@spotData$Cy5, normalized_bg_patient2@spotData$Cy5, 
        normalized_bg_patient3@spotData$Cy5,
        ylim = c(0, 5000),
        xlab = 'Chips',
        ylab = 'Normalized intensities',
        main = 'Cy5 Normalized intensities')
axis(1, at=1:3, labels=c('Chip 1', 'Chip 2', 'Chip 3'))

## Centering

# normalized_bg_patient1@spotData$Cy3 <- scale(normalized_bg_patient1@spotData$Cy3) 
# normalized_bg_patient1@spotData$Cy5 <- scale(normalized_bg_patient1@spotData$Cy5)
# 
# normalized_bg_patient2@spotData$Cy3 <- scale(normalized_bg_patient2@spotData$Cy3) 
# normalized_bg_patient2@spotData$Cy5 <- scale(normalized_bg_patient2@spotData$Cy5)
# 
# 
# normalized_bg_patient3@spotData$Cy3<- scale(normalized_bg_patient3@spotData$Cy3) 
# normalized_bg_patient3@spotData$Cy5 <- scale(normalized_bg_patient3@spotData$Cy5)

## colocar num data frame
# data_centered <- data.frame(centered_patient1_cy3, 
#                             centered_patient1_cy5,
#                             centered_patient2_cy3,
#                             centered_patient2_cy5,
#                             centered_patient3_cy3,
#                             centered_patient3_cy5)

## Zscore

zscores1 <- Zscore(normalized_bg_patient1, type = 'ma')
zscores2 <- Zscore(normalized_bg_patient2, type = 'ma')
zscores3 <- Zscore(normalized_bg_patient3, type = 'ma')

Zscore.plot(zscores1)
Zscore.plot(zscores2)
Zscore.plot(zscores3)

## Package limma precisa de uma matrix com os cy3 e com os cy5


data_cy3 <- data.frame('Cy3.1' = normalized_bg_patient1@spotData$Cy3,
                       'Cy3.2' = normalized_bg_patient2@spotData$Cy3,
                       'Cy3.3' = normalized_bg_patient3@spotData$Cy3)
data_cy3 <- as.matrix(data_cy3)
data_cy5 <- data.frame('Cy5.1' = normalized_bg_patient1@spotData$Cy5,
                       'Cy5.2' = normalized_bg_patient2@spotData$Cy5,
                       'Cy5.3' = normalized_bg_patient3@spotData$Cy5)
data_cy5 <- as.matrix(data_cy5)


data <- read.maimages(files = c(path1, path2, path3))










obj <- list(M = '', A = '')
obj
class(obj) <- 1
obj$M <- data_cy3
obj$A <- data_cy5



fit <- lmFit(obj)
b_fit <- eBayes(fit)
b_fit

table <- topTable(b_fit,number=20,adjust="BH")
table
