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

## zscores
zscores1 <- Zscore(normalized_bg_patient1, type = 'ma')
zscores2 <- Zscore(normalized_bg_patient2, type = 'ma')
zscores3 <- Zscore(normalized_bg_patient3, type = 'ma')

Zscore.plot(zscores1)
Zscore.plot(zscores2)
Zscore.plot(zscores3)

## Centering com log, como nos slides?
## comentado pois estava a causar problemas com o calculo dos Z scores




normalized_bg_patient1@spotData$Cy3 <- scale(normalized_bg_patient1@spotData$Cy3, 2) 
normalized_bg_patient1@spotData$Cy5 <- scale(normalized_bg_patient1@spotData$Cy5, 2)
#
normalized_bg_patient2@spotData$Cy3 <- scale(normalized_bg_patient2@spotData$Cy3, 2) 
normalized_bg_patient2@spotData$Cy5 <- scale(normalized_bg_patient2@spotData$Cy5, 2)

normalized_bg_patient3@spotData$Cy3<- scale(normalized_bg_patient3@spotData$Cy3, 2) 
normalized_bg_patient3@spotData$Cy5 <- scale(normalized_bg_patient3@spotData$Cy5, 2)



## Package limma precisa de um matrix like data com os cy3 e com os cy5


data <- data.frame('Cy3.1' = normalized_bg_patient1@spotData$Cy3,
                       'Cy3.2' = normalized_bg_patient2@spotData$Cy3,
                       'Cy3.3' = normalized_bg_patient3@spotData$Cy3,
                       'Cy5.1' = normalized_bg_patient1@spotData$Cy5,
                       'Cy5.2' = normalized_bg_patient2@spotData$Cy5,
                       'Cy5.3' = normalized_bg_patient3@spotData$Cy5)

## O design muda com a mudança da cor, logo 3x -1 e 3x 1, de acordo com os slides

fit <- lmFit(data, design = c(-1, -1, -1, 1, 1, 1))
b_fit <- eBayes(fit)
topTable(b_fit, genelist = normalized_bg_patient1@spotData$Id)

## valores p muito altos, algo não está bem ?!

volcanoplot(b_fit)







