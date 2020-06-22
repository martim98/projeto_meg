## dados chips ##

file_name<- "data/chip1.txt"

### packages ###

library(genArise)

## read.spot para criar um objeto, como diz no pdf desta bibiloteca ##
## http://www.bioconductor.org/packages/release/bioc/vignettes/genArise/inst/doc/genArise.pdf##
## O Ven é o canal vermelho, logo será o cy5

spot <- read.spot(file_name, cy3 = 2, cy5 = 3,
                   bg.cy3 = 4, bg.cy5 = 5, ids = 1, header = T)

## extrair spot data ##
spot_data <- attr(spot, "spotData")

spot_data

## graficos de intensidade sem normalização
plot(spot_data$Cy3, spot_data$Cy5, pch = 19, main = 'Raw Cy3 vs Cy5',
     xlab = "Cy3", ylab = "Cy5")



## ma plot, ri plot e o gráfico da genArise dos das intensidades##

ma.plot(spot)
ri.plot(spot)
cys.plot(spot)


## correção de background de acordo com a genArise ##
new_backg <- bg.correct(spot)
cys.plot(new_backg)

par(mfrow = c(1, 2))
cys.plot(spot)
cys.plot(new_backg)




## comparação dos dois gráficos
cys.plot(new_backg)
cys.plot(spot)


## comparação sem a correção do background ##
par(mfrow = c(1, 2))
ri.plot(spot)
ri.plot(new_backg)
data(Simon)


## c) 
# transformar os dados com log_2

spot_norm <- spot
spot_norm@spotData$Cy3 <- log(spot_norm@spotData$Cy3, 2)
spot_norm@spotData$Cy5 <- log(spot_norm@spotData$Cy5, 2)
spot_norm@spotData$BgCy3 <- log(spot_norm@spotData$BgCy3, 2)
spot_norm@spotData$BgCy5 <- log(spot_norm@spotData$BgCy5, 2)



# Mesmos gráficos mas normalizados com log_2
spot_data_norm <- attr(spot_norm, "spotData")

plot(spot_data_norm$Cy3, spot_data_norm$Cy5, pch = 19, main = 'Log 2 Cy3 vs Cy5',
     xlab = "Cy3", ylab = "Cy5")


ma.plot(spot_norm)
ri.plot(spot_norm)
cys.plot(spot_norm)


new_backg_norm <- bg.correct(spot_norm)
ri.plot(new_backg_norm)



## calcular z scores com os dados log_2

s.spot_ri <- Zscore(spot_norm,"ri")
s.spot_ma <- Zscore(spot_norm,"ma")

Zscore.plot(s.spot_ri)
Zscore.plot(s.spot_ma)

