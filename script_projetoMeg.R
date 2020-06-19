## dados chips ##
library(data.table)

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



## ma plot, ri plot e o gráfico da genArise dos das intensidades log_2 ##

ma.plot(spot)
ri.plot(spot)
cys.plot(spot)


## correção de background de acordo com a genArise ##
new_backg <- bg.correct(spot)
ri.plot(new_backg)

## comparação sem a correção do background ##
par(mfrow = c(1, 2))
ri.plot(spot)
ri.plot(new_backg)
data(Simon)

# normalizar os dados com log_2













