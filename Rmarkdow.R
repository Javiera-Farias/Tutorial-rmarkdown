#Rmarkdown
#https://ourcodingclub.github.io/tutorials/rmarkdown/
#Diego González
#04/10/2023

#Packages----
install.packages("permute")
install.packages("indicspecies")
install.packages("readxl")
install.packages("dplyr")
install.packages("tidyverse")
install.packages("car")
install.packages("ggpubr")
install.packages("ggplot2")
install.packages("vegan")
install.packages("kableExtra")
install.packages("glue")
#Library----
library(permute)
library(indicspecies)
library(readxl)
library(dplyr)
library(tidyverse)
library(car)
library(ggpubr)
library(ggplot2)
library(vegan)
library(kableExtra)
library(glue)

#Setwd----
diversidad <- read_excel("C:\\TutorialRmarkdow\\Matriz_modificada.xlsx")
view(diversidad)
choose.files()
#Tablalarga----
data_largo <- diversidad %>% 
  gather(key = "Cuadrante", value = "Abundancia", -ID, -Familia, -PGF)

Tabla <- data.frame(
  data_largo
)

Tabla$ladera <- substr(Tabla$Cuadrante, 1, 1) 

kable(Tabla, caption = "Tabla de datos") %>%
  kable_styling()

#PCoA----
# La funcion tapply se aplica a "Tabla", sobre la matriz abundancia.
#Define los factores (o grupos) en función de los cuales se quiere subdividir Tabla$Abundancia (site,species)
#El resultado (nevados_commat) será una matriz donde las filas representan los diferentes cuadrantes y las columnas representan las diferentes especies (ID). Cada entrada en esta matriz es la suma de abundancias para una combinación específica de cuadrante y especie.

nevados_commat<-tapply(Tabla$Abundancia,INDEX=list(site=Tabla$Cuadrante,species=Tabla$ID),FUN=sum)
dim(nevados_commat)

#calcular las distancias (disimilitudes) entre las filas de nevados_commat usando el método de Bray-Curtis.
distancias <- vegdist(nevados_commat, method = "bray")

# k es el número de dimensiones, usualmente 2 para una visualización 2D
pcoa_resultado <- cmdscale(distancias, eig = TRUE, k = 2) 

# Seleccionar el primer eje de coordenadas para asignar colores
colores_temp <- as.data.frame(pcoa_resultado$points) 

colores <- as.data.frame(substr(rownames(colores_temp), 1, 1)) 

colnames(colores) <- "color"
colores$color2[colores$color=="N"] <- "#d95f02"
colores$color2[colores$color=="O"] <- "#1b9e77" 
colores$color2[colores$color=="S"] <- "#7570b3" 

plot(pcoa_resultado$points[,1], pcoa_resultado$points[,2], xlab = "PCoA1", ylab = "PCoA2", main = "Análisis PCoA", col=colores$color2, pch=16, cex=1.5)
text(pcoa_resultado$points[,1], pcoa_resultado$points[,2], labels = rownames(nevados_commat), cex = 0.7, pos=4)
#PCoA2----
pcoa <- cmdscale(distancias, k=2, eig=TRUE, add=TRUE)

position <- pcoa$points

colnames(position) <- c("pcoa1", "pcoa2")

percent_explained <- 100 * pcoa$eig / sum(pcoa$eig)

round(percent_explained[1:2], digits = 1)
pretty_pe <- format(round(percent_explained[1:2], digits = 1), nsmall=1, trim=TRUE)
labs <- c(glue("PCo 1 ({pretty_pe[1]}%)"),
          glue("PCo 2 ({pretty_pe[2]}%)"))

position
as_tibble(position, rownames="samples")
# Convertir a tibble y extraer solo la primera letra de la columna "samples"
position_tibble <- as_tibble(position, rownames = "samples") %>%
  mutate(samples = str_sub(samples, 1, 1))
ggplot(position_tibble, aes(x=pcoa1, y=pcoa2, color=samples, main = "Analisis PCoA1")) +
  geom_point() +
  labs(x=labs[1], y=labs[2]) +
  theme_minimal() +
  scale_color_manual(values=c(N="#d95f02", S="#7570b3", O="#1b9e77"))
