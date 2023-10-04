library(stringr)
library(readxl)
library(plyr)
library(tidyverse)

iq <- read_excel("Dados_mundiais.xlsx", sheet = "IQ")
cost <- read_excel("Dados_mundiais.xlsx", sheet = "cost")
life <- read_excel("Dados_mundiais.xlsx", sheet = "life")
paises <- read_csv("paises.csv")
paises <- paises[,c(2,7)]

df <- merge(iq,cost,by="Country")
df <- merge(df, life,by="Country")
df <- merge(df, paises,by.x="Country",by.y="name")

nrow(df)
df[,5]<- as.numeric(str_replace(as.character(df[,5]),"\\.",''))
df[,6] <- as.numeric(str_replace(df[,6], ' °C', ''))
df[,8] <- str_replace(df[,8], ' USD', '')
df[,8] <- as.numeric(str_replace(df[,8], ',', ''))
df[,10] <- as.numeric(str_replace(df[,10], ' years', ''))
df[,11] <- as.numeric(str_replace(df[,11], ' years', ''))
df[,12] <- as.numeric(str_replace(df[,12], ' ‰', ''))
df[,13] <- as.numeric(str_replace(df[,13], ' ‰', ''))




df <- df %>% 
  group_by(continente) %>% 
  mutate(Continente = case_when(
    continente == "África" ~ "África",
    continente == "América" ~ "América",
    continente == "Asia" ~ "Asia",
    continente == "Australia y Oceanía" ~ "Oceania",
    continente == "Europa" ~ "Europa"
  ))

df <- df[,-c(2,14)]
df <- data.frame(df)

df <- df |>
  mutate(IQ = as.numeric(IQ),
         Ø.Income=as.numeric(Ø.Income),
         Education.expenditure..per.inhabitant = as.numeric(Education.expenditure..per.inhabitant),
         Ø.Daily.maximum.temperature = as.numeric(Ø.Daily.maximum..temperature),
         Cost.index = as.numeric(Cost.index),
         Ø.Monthly.income = as.numeric(Ø.Monthly.income),
         Purchasing.power.index = as.numeric(Purchasing.power.index),
         Life.expectancy.males = as.numeric(Life.expectancy.males),
         Life.expectancy.females = as.numeric(Life.expectancy.females),
         Birth.rate = as.numeric(Birth.rate),
         Death.rate = as.numeric(Death.rate)
         )

glimpse(df)

#names(df) <- c("País", "QI", "Renda média", "Despesa educação", "Temperatura máxima média diária",
#               "Índice de custo de vida", "Renda média mensal", "Índice de poder de compra", 
#               "Expectativa de vida masculina", "Expectativa de vida feminina", "Taxa de natalidade", 
#               "Taxa de mortalidade", "Continente")

names(df) <- c("Pais", "QI", "Renda", "Educacao", "TempMax", "Custo", "RendaMensal", "PoderCompra", 
               "ExpectHomens", "ExpectMulheres", "Natalidade", "Mortalidade", "Continente")


View(df)
RendaMensal <- df$RendaMensal
Paises <- df$Pais
df <- df[-c(1,3,6,7,8,length(df))]
df <- data.frame(Paises,RendaMensal,df)
write.csv(df,"dados_mundiais.csv")
