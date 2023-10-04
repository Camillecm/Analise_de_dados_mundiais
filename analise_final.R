library(tidyverse)
library(GGally)
library(corrplot)
library(patchwork)
library(olsrr)
library(gridExtra)
library(gmodels)
library(lindia)
library(qqplotr)
library(dplyr)
library(lmtest)
library(MASS)
library(maps)

df <- read.csv("dados_mundiais.csv")[c(-54,-60),-1]

df <- df|>
  mutate(Continente = as.factor(Continente))

df$Continente = relevel(df$Continente,"Europa")

df1 <- subset(df, df$Continente != "Oceania")

df1<- df1%>%relocate(RendaMensal, .before=Continente)


sqrt(diag(var(df[2:9])))/colMeans(df[2:9],na.rm = TRUE)
summary(df[2:9])
var(df[2:9],na.rm = TRUE)
var(df[2:9])


Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

apply(df[2:9], 2, Mode)



ggpairs(df1[,2:(length(df)-1)], 
        aes(color=df1$Continente),
        lower  = list(continuous = "blank") ,
        diag  = list(continuous = "blankDiag"))

ggpairs(df1[,2:(length(df)-1)])


ggpairs(df[,2:(length(df)-1)],lower  = list(continuous = "blank"),
        upper  = list(continuous = "points") ,
        diag  = list(continuous = "blankDiag"))

df1<- df
df1<- df1 %>%
  relocate(RendaMensal, .before=Continente)
df1[c(3,58),'Continente'] = 'Asia'

boxp <- function(y,name){
  ggplot(df)+
    geom_boxplot(aes(y=df[,y+1],x=Continente,color=Continente))+
    labs(y=name)+
    theme(legend.position="none")+
    scale_color_manual("Continente",breaks=NULL,
                       values=c("#ff6961","#77dd77","#e69c0d",
                                "#84b6f4","#cd9cb2"))+
    scale_x_discrete("",breaks = NULL) +
    theme_bw()
}

boxp1 <- function(y,name){
  ggplot(df)+
    geom_boxplot(aes(y=df[,y+1],x=Continente,color=Continente))+
    labs(y=name)+
    theme(legend.position="none")+
    scale_color_manual("Continente",
                       values=c("#ff6961","#77dd77","#e69c0d",
                                "#84b6f4","#cd9cb2"))+
    scale_x_discrete("",breaks = NULL) +
    theme_bw()
}

boxp(2,"QI") / boxp(3,"Despesas com educação") | 
  boxp(4,"Temperatura máxima diária") /
  boxp(1,"Renda mensal") | boxp(5,"Expectativa de vida masculina") /
  boxp(6,"Expectativa de vida feminina") | boxp1(7,"Taxa de natalidade") /
  boxp1(8,"Taxa de mortalidade")


# Transformações
fit_orig <-  lm(RendaMensal~., data = df1[,-1])
bc <- boxcox(RendaMensal~.,data=df1[,-1],ylab='Log-verossimilhança')
plot(bc)
bc$x[which.max(bc$y)]


# Modelo
fit <- lm(log(RendaMensal)~., data = df1[,-1])
summary(fit)
shapiro.test(fit$residuals)
gqtest(fit)
dwtest(fit)

MyAnova <- function(modelo, nd){
  m1 <- modelo
  np <- dim(anova(m1))[1]
  SQReg <- round(sum(anova(m1)$"Sum Sq"[1:(np-1)]), nd)
  glReg <- sum(anova(m1)$"Df"[1:(np-1)])
  SQRes <- round(anova(m1)$"Sum Sq"[np], nd)
  glRes <- anova(m1)$"Df"[np]
  SQTotal <- round(SQReg + SQRes, nd)
  glTotal <- glReg + glRes
  QMReg <- round(SQReg/glReg, nd)
  QMRes <- round(SQRes/glRes, nd)
  MyF <- round(QMReg/QMRes, nd)
  vpF <- ifelse(pf(MyF, glReg, glRes, lower.tail = F) < 0.0001, "<0.001", 
                roud(pf(MyF, glReg, glRes, lower.tail = F), nd))
  ncolunas <- c("Fonte de Variação", "SQ", "gl", "F", "valor p") 
  Tanova <- data.frame(FV = c("Regressão",
                              "Resíduos",
                              "Total"),
                       gl = c(glReg, glRes, glTotal),
                       SQ = c(SQReg, SQRes, SQTotal),
                       QM = c(QMReg, QMRes, " "),
                       Est.F = c(MyF, " ", " "),
                       valor.p = c(vpF, " ", " ")
  )
  Tanova
}

MyAnova(fit,2)

# selecao de modelos
summary(step(fit))
summary(step(fit, k = log(nrow(df))))

ols_step_best_subset(fit)
ols_step_backward_aic(fit)

all <- data.frame(ols_step_all_possible(fit))
arrange(all,desc(adjr))
head(arrange(all,desc(adjr)))

# Teste F parcial
anova(fit,lm(log(RendaMensal)~QI+Educacao+TempMax+ExpectHomens+ExpectMulheres+
               Natalidade+Mortalidade, data=df1[,-1]))
anova(fit,lm(log(RendaMensal)~QI+Educacao+TempMax+ExpectHomens+ExpectMulheres+
               Natalidade+Continente, data=df1[,-1]))
anova(fit,lm(log(RendaMensal)~QI+Educacao+TempMax+ExpectHomens+ExpectMulheres+
               Mortalidade+Continente, data=df1[,-1]))
anova(fit,lm(log(RendaMensal)~QI+Educacao+TempMax+ExpectHomens+
               Natalidade+Mortalidade+Continente, data=df1[,-1]))
anova(fit,lm(log(RendaMensal)~QI+Educacao+TempMax+ExpectMulheres+
               Natalidade+Mortalidade+Continente, data=df1[,-1]))
anova(fit,lm(log(RendaMensal)~QI+Educacao+ExpectHomens+ExpectMulheres+
               Natalidade+Mortalidade+Continente, data=df1[,-1]))
anova(fit,lm(log(RendaMensal)~QI+TempMax+ExpectHomens+ExpectMulheres+
               Natalidade+Mortalidade+Continente, data=df1[,-1]))
anova(fit,lm(log(RendaMensal)~Educacao+TempMax+ExpectHomens+ExpectMulheres+
               Natalidade+Mortalidade+Continente, data=df1[,-1]))

# sai expect feminina
fit2 <- lm(log(RendaMensal)~QI+Educacao+TempMax+ExpectHomens+ 
           Natalidade+Mortalidade+Continente, data=df[,-1])

anova(fit2,lm(log(RendaMensal)~QI+Educacao+TempMax+ExpectHomens+
               Natalidade+Mortalidade, data=df1[,-1]))
anova(fit2,lm(log(RendaMensal)~QI+Educacao+TempMax+ExpectHomens+
               Natalidade+Continente, data=df1[,-1]))
anova(fit2,lm(log(RendaMensal)~QI+Educacao+TempMax+ExpectHomens+
               Mortalidade+Continente, data=df1[,-1]))
anova(fit2,lm(log(RendaMensal)~QI+Educacao+TempMax+
               Natalidade+Mortalidade+Continente, data=df1[,-1]))
anova(fit2,lm(log(RendaMensal)~QI+Educacao+ExpectHomens+
               Natalidade+Mortalidade+Continente, data=df1[,-1]))
anova(fit2,lm(log(RendaMensal)~QI+TempMax+ExpectHomens+
               Natalidade+Mortalidade+Continente, data=df1[,-1]))
anova(fit2,lm(log(RendaMensal)~Educacao+TempMax+ExpectHomens+
               Natalidade+Mortalidade+Continente, data=df1[,-1]))
# sai temp max
fit3 <- lm(log(RendaMensal)~QI+Educacao+ExpectHomens+ Natalidade
         +Mortalidade+Continente, data=df[,-1])

anova(fit3,lm(log(RendaMensal)~QI+Educacao+ExpectHomens+
                Natalidade+Mortalidade, data=df1[,-1]))
anova(fit3,lm(log(RendaMensal)~QI+Educacao+ExpectHomens+
                Natalidade+Continente, data=df1[,-1]))
anova(fit3,lm(log(RendaMensal)~QI+Educacao+ExpectHomens+
                Mortalidade+Continente, data=df1[,-1]))
anova(fit3,lm(log(RendaMensal)~QI+Educacao+
                Natalidade+Mortalidade+Continente, data=df1[,-1]))
anova(fit3,lm(log(RendaMensal)~QI+ExpectHomens+
                Natalidade+Mortalidade+Continente, data=df1[,-1]))
anova(fit3,lm(log(RendaMensal)~Educacao+ExpectHomens+
                Natalidade+Mortalidade+Continente, data=df1[,-1]))

# sai mortalidade
fit4 <- lm(log(RendaMensal)~QI+Educacao+ExpectHomens+ Natalidade+Continente, 
         data=df[,-1])

anova(fit4,lm(log(RendaMensal)~QI+Educacao+ExpectHomens+
                Natalidade, data=df1[,-1]))
anova(fit4,lm(log(RendaMensal)~QI+Educacao+ExpectHomens+Continente, 
              data=df1[,-1]))
anova(fit4,lm(log(RendaMensal)~QI+Educacao+
                Natalidade+Continente, data=df1[,-1]))
anova(fit4,lm(log(RendaMensal)~QI+ExpectHomens+
                Natalidade+Continente, data=df1[,-1]))
anova(fit4,lm(log(RendaMensal)~Educacao+ExpectHomens+
                Natalidade+Continente, data=df1[,-1])) 

#Sai Natalidade
fit5 <- lm(log(RendaMensal)~QI+Educacao+ExpectHomens+Continente, data=df1[,-1])

anova(fit5,lm(log(RendaMensal)~QI+Educacao+ExpectHomens, data=df1[,-1]))
anova(fit5,lm(log(RendaMensal)~QI+Educacao+Continente, data=df1[,-1]))
anova(fit5,lm(log(RendaMensal)~QI+ExpectHomens+Continente, data=df1[,-1]))
anova(fit5,lm(log(RendaMensal)~Educacao+ExpectHomens+Continente, data=df1[,-1])) 

summary(fit5)


fit_51 <- lm(Natalidade~QI+Educacao+ExpectHomens+Continente, data=df1[,-1])

ggplot()+
  aes(x = fit_51$residuals, y = fit5$residuals) +
  geom_point()+
  geom_smooth(method = lm,se = FALSE) +
  labs(x = "Resíduos Natalidade|.",y="Resíduos log(Renda)|.")+
  theme_bw()

summary(lm(fit5$residuals~fit_51$residuals))


# Teste de interação no modelo
fit6 <- lm(log(RendaMensal)~(QI+Educacao+ExpectHomens)*Continente, data=df1[,-1])
summary(fit6)
anova(lm(log(RendaMensal)~QI+Educacao+ExpectHomens+Continente,data=df1[,-1]),fit6)

# Modelo selecionado
model1 <- lm(log(RendaMensal)~QI+Educacao+ExpectHomens+Continente,data=df1[,-1])
model2<- fit6

# Analise dos residuos
shapiro.test(model1$residuals)
gqtest(model1)
bptest(model1)
dwtest(model1)

shapiro.test(model2$residuals)
gqtest(model2)
bptest(model2)
dwtest(model2)
b <- ggplot(model1)+
  geom_point(aes(x=.fitted,y=studres(model1)))+
  geom_hline(yintercept = c(2,-2), colour = "black",linetype = "dashed")+
  geom_hline(yintercept = c(0), colour = "red",linetype = "dashed")+
  scale_color_discrete("Continente")+
  labs(x="Valores preditos \n A",y="Resíduo Studentizado")+
  theme_bw()

a <- ggplot(model1)+
  aes(sample=studres(model1))+
  geom_qq_line(color="red")+
  #geom_qq_band(fill="lightblue")+
  geom_qq()+
  labs(x="Quantis teóricos \n B", y="Resíduo Studentizado")+
  theme_bw()

b | a


b <- ggplot(model2)+
  geom_point(aes(x=.fitted,y=studres(model2)))+
  geom_hline(yintercept = c(2,-2), colour = "black",linetype = "dashed")+
  geom_hline(yintercept = c(0), colour = "red",linetype = "dashed")+
  scale_color_discrete("Continente")+
  labs(x="Valores preditos \n A",y="Resíduo Studentizado")+
  theme_bw()

a <- ggplot(model2)+
  aes(sample=studres(model2))+
  geom_qq_line(color="red")+
  #geom_qq_band(fill="lightblue")+
  geom_qq()+
  labs(x="Quantis teóricos \n B", y="Resíduo Studentizado")+
  theme_bw()

b | a

# Análise de Multicolinearidade

model <- model1

1-1/vif(model1)

# Analise de observações nao usuais
model <- model1
#leverage
g <- ols_prep_rstudlev_data(model)
d <- g$levrstud
d$txt <- ifelse(d$color == "normal", NA, d$obs)
f <- d[d$color == "outlier", c("obs", "leverage", "rstudent")]
colnames(f) <- c("observation", "leverage", "stud_resid")
a <- ggplot(d, aes(leverage, rstudent, label = txt)) +
  geom_point(aes(colour = fct_color)) +
  scale_color_manual(labels = c("normal","ponto de alavanca","outlier",
                                "outlier e ponto de alavanca"),
                     values = c("black","blue", "red", "green")) +
  xlim(g$minx, g$maxx) +
  ylim(g$miny, g$maxy) +
  labs(colour = "Observação", x = "Leverage", y = "Resíduo estudentizado") +
  geom_hline(yintercept = c(2,-2), colour = "black",linetype = "dashed") +
  geom_vline(xintercept = g$lev_thrsh, colour = "black",linetype = "dashed") +
  geom_text(vjust = -1, size = 3, colour = "black") +
  theme_bw()



# Distancia de Cook
k <- ols_prep_cdplot_data(model)
d <- ols_prep_outlier_obs(k)
f <- ols_prep_cdplot_outliers(k)
b <- ggplot(d, aes(x = obs, y = cd, label = txt)) +
  geom_bar(width = 0.2, stat = "identity",
           aes(fill = fct_color)) +
  scale_fill_manual(values = c("black", "black"),breaks=NULL) +
  ylim(0, k$maxx) +
  labs(x = "Observação", y = "Distância de Cook") +
  geom_hline(yintercept = 0, colour = "white") +
  geom_hline(yintercept = k$ts, colour = "red",linetype = "dashed") +
  geom_text(hjust = -0.2, nudge_x = 0.05, size = 3, na.rm = TRUE, color="red") +
  theme_bw()

b | a


ols_dfbetas <- function (model, print_plot = TRUE) 
{  obs <- NULL
txt <- NULL
dfb <- dfbetas(model)
n <- nrow(dfb)
np <- ncol(dfb)
threshold <- 2/sqrt(n)
myplots <- list()
outliers <- list()
for (i in 2:seqs) {
  dbetas <- dfb[, i]
  df_data <- data.frame(obs = seq_len(n), dbetas = dbetas)
  d <- ols_prep_dfbeta_data(df_data, threshold)
  f <- ols_prep_dfbeta_outliers(d)
  p <- ggplot(d, aes(x = obs, y = dbetas, 
                     label = txt, ymin = 0, ymax = dbetas)) + 
    geom_linerange(colour = "black") + 
    geom_hline(yintercept = c(0, threshold, -threshold), 
               colour = "red") +
    labs(x="",y="") + 
    ggtitle(paste(colnames(dfb)[i])) + 
    theme_bw()+
    geom_text(hjust = -0.2, nudge_x = 0.15, size = 3, 
              colour = "black", na.rm = TRUE) 
  myplots[[i]] <- p
  outliers[[i]] <- f
}
if (print_plot) {
  marrangeGrob(myplots, nrow = 2, ncol = 3, top = quote(paste("")),
               left="DFBETAS",bottom="Observações")
}
}

ols_dfbetas(model)


#DFFITs 
dbetas <- NULL
obs <- NULL
txt <- NULL
dffitsm <- unlist(dffits(model))
k <- length(coef(model))
n <- nrow(df)
dffits_t <- sqrt(k/n) * 2
title <- names(model.frame(model))[1]
dfits_data <- data.frame(obs = seq_len(n), dbetas = dffitsm)
d <- ols_prep_dfbeta_data(dfits_data, dffits_t)
f <- ols_prep_dfbeta_outliers(d)
a <- ggplot(d, aes(x = obs, y = dbetas, label = txt, ymin = 0, 
                   ymax = dffitsm)) + 
  geom_linerange(colour = "black") + 
  geom_hline(yintercept = c(0, dffits_t, -dffits_t), colour = "red") + 
  labs(x = "Observação", y = "DFFITS") +
  theme_bw()+
  geom_text(hjust = -0.2, nudge_x = 0.15, size = 4, colour = "black", 
            na.rm = TRUE) 
#COVRATIOs
cr = data.frame(covratio(model))
cr['obs'] = seq(1,nrow(df))
limite = (3*length(model$coefficients) ) / (nrow(df)-length(model$coefficients))

b <- ggplot(data=cr,aes(x = obs, y = (1-covratio.model.), label = obs))+
  geom_point() +
  geom_hline(yintercept = c(limite, -1*limite),colour = "red")+
  geom_text(aes(label=ifelse(abs(1-covratio.model.)>limite,as.character(obs),'')),
            hjust=-0.2,vjust=0) + 
  labs(x="Observação",y="COVRATIO")+
  theme_bw() 
a | b
