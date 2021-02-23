library("Biodem")
library(readxl)
library(markovchain)
install.packages("writexl")
library("writexl")

getwd()
dir()
directorio<-choose.dir()
setwd(directorio)
dir()

df<-read_excel("Dataframe.xlsx")

#######################################################################

#Procesamiento de la data
seq_meli <- df$`secuencia meli`
seq_meli <- seq_meli[is.na(seq_meli) == F]
seq_meli[seq_meli == 3] <- median(seq_meli)

seq_irsa <- df$`secuencia irsa`
seq_irsa <- seq_irsa[is.na(seq_irsa) == F]
seq_irsa[seq_irsa == 3] <- median(seq_irsa)

seq_tenaris <- df$`secuencia tenaris`
seq_tenaris <- seq_tenaris[is.na(seq_tenaris) == F]
seq_tenaris[seq_tenaris == 3] <- median(seq_tenaris)

#######################################################################

#Verifico si las secuencias empíricas cumplen la propiedad de Markov
#Test de independencia Chi cuadrado
verifyMarkovProperty(seq_meli, verbose = TRUE)
verifyMarkovProperty(seq_irsa, verbose = TRUE)
verifyMarkovProperty(seq_tenaris, verbose = TRUE)
#Las tres secuencias cumplen la propiedad de Markov

#######################################################################
#######################################################################
#Se crean tablas de contingencia para hacer el test de homogeneidad

mc_seq<-function(x){createSequenceMatrix(
  stringchar = x,
  toRowProbs = F,
  sanitize = F,
)
}

obs_trans_meli <- split(seq_meli, ceiling(seq_along(seq_meli)/31))
list_meli <- vector('list', length(obs_trans_meli))
for (i in seq_along(obs_trans_meli)){
  list_meli [[i]] <- mc_seq(obs_trans_meli[[i]])
};rm(i)

obs_trans_irsa <- split(seq_irsa, ceiling(seq_along(seq_irsa)/31))
list_irsa <- vector('list', length(obs_trans_irsa))
for (i in seq_along(obs_trans_irsa)){
  list_irsa[[i]] <- mc_seq(obs_trans_irsa[[i]])
};rm(i)

obs_trans_tenaris <- split(seq_tenaris, ceiling(seq_along(seq_tenaris)/31))
list_tenaris <- vector('list', length(obs_trans_tenaris))
for (i in seq_along(obs_trans_tenaris)){
  list_tenaris[[i]] <- mc_seq(obs_trans_tenaris[[i]])
};rm(i)
########################################################################

#Verifico homogeneidad con las tablas de contingencia
verifyHomogeneity(list_meli,verbose = T)
verifyHomogeneity(list_irsa,verbose = T)
verifyHomogeneity(list_tenaris, verbose=T)
#Se verifica homogeneidad para cada acción dada las transiciones observadas
#Observaciones de transición por mes

assessStationarity(seq_meli,1,verbose=T)
assessStationarity(seq_irsa,1,verbose=T)
assessStationarity(seq_tenaris,1,verbose=T)
#########################################################################

#Se estiman las matrices de transición para cada acción
mcf_meli = markovchainFit(
  seq_meli,
  method = "mle",
  byrow = T,
  name = "Markov Chain Meli",
  confidencelevel = 0.95,
  sanitize = FALSE,
)
mct_meli <- mcf_meli$estimate

mcf_irsa = markovchainFit(
  seq_irsa,
  method = "mle",
  byrow = T,
  name = "Markov Chain IRSA",
  confidencelevel = 0.95,
  sanitize = FALSE,
)
mct_irsa <- mcf_irsa$estimate

mcf_tenaris = markovchainFit(
  seq_tenaris,
  method = "mle",
  byrow = T,
  name = "Markov Chain Tenaris",
  confidencelevel = 0.95,
  sanitize = FALSE,
)
mct_tenaris <- mcf_tenaris$estimate

######################################################################

#Resumen de cada cadena de Markov

summary(mct_meli)
summary(mct_irsa)
summary(mct_tenaris)
is.regular(mct_meli)

steadyStates(mct_meli)
steadyStates(mct_irsa)
steadyStates(mct_tenaris)

meanRecurrenceTime(mct_meli)
meanRecurrenceTime(mct_irsa)
meanRecurrenceTime(mct_tenaris)

M_meli <- matrix(c(transitionProbability(mct_meli,t0 = "1",t1 = "1"),transitionProbability(mct_meli,t0 = "1",t1 = "2"),transitionProbability(mct_meli,t0 = "2",t1 = "1"),transitionProbability(mct_meli,t0 = "2",t1 = "2")),ncol=2)
M_irsa <- matrix(c(transitionProbability(mct_irsa,t0 = "1",t1 = "1"),transitionProbability(mct_irsa,t0 = "1",t1 = "2"),transitionProbability(mct_irsa,t0 = "2",t1 = "1"),transitionProbability(mct_irsa,t0 = "2",t1 = "2")),ncol=2)
M_tenaris <- matrix(c(transitionProbability(mct_tenaris,t0 = "1",t1 = "1"),transitionProbability(mct_tenaris,t0 = "1",t1 = "2"),transitionProbability(mct_tenaris,t0 = "2",t1 = "1"),transitionProbability(mct_tenaris,t0 = "2",t1 = "2")),ncol=2)

sim_meli <- markovchainSequence(
  10000,
  mct_meli,
  include.t0 = FALSE,
  useRCpp = TRUE
)

table(sim_meli)/length(sim_meli)
steadyStates(mct_meli)

sim_irsa <- markovchainSequence(
  10000,
  mct_irsa,
  include.t0 = FALSE,
  useRCpp = TRUE
)

table(sim_irsa)/length(sim_irsa)
steadyStates(mct_irsa)

sim_tenaris <- markovchainSequence(
  10000,
  mct_tenaris,
  include.t0 = FALSE,
  useRCpp = TRUE
)

table(sim_tenaris)/length(sim_tenaris)
steadyStates(mct_tenaris)

sim_data <- data.frame(sim_meli,sim_irsa,sim_tenaris)
sim_data = as.data.frame(sapply(sim_data, as.numeric))

write_xlsx(sim_data,"C:\\Users\\agust\\Desktop\\Final estadiística actuarial\\simulaciones.xlsx")


