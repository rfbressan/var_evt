## Artigo a ser proposto na disciplina MEFCA. Baseado no método de McNeil2000
## Filtar os retornos diários de vários índices primeiramente por um modelo 
## ARMA-GARCH para em seguida ajustar uma distribuição GPD aos resíduos.
## Com esta modelagem é possível estimar os valores de VaR e ES
# Indices de bolsas utilizados
# BVSP - Bovespa Brasil
# MERV - Merval Argentina
# IPSA - IPSA Chile
# MXX - IPC México
# GSPC - SP500 EUA
# GSPTSE - SP/TSX Canada
# Dados entre 31/08/2005 a 31/08/2017

#' # Revisão do artigo para a revista Empirical Economics
#' 
#' ## Ativos
#' 
#' O parecerista recomendou utilizar outros ativos ou índices que não somente
#' de ações. Podemos considerar as seguintes classes de ativos e buscar seus 
#' melhores benchmarks
#' 
#' * Stocks 
#'     - FTSE Global All Cap Index (GEISAC.FGI)
#'     - iShares MSCI World ETF (URTH) 2012-
#'     - MSCI World
#' * Commodities
#'     - MSCI World Commodity Producers Price Return USD
#' * Infrastructure
#'     - MSCI World Infrastructure Price Return USD
#' * Currency pairs
#' * Cryptocurrencies (bitcoin)
#'     - Bitcoin USD (BTC-USD)
#' * Alternative assets
#'     - 
#' * Long-term interest rates
#'     - Treasury Yield 30 Years (^TYX)
#' * Real Estate
#'     - MSCI US REIT
#'     - MSCI Sector Real Estate Indexes
#'     - MSCIUSREIT INDEX (^RMZ)
#'     
#' ## Códigos BBG
#' 
#' MSCI   BBG       Nome
#' 990100	MSDUWI	  MSCI World USD
#' 128456	RMZ	   	  MSCI US REIT INDEX 
#' 126174	MXWO0CMP	World Commodity Producers Price Return USD
#' 129867	MXWO0INF	World Infrastructure Price Return USD

#' 
#' 
#' ## Períodos in-sample e out-of-sample
#' 
#' Evitar terminar o período in-sample exatamente no meio da crise financeira
#' de 2008/09. Pode ser interessante fazer o estudo com 2 backtests diferentes.
#' Um primeiro com período in-sample terminando em 31/12/2007 e backtest indo 
#' até final de 2010. E outro com período in-sample a partir de 01/07/2009 até 
#' 31/12/2013 e out-of-sample a partir daí. 

# Inicio ------------------------------------------------------------------
packages <- c("fExtremes", 
              "rugarch", 
              "PerformanceAnalytics", 
              "xtable", 
              "tidyverse", 
              "broom", 
              "purrr", 
              "gridExtra", 
              "WeightedPortTest")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if (length(new.packages)) install.packages(new.packages)

#library(fGarch)
#library(fBasics)
#library(QRM)
#library(timeSeries)
library(xtable)
library(WeightedPortTest)
#library(CADFtest) # Teste Dickey-Fuller
#library(xts)
library(PerformanceAnalytics)
library(fExtremes)
library(rugarch)
library(tidyverse)
library(broom)
library(gridExtra)
library(kableExtra)
source("./R/artigo_fun.R") # Carrega a funcao roll_fit para fazer o backtest

# AMOSTRA COM DADOS A PARTIR DE 31-08-2000
# 
start <- as.Date("2002-12-31")
end <- as.Date("2008-12-31")
backstart <- end + 1
u_quant <- 0.92 # quantile for treshold u
options(xtable.booktabs = TRUE) # utilizar o pacote booktabs no Latex
options(xtable.tabular.environment = "longtable") # Utiliza o pacote longtable no Latex
options(xtable.floating = FALSE) # Nao pode ser utilizado em conjunto com longtable

list.returns <- function(asset, start) {
  tb <- read_csv(paste0("./input/artigo-", asset, ".csv"), 
                 col_types = cols_only(Date = col_date(), `Adj Close` = col_double()))
  prices <- xts(tb$`Adj Close`, order.by = tb$Date)[paste0(start, "/")]
  colnames(prices) <- "close"
  return(na.omit(Return.calculate(prices, method = "log")))
}
# Gera um tible com uma coluna com o codigo do ativo, a serie de retornos - ts e
# o nome do indice - id_name
assets <- c("BVSP", "GSPC", "GSPTSE", "IPSA", "MERV", "MXX")
lista <- lapply(assets, list.returns, start)
names(lista) <- assets
assets.tbl <- enframe(lista) %>% 
  bind_cols(tibble(id_name = c("IBovespa", "S&P500", "S&P TSE", "IPSA", "Merval", "IPC")))

colnames(assets.tbl) <- c("indice", "ts", "id_name")
# Qual o tamanho da janela in sample?
insample <-  enframe(map_dbl(lista, ~ndays(.x[paste0("/", end)])))
names(insample) <- c("indice", "insample")
# Remove a variavel lista e assets que agora sao desnecessarios
rm(lista)

# Estatisticas descritivas retornos-----------------------------------------
df.descritivas <- assets.tbl %>% 
  transmute(id_name = id_name,
            media = map_dbl(ts, ~ mean(.x)),
            mediana = map_dbl(ts, ~median(.x)),
            maximo = map_dbl(ts, ~max(.x)),
            minimo = map_dbl(ts, ~min(.x)),
            desvp = map_dbl(ts, ~sd(.x)),
            assim = map_dbl(ts, ~skewness(.x)),
            curtose = map_dbl(ts, ~kurtosis(.x)),
            jbstat = map_dbl(ts, ~jarqueberaTest(as.timeSeries(.x))@test$statistic),
            jbpvalue = map_dbl(ts, ~jarqueberaTest(as.timeSeries(.x))@test$p.value),
            q10stat = map_dbl(ts, ~Weighted.Box.test(.x, lag = 10, type = "Ljung-Box")$statistic),
            q10pvalue = map_dbl(ts, ~Weighted.Box.test(.x, lag = 10, type = "Ljung-Box")$p.value),
            q2_10stat = map_dbl(ts, ~Weighted.Box.test(.x, lag = 10, type = "Ljung-Box", sqrd.res = TRUE)$statistic),
            q2_10pvalue = map_dbl(ts, ~Weighted.Box.test(.x, lag = 10, type = "Ljung-Box", sqrd.res = TRUE)$p.value),
            nobs = map_int(ts, ~as.integer(length(.x)))) %>% 
  gather(key = stat_name, value = stat_value, -id_name, factor_key = TRUE) %>% 
  spread(key = id_name, value = stat_value)

df.descritivas$stat_name <- c("Média", "Mediana", "Máximo", "Mínimo", "Desvp", "Assimetria", "Curtose exc.",
                              "Jarque-Bera", "", "Q(10)", "", "$Q^2(10)$", "", "N.obs")
colnames(df.descritivas)[1] <- "Descritivas"

# Cria o xtable
ncol_x <- ncol(df.descritivas) + 1
m5 <- matrix(rep(5, ncol_x*7), nrow = 7)
jb <- rep(2, ncol_x)
p <- rep(5, ncol_x)
q <- rep(4, ncol_x)
q2 <- jb
obs <- rep(0, ncol_x)
digits <- rbind(m5, jb, p, q, p, q2, p, obs)
cap <- paste("Estatísticas descritivas dos retornos (amostra completa de",
           format(start+1, "%d/%m/%Y"), "a 30/08/2017).")
tab1 <- xtable(df.descritivas, 
               caption = cap,
               digits = digits,
               label = "tab:descritivas",
               auto = TRUE)
print.xtable(tab1, 
             file = "./tables/artigo-tab-descritivas.tex",
             caption.placement = "top",
             table.placement = "H",
             sanitize.colnames.function = NULL,
             sanitize.text.function = function(x) {x},
             include.rownames = FALSE)

# Graficos retornos------------------------------------------------------

pdf(file = "./figs/artigo-retornos.pdf",
    width = 7,
    height = 9.5,
    colormodel = "grey")
op <- par(mfrow = c(3, 2))
for (i in seq_len(dim(assets.tbl)[1])) {
  print(
    plot(assets.tbl$ts[[i]], 
       #ylim = c(-0.1, 0.1),
       lwd = 1,
       grid.ticks.on = "years",
       major.ticks = "years",
       main = assets.tbl$indice[i],
       yaxis.right = FALSE,
       format.labels = "%Y")
  )
}
par(op)
dev.off() # fecha o arquivo pdf

## Teste de estacionariedade das series de retornos
# Teste Dickey-Fuller encontrado no pacote CADFtest
# Delta_y = mu + theta*t + delta*y_t-1
# adf <- lapply(assets.tbl$ts, CADFtest)
# names(adf) <- assets.tbl$indice
# adf <- enframe(adf)
# adf <- adf %>%
#   mutate(teste = map(value, ~tidy(.x)))
# adf$value <- NULL
# adf <- unnest(adf)
# adf

# Teste para graficos QQ normal
jpeg(filename = "./figs/artigo-qqplots.jpeg",
     width = 640, height = 800, quality = 100)
op <- par(mfrow = c(3,2),
          mar = c(4, 3, 3, 2))
for(i in 1:dim(assets.tbl)[1]){
  qqnormPlot(assets.tbl$ts[[i]], labels = FALSE, title = FALSE, mtext = FALSE,
             main = assets.tbl$id_name[[i]], 
             xlab = "")
             #ylab = "Amostra")
}
par(op)
dev.off()

# Modelo eGARCH in Sample-----------------------------------------------------------
## Ja as distribuicoes de zt a normal e t-Student nao apresentam bom fit
## A Johnson, GED, NIG, SkewStudent e a Ghyp sao melhores
## Lembrando, o modelo das perdas é AR(1) e a volatilidade é eGARCH(2,1)
## L_t=mu_t+e_t      mu_t=mu+ar1*mu_t-1+e_t
## e_t=sigma_t*z_t   ln(sigma^2_t)=omega+alpha1*z_t-1+gamma1(|z_t-1|-E[|zt-1|])+beta1*ln(sigma^2_t-1)
## LEMBRAR: ts contem os retornos e não as perdas!

# Uma especificacao para cada ativo
ruspec <- ugarchspec(mean.model = list(armaOrder = c(1,0)),
                     variance.model = list(model = "eGARCH", garchOrder = c(2,1)),
                     distribution.model = "norm")
garch.specs <- replicate(length(assets), ruspec)
names(garch.specs) <- assets
garch.specs <- enframe(garch.specs)
colnames(garch.specs) <- c("indice", "spec")
## Modelando as PERDAS!! parametro eh loss para a funcao ugarchfit
# Deixamos um numero outsample para fazer o backtesting. Diferente para cada ativo
# garch.models vai conter o modelo eGARCH dentro da amostra. Apresentar os 
# parametros e seus erros padrao robustos
# Depois apresentar novamente estatisticas como JB, Q e Q^2 para os erros padronizados
garch.models <- assets.tbl[,1:3] %>% 
  inner_join(insample, by = "indice") %>% 
  inner_join(garch.specs, by = "indice") %>% 
  mutate(loss = map(ts, ~-.x),
         loss_in = map2(ts, insample, ~-.x[1:.y]),
         garch_fit = map2(spec, loss_in, ~ugarchfit(.x, 
                                                    .y,
                                                    solver = "hybrid")),
         ts = NULL)

# Mostra a convergencia para cada modelo
#lapply(garch.models$garch_fit, convergence)

# Sumarios dos modelos GARCH
#show(garch.models$garch_fit[[1]])

## Construindo a tabela com os parametros estimados do eGARCH in sample
garch.models.par <- garch.models %>% 
  transmute(id_name = id_name,
            mu = map(.$garch_fit, ~.x@fit$robust.matcoef["mu", c(1, 4)]), # Estimativa e P-valor
            ar1 = map(.$garch_fit, ~.x@fit$robust.matcoef["ar1", c(1, 4)]),
            omega = map(.$garch_fit, ~.x@fit$robust.matcoef["omega", c(1, 4)]),
            alpha1 = map(.$garch_fit, ~.x@fit$robust.matcoef["alpha1", c(1, 4)]),
            alpha2 = map(.$garch_fit, ~.x@fit$robust.matcoef["alpha2", c(1, 4)]),
            beta1 = map(.$garch_fit, ~.x@fit$robust.matcoef["beta1", c(1, 4)]),
            gamma1 = map(.$garch_fit, ~.x@fit$robust.matcoef["gamma1", c(1, 4)]),
            gamma2 = map(.$garch_fit, ~.x@fit$robust.matcoef["gamma2", c(1, 4)])) %>% 
  gather(key = stat_name, value = stat_value, -id_name, factor_key = TRUE) %>% 
  spread(key = id_name, value = stat_value) %>% 
  unnest() 

  #colnames(garch.models.par) <- garch.models$id_name
garch.models.par$stat_name <- c("$\\mu$", "",
                                "$\\phi_1$", "",
                                "$\\omega$", "",
                                "$\\alpha_1$", "",
                                "$\\alpha_2$", "",
                                "$\\beta_1$", "",
                                "$\\gamma_1$", "",
                                "$\\gamma_2$", "")

colnames(garch.models.par)[1] <- "Parâmetros"
# Xtable
cap <- paste("Par\\^ametros estimados do modelo eGARCH. Valores p apresentados de acordo 
com erros padrão robustos e valores menores que 0,01 não são mostrados. (Período 
             dentro da amostra entre",
             format(start+1, "%d/%m/%Y"), "a",
             format(end, "%d/%m/%Y"), ").")
tab2 <- xtable(garch.models.par, 
               caption = cap,
               digits = 5,
               label = "tab:garchcoef",
               auto = TRUE)
print.xtable(tab2, 
             file = "./tables/artigo-tab-garchcoef.tex",
             caption.placement = "top",
             table.placement = "H",
             sanitize.colnames.function = NULL,
             sanitize.text.function = function(x) {x},
             include.rownames = FALSE)

# Gerar 6 figuras com estes 4 graficos ACF
for(i in 1:dim(garch.models)[1]) {
  pdf(file = paste0("./figs/artigo-acf-", garch.models$id_name[i], ".pdf"),
       width = 7, height = 7, colormodel = "grey")
  op <- par(mfrow=c(2,2))
  plot(garch.models$garch_fit[[i]], which = 4)
  plot(garch.models$garch_fit[[i]], which = 5)
  plot(garch.models$garch_fit[[i]], which = 10)
  plot(garch.models$garch_fit[[i]], which = 11)
  par(op)
  dev.off()
}
file.rename(c("./figs/artigo-acf-S&P500.pdf", "./figs/artigo-acf-S&P TSE.pdf"), 
            c("./figs/artigo-acf-SP500.pdf", "./figs/artigo-acf-SP-TSE.pdf"))
## Estatisticas modelo eGARCH in sample
# JB, Q, Q^2 e Sign bias para os residuos padronizados
garch.models.stats <- garch.models %>%
  transmute(id_name = id_name,
            resid_z = map(.$garch_fit, ~residuals(.x, standardize = TRUE)),
            curtose = map_dbl(resid_z, ~kurtosis(.x)),
            jbstat = map_dbl(resid_z, ~jarqueberaTest(as.timeSeries(.x))@test$statistic),
            jbpvalue = map_dbl(resid_z, ~jarqueberaTest(as.timeSeries(.x))@test$p.value),
            q10stat = map_dbl(resid_z, ~Weighted.Box.test(.x, lag = 10, 
                                                      type = "Ljung-Box")$statistic),
            q10pvalue = map_dbl(resid_z, ~Weighted.Box.test(.x, lag = 10, 
                                                          type = "Ljung-Box")$p.value),
            q2_10stat = map_dbl(resid_z, ~Weighted.Box.test(.x, lag = 10, 
                                                            type = "Ljung-Box", 
                                                            sqrd.res = TRUE)$statistic),
            q2_10pvalue = map_dbl(resid_z, ~Weighted.Box.test(.x, lag = 10, 
                                                            type = "Ljung-Box", 
                                                            sqrd.res = TRUE)$p.value),
            sign_bias = map_dbl(.$garch_fit, ~signbias(.x)[1,1]),
            sign_pvalue = map_dbl(.$garch_fit, ~signbias(.x)[1,2])) %>% 
  mutate(resid_z = NULL) %>% 
  gather(key = stat_name, value = stat_value, -id_name, factor_key = TRUE) %>% 
  spread(key = id_name, value = stat_value)

garch.models.stats$stat_name <- c("Curtose exc.",
                                  "Jarque-Bera", "",
                                  "Q(10)", "",
                                  "$Q^2(10)$", "",
                                  "Sign-bias", "")
colnames(garch.models.stats)[1] <- "Estatística"
# Xtable
cap <- paste("Estatísticas de diagnóstico para o modelo eGARCH. 
               (Período dentro da amostra entre",
             format(start+1, "%d/%m/%Y"), "a",
             format(end, "%d/%m/%Y"), ").")

tab3 <- xtable(garch.models.stats, 
               caption = cap,
               digits = 5,
               label = "tab:garchstats",
               auto = TRUE)
print.xtable(tab3, 
             file = "./tables/artigo-tab-garchstats.tex",
             caption.placement = "top",
             table.placement = "H",
             sanitize.colnames.function = NULL,
             sanitize.text.function = function(x) {x},
             include.rownames = FALSE)

# Modelo EVT para os residuos padronizados in Sample --------------------------------

## Os residuos podem ser retirados atraves do metodo residuals com a opcao "standardize=T"
## os valores zq (quantile) e sq (shortfall)
evt.models <- garch.models %>% 
  transmute(indice = indice,
            id_name = id_name,
            loss_in = loss_in,
            insample = insample,
            resid_z = map(garch_fit, ~coredata(residuals(.x, standardize = TRUE))),
            mut = map(garch_fit, ~coredata(fitted(.x))),
            sigmat = map(garch_fit, ~coredata(sigma(.x))),
            gpdfit = map(resid_z, ~gpd(.x, threshold = quantile(.x, u_quant))),
            Nu = map_int(gpdfit, ~.x$n.exceed),
            u = map_dbl(gpdfit, ~.x$threshold),
            xi = map_dbl(gpdfit, ~.x$par.ests[1]),
            beta = map_dbl(gpdfit, ~.x$par.ests[2]),
            xi_se = map_dbl(gpdfit, ~.x$par.ses[1]),
            beta_se = map_dbl(gpdfit, ~.x$par.ses[2]),
            risk = map(gpdfit, ~riskmeasures(.x, c(0.975, 0.990))),
            zq975= map_dbl(risk, ~.x[1, "quantile"]),
            zq990= map_dbl(risk, ~.x[2, "quantile"]),
            sq975= map_dbl(risk, ~.x[1, "sfall"]),
            sq990= map_dbl(risk, ~.x[2, "sfall"]))
# ## Teste com o pacote evir
# testez <- coredata(residuals(garch.models$garch_fit[[1]], standardize = TRUE))
# teste_evir <- gpd(testez, threshold = quantile(testez, 0.925)) # Mesmos valores do fExtremes
# meplot(testez, type = "l")
# shape(testez, models = 10, start = 90, end = 150)
# ## Teste com o pacote evd
# teste_evd <- fpot(testez, quantile(testez, 0.95))
# fitted.values(teste_evd)
# std.errors(teste_evd)  # Iguais ao fExtremes e evir
# mrlplot(testez)
evtcoef <- evt.models %>% 
  transmute(id_name = id_name,
            insample = insample,
            u = u,
            Nu = Nu,
            xi = xi,
            xi_se = xi_se,
            beta = beta,
            beta_se = beta_se,
            zq975 = zq975,
            zq990 = zq990) %>% 
  gather(key = stat_name, value = stat_value, -id_name, factor_key = TRUE) %>% 
  spread(key = id_name, value = stat_value)

evtcoef$stat_name <- c("N.obs.", "Limiar $u$", "Num.exc. $N_u$", "Par. forma $\\xi$", "Erro padrão",
                       "Par. escala $\\psi$", "Erro padrão", "Quantil $z_{0.975}$", "Quantil $z_{0.990}$")
colnames(evtcoef)[1] <- ""
# Xtable
cap <- paste("Parâmetros estimados para o modelo EVT dos resíduos padronizados. 
               (Período dentro da amostra entre",
             format(start+1, "%d/%m/%Y"), "a",
             format(end, "%d/%m/%Y"), ").")

tab4 <- xtable(evtcoef, 
               caption = cap,
               digits = 5,
               label = "tab:evtcoef",
               auto = TRUE)
print.xtable(tab4, 
             file = "./tables/artigo-tab-evtcoef.tex",
             caption.placement = "top",
             table.placement = "H",
             sanitize.colnames.function = NULL,
             sanitize.text.function = function(x) {x},
             include.rownames = FALSE)

## Gráficos para analisar a qualidade do gpdFit
pdf(file = "./figs/artigo-gpdfit.pdf",
    width = 7, height = 8,
    colormodel = "grey")
  op <- par(mfrow=c(3,2))
  for(i in seq_len(dim(evt.models)[1])){
    plot(evt.models$gpdfit[[i]], main = evt.models$id_name[i])
  }
  par(op)
dev.off()

# for(i in seq_len(dim(garch.models)[1])) {
#   jpeg(filename = paste0("./figs/artigo-evtgof-", garch.models$indice[i], ".jpeg"),
#        width = 800, height = 800, quality = 100)
#   op <- par(mfrow=c(2,2))
#   plot(evt.models$gpdfit[[i]], which='all')
#   par(op)
#   dev.off()
# }

# Reconstruindo o VaR e o ES condicionais in Sample ---------------------------------

## VaR: xq_t = mu_t+1 + sigma_t+1*zq
## ES: Sq_t = mu_t+1 + sigma_t+1*sq
# riskmeasures <- evt.models %>%
#   transmute(indice = indice,
#             id_name = id_name,
#             loss_in = loss_in,
#             VaR975 = pmap(., ~(..6+..7*..15)), ## VaR = mu_t+1 + sigma_t+1*zq
#             VaR990 = pmap(., ~(..6+..7*..16)),
#             ES975 = pmap(., ~(..6+..7*..17)),  ## ES = mu_t+1 + sigma_t+1*sq
#             ES990 = pmap(., ~(..6+..7*..18)))
# %>% 
#   mutate(out_VaR975 = pmap(., ~..6[c((..5+1):(..5+..4-1))]), #out_VaR = VaR[(n_old+1):(n_old+out-1)]
#          out_VaR990 = pmap(., ~..7[c((..5+1):(..5+..4-1))]),
#          out_ES975 = pmap(., ~..8[c((..5+1):(..5+..4-1))]),  #out_ES = ES[(n_old+1):(n_old+out-1)]
#          out_ES990 = pmap(., ~..9[c((..5+1):(..5+..4-1))]),
#          out_loss = map2(loss, n_old, ~coredata(.x[-c(1:(.y+1))])[, 1, drop = TRUE])) #out_los = loss[-c(1:n_old+1)]
# Plotando os valores dentro da amostra
# plot_risks <- function(loss, VaR, ES, id_name) {
#   xts <- merge(loss = loss, VaR = VaR, ES = ES)
#   plot <- ggplot(xts, aes(x = Index))+
#     geom_line(aes(y = loss), color = "black")+
#     geom_line(aes(y = VaR), color = "red")+
#     geom_line(aes(y = ES), color = "darkgreen")+
#     labs(x = "", y = "", title = id_name)
#   return(plot)
# }
# VaR_plots <- riskmeasures %>% 
#   transmute(VaR975_plot = pmap(., ~plot_risks(..3, ..4, ..6, ..2)),
#             VaR990_plot = pmap(., ~plot_risks(..3, ..5, ..7, ..2)))
# VaR_plots$VaR975_plot[[1]]+
#   coord_cartesian(xlim = c(as.Date("2011-08-31"), 
#                            as.Date("2014-08-31")))
# # Verifica quantas violacoes
# sum(riskmeasures$loss_in[[1]] > riskmeasures$VaR975[[1]])
# 
# grid.arrange(grobs = VaR_plots$VaR975_plot)
# grid.arrange(grobs = VaR_plots$VaR990_plot)

# Backtesting com refit ----------------------------------------------
# Monta o tibble para as estimacoes out of sample
assets_os.tbl <- garch.models %>% 
  transmute(indice = indice,
            id_name = id_name,
            loss = loss,
            window.size = map_int(loss, ~as.integer(ndays(.x[paste0("/", end)]))),
            n.roll = map2_int(loss, window.size, ~length(.x)-.y),
            spec = spec)
# Periodo fora da amostra + 1, para fazer as comparacoes entre
# VaRt e realizado t+1
# realized tera n.roll-1 observacoes devido ao deslocamento das observacoes de VaR e ES
# para serem comparados.
# As series de VaR e ES no tibble de risco tambem terao n.roll-1 observacoes e estarao nas 
# datas corretas para COMPARACAO (a medida foi calculada no dia anterior)
realized <- assets_os.tbl %>% 
  transmute(indice = indice,
            real = pmap(., ~..3[(..4+2):(..4+..5)])) # real = loss[(window.size+2):(window.size+n.roll)]
saveRDS(realized, 
        file = paste0("./output/", format(Sys.Date(), "%Y-%m-%d"), "realized.rds"))
# Testes estatisticos para o VaR ------------------------------------------
# Quais testes fazer?
# VaRTest possui 2 testes, incondicional de Kupiec1995 e condicional de Christoffersen2001
# que nao esta convergindo
# Existe tambem o teste de Christoffersen2004 de tempo entre as violacoes
# Pode ser feito um teste do tipo Model Confidence Set - MCS - com a funcao VaRloss e 
# mcsTest

######### TESTE PARA OS 6 INDICES ###
# roll_fit(data, window.size, n.roll, spec, models)
# de assets_os.tbl temos os seguintes numeros
# 1: indice 
# 2: id_name
# 3: loss 
# 4: window.size 
# 5: n.roll 
# 6: spec 
models <- c("cevt", "cnorm", "ct", "uevt", "unorm", "ut", "riskmetrics")
# teste_assets_os <- readRDS("./input/teste_assets_os.rds") # Copia dados de teste
# teste_realized <- readRDS("./input/teste_realized.rds")
# teste_os_roll <- teste_assets_os %>% 
#   transmute(indice = indice,
#             id_name = id_name,
#             roll.fit = pmap(., ~roll_fit(..3, ..4, ..5, ..6, models)))

## ATENCAO! Aqui eh a rotina de estimacao do backtesting!! Rodar somente na AWS!!

# f <- file("log.txt", open = "wt")
# sink(f) # Inicia o log no arquivo
# sink(f, type = "message") # Inclusive mensagens de erro e avisos
# cat("\nInicio do map roll_fit:", as.character(Sys.time()))
# os_roll.tbl <- assets_os.tbl %>%
#   transmute(indice = indice,
#             id_name = id_name,
#             roll.fit = pmap(., ~roll_fit(..3, ..4, ..5, ..6, models)))
# cat("\nFim do map roll_fit:", as.character(Sys.time()))
# saveRDS(os_roll.tbl,
#         file = paste0("./output/", format(Sys.Date(), "%Y-%m-%d"), "os_roll_tbl.rds"))
# sink(type = "message")
# sink() # Finaliza o log

os_roll.tbl <- readRDS(file = "./output/2018-01-02os_roll_tbl.rds")
realized <- readRDS(file = "./output/2018-01-01realized.rds")

os_roll_unnest <- os_roll.tbl %>% 
  unnest() %>% 
  mutate(risk.tbl = map(roll, ~.x$risk.tbl),
         param.xts = map(roll, ~.x$param.xts))

# Extrai a evolucao dos parametros de cada modelo no tempo
param.tbl <- os_roll_unnest %>% 
  transmute(indice = indice,
            id_name = id_name,
            model_type = model_type,
            param.xts = param.xts)
# Extrai a evolucao das medidas de risco de cada modelo, para cada cobertura, no tempo
os_risk.tbl <- os_roll_unnest %>% 
  transmute(indice = indice,
            id_name = id_name,
            model_type = model_type,
            risk.tbl = risk.tbl) %>% 
  unnest()
format(object.size(os_risk.tbl), units = "Kb") # Verifica o tamanho do objeto

# Grafico VaR e violacoes -------------------------------------------------

# Plota um grafico da evolucao do VaR e das perdas realizadas
cevt99 <- subset(os_risk.tbl, 
              subset = (indice == "GSPC" & model_type == "cevt" & coverage == 0.01),
              select = VaR.xts)$VaR.xts[[1]]
uevt99 <- subset(os_risk.tbl, 
                  subset = (indice == "GSPC" & model_type == "uevt" & coverage == 0.01),
                  select = VaR.xts)$VaR.xts[[1]]
real <- subset(realized,
               subset = indice == "GSPC",
               select = real)$real[[1]]

plot(real, 
     ylim = c(-0.1, 0.1),
     lwd = 1,
     grid.ticks.on = "years",
     main = "S&P500 EVT condicional vs incondicional")
lines(cevt99, 
      col = "red",
      lty = "dotted")
# Abre o arquivo PDF
pdf(file = "./figs/artigo-sp500evt.pdf",
    width = 7,
    height = 7,
    colormodel = "grey")
lines(uevt99,
      col = "darkgreen",
      lty = "dashed")
dev.off() # fecha o arquivo pdf

## Plot do VaR e violacoes
teste_varplot <- os_risk.tbl %>% 
  subset(subset = (indice == "BVSP" & 
                     (model_type == "cevt" | model_type == "unorm") & 
                     coverage == 0.01)) %>% 
  left_join(realized, by = "indice")
# Abre o arquivo PDF
pdf(file = "./figs/artigo-varviol.pdf",
    width = 7,
    height = 7)
op <- par(mfrow = c(2,1))
VaRplot(teste_varplot$coverage, -teste_varplot$real[[1]], -teste_varplot$VaR.xts[[1]])
title(main = "Violações ao VaR. Modelo CEVT99% BVSP")
VaRplot(teste_varplot$coverage, -teste_varplot$real[[2]], -teste_varplot$VaR.xts[[2]])
title(main = "Violações ao VaR. Modelo UNORM99% BVSP")
par(op)
dev.off() # fecha o arquivo pdf

## Tabela com os percentuais de violacoes
varviolations.tbl <- os_risk.tbl %>% 
  left_join(realized, by = "indice") %>% 
  mutate(violations = map2_dbl(VaR.xts, real, ~100*sum(.y > .x)/length(.y)),
         cov = 100*coverage) %>% 
  select(id_name, cov, model_type, violations) %>% 
  spread(key = id_name, value = violations) %>% 
  mutate(model_type = map_chr(model_type, ~switch(.x,
                                              cevt = "EVT Condicional",
                                              cnorm = "Normal Condicional",
                                              ct = "t-Student Condicional",
                                              uevt = "EVT Incond. Filtrada",
                                              unorm = "Normal Incondicional",
                                              ut = "t-Student Incondicional",
                                              riskmetrics = "RiskMetrics") # Fim do switch
  )) # Fim do map_chr e mutate
colnames(varviolations.tbl)[2] <- "Modelo"
varviolations.tbl <- add_row(varviolations.tbl, Modelo = "Cobertura = 1\\%", 
                             .before = 1)
varviolations.tbl <- add_row(varviolations.tbl, Modelo = "Cobertura = 2.5\\%", 
                             .after = ceiling(dim(varviolations.tbl)[1]/2))
varviolations.tbl$cov <- NULL # Retira a coluna cov, que nao eh mais necessaria

# Xtable
cap <- paste("Percentual de violações. (Período fora da amostra entre",
             format(backstart+1, "%d/%m/%Y"), "e 30/08/2017).")

tab5 <- xtable(varviolations.tbl, 
               caption = cap,
               digits = 2,
               label = "tab:varviol",
               auto = TRUE)
print.xtable(tab5, 
             file = "./tables/artigo-tab-varviol.tex",
             caption.placement = "top",
             table.placement = "H",
             sanitize.colnames.function = NULL,
             sanitize.text.function = function(x) {x},
             include.rownames = FALSE)

## Faz os testes de VaR para os 6 indices
## VaRTest(alpha = 0.05, actual, VaR, conf.level = 0.95)
## var_test(cover = 0.025, loss, var, conf.level = 0.95) Com bootstrap para cc
## vartest(cover = 0.025, loss, var, conf.level = 0.95) engloba uc e DurTest
## LRdur = 2*(uLL-rLL)
vartest.tbl <- os_risk.tbl %>% 
  left_join(realized, by = "indice") %>% 
  transmute(indice = indice,
            id_name = id_name,
            model_type = model_type,
            coverage = coverage,
            VaRtest = pmap(., ~vartest(..4, -coredata(..7), -coredata(..5)))) %>% 
  unnest() %>% 
  mutate(dur.LR = 2*(uLL - rLL)) %>% 
  select(id_name, coverage, model_type, uc.LRstat, uc.LRp, dur.LR, LRp) %>% 
  mutate(uc.LRp = ifelse(uc.LRp > 0.05, NA_real_, uc.LRp),
         LRp = ifelse(LRp > 0.05, NA_real_, LRp)) %>% 
  gather(key = stat_name, value = stat_value, -c(id_name, coverage, model_type), factor_key = TRUE) %>% 
  spread(key = id_name, value = stat_value)
levels(vartest.tbl$stat_name) <- c("LRuc", "LRuc p-valor", "LRdur", "LRdur p-valor")
# Tabela sumarizando os p-valores significativos a 5%
vartest_suma <- vartest.tbl %>% 
  dplyr::filter(str_detect(stat_name, "p-valor")) %>% 
  gather(key = indice, value = p_valor, -c(coverage, model_type, stat_name)) %>% 
  group_by(coverage, model_type, stat_name) %>% 
  summarise(n = sum(p_valor <= 0.05, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(cov_test = paste0(as.character(coverage), stat_name)) %>% 
  select(-c(coverage, stat_name)) %>% 
  spread(key = cov_test, value = n)

colnames(vartest_suma) <- c("Modelo", "LRdur", "LRuc", "LRdur", "LRuc")
colnames(vartest.tbl)[2:3] <- c("Modelo", "Estatística")
  
# Monta as tabelas para o Xtable
vartest.tbl <- add_row(vartest.tbl, Modelo = "Cobertura 1\\%", 
                       .before = 1)
vartest.tbl <- add_row(vartest.tbl, Modelo = "Cobertura 2.5\\%",
                       .after = ceiling(dim(vartest.tbl)[1]/2))

# Xtable
cap <- paste("Testes estatísticos para o VaR. Teste incondicional de Kupiec, \\emph{LRuc}, e teste de
             independência por duração de Christoffersen e Pelletier, \\emph{LRdur}. Os modelos testados
são: EVT condicional (cevt), Normal condicional (cnorm), t-Student condicional (ct), Riskmetrics 
(riskmetrics), EVT incondicioanl (uevt), Normal incondicional (unorm) e t-Student incondicional (ut).
Valores p maiores que 0,05 foram omitidos. (Período fora da amostra entre", format(backstart+1, "%d/%m/%Y"), "e 30/08/2017).")

tab6 <- xtable(vartest.tbl[,-1], 
               caption = cap,
               digits = 2,
               label = "tab:vartest",
               auto = TRUE)
print.xtable(tab6, 
             file = "./tables/artigo-tab-vartest.tex",
             caption.placement = "top",
             table.placement = "H",
             sanitize.colnames.function = NULL,
             sanitize.text.function = function(x) {x},
             include.rownames = FALSE,
             tabular.environment="longtable")

# Xtable vartest_suma
cap <- paste("Sumário para o número de rejeições das hipóteses nulas de um modelo 
corretamente especificado. Nível de confiança a 95\\%. De seis índices com 
dois testes, resulta em um total de doze rejeições possíveis. 
(Período fora da amostra entre", format(backstart+1, "%d/%m/%Y"), "e 30/08/2017).")

suma_tex <- knitr::kable(vartest_suma, 
                         format = "latex", 
                         booktabs = TRUE, 
                         caption = cap, 
                         align = "c") %>% 
  add_header_above(c("", "Cobertura 1%" = 2, "Cobertura 2.5%" = 2)) %>% 
  kable_styling(latex_options = "HOLD_position")
write(suma_tex, "./tables/artigo-tab-vartest_suma.tex")

# tab7 <- xtable(vartest_suma, 
#                caption = cap,
#                digits = 2,
#                label = "tab:vartest_suma",
#                auto = TRUE)
# print.xtable(tab7, 
#              file = "./tables/artigo-tab-vartest_suma.tex",
#              caption.placement = "top",
#              table.placement = "H",
#              sanitize.colnames.function = function(x) {x},
#              sanitize.text.function = function(x) {x},
#              include.rownames = FALSE)
###### Teste da funcao var_test
# teste_risk <- os_risk.tbl %>% 
#   subset(subset = (indice == "IPSA" & model_type == "cevt" & coverage == 0.01)) %>% 
#   left_join(realized, by = "indice")
# var_test.tbl <- var_test(cover = teste_risk$coverage,
#                          loss = unlist(teste_risk$real), 
#                          var = unlist(teste_risk$VaR.xts))

# Plot do VaR e violacoes
teste_varplot <- os_risk.tbl %>% 
  subset(subset = (indice == "IPSA" & model_type == "unorm" & coverage == 0.01)) %>% 
  left_join(realized, by = "indice")
VaRplot(teste_varplot$coverage, -teste_varplot$real[[1]], -teste_varplot$VaR.xts[[1]])

## MCS para os VaR atraves da funcao VaRLoss
# significancia a 5% 
# a sequencia dos modelos eh igual a variavel "models"
mcs_test <- os_risk.tbl %>% 
  left_join(realized, by = "indice") %>% 
  transmute(indice = indice,
            id_name = id_name,
            model_type = model_type,
            coverage = coverage,
            VaRloss = pmap(., ~VaRloss(..4, -coredata(..7), -coredata(..5)))) %>% # Troca o sinal!!
  group_by(indice, id_name, coverage) %>% 
  summarise(loss_matrix = list(do.call(cbind, VaRloss))) %>% 
  mutate(mcs_test = map(loss_matrix, ~mcsTest(.x, alpha = 0.05, nboot = 1000))) %>% 
  ungroup()

models_look <- tibble(num = as.numeric(1:7), models = models)
mcs.tbl <- mcs_test %>%
  transmute(id_name = id_name,
            coverage = coverage,
            mcs_test = mcs_test) %>% 
  mutate(pvals = map(mcs_test, ~.x$pvalsR),
         num = map(mcs_test, ~c(.x$excludedR, .x$includedR)),
         mcs_test = NULL) %>% 
  unnest() %>% 
  left_join(models_look, by = "num") %>% 
  select(coverage, models, id_name, pvals) %>% 
  spread(key = id_name, value = pvals)
colnames(mcs.tbl)[2] <- "Modelo"

# Sumario de exclusoes
mcs_suma <- mcs.tbl %>%
  gather(key = id_name, value = pvals, -c(coverage, Modelo)) %>% 
  group_by(Modelo, coverage) %>% 
  summarise(excl = sum(pvals < 0.05)) %>% 
  spread(key = coverage, value = excl)
colnames(mcs_suma) <- c("Modelo", "Cobertura 1\\%", "Cobertura 2.5\\%")

# Monta as tabelas para o Xtable
mcs.tbl <- add_row(mcs.tbl, Modelo = "Cobertura 1\\%", 
                       .before = 1)
mcs.tbl <- add_row(mcs.tbl, Modelo = "Cobertura 2.5\\%",
                       .after = ceiling(dim(mcs.tbl)[1]/2))

# Xtable mcs.tbl
cap <- paste("Teste MCS de Hansen et al. Apresentados os valores p para cada um dos modelos
             dado um nível de cobertura do VaR e índices. Um valor p abaixo do nível de significância
             exclui o modelo do conjunto superior. (Período fora da amostra entre", 
             format(backstart+1, "%d/%m/%Y"), "e 30/08/2017).")

tab8 <- xtable(mcs.tbl[,-1], 
               caption = cap,
               digits = 2,
               label = "tab:mcs",
               auto = TRUE)
print.xtable(tab8, 
             file = "./tables/artigo-tab-mcs.tex",
             caption.placement = "top",
             table.placement = "H",
             sanitize.colnames.function = NULL,
             sanitize.text.function = function(x) {x},
             include.rownames = FALSE)

# Xtable mcs_suma
cap <- paste("Sumário com as exclusões do conjunto de modelos superiores no teste MCS.
             Os modelos condicionais se revelam aqueles com o menor número de exclusões. 
             (Período fora da amostra entre", format(backstart+1, "%d/%m/%Y"), "e 30/08/2017).")

tab9 <- xtable(mcs_suma, 
               caption = cap,
               digits = 2,
               label = "tab:mcs_suma",
               auto = TRUE)
print.xtable(tab9, 
             file = "./tables/artigo-tab-mcs_suma.tex",
             caption.placement = "top",
             table.placement = "H",
             sanitize.colnames.function = function(x) {x},
             sanitize.text.function = function(x) {x},
             include.rownames = FALSE)

# Testes estatisticos para o ES -------------------------------------------
# ATENCAO!! Este teste so faz sentido ser realizados apos os teste de VaR
# aprovarem o modelo
# ESTest(alpha = 0.05, actual, ES, VaR, conf.level = 0.95, boot = FALSE, n.boot = 1000)
# Pelo codigo fonte nao parece ser bem o teste implementado pelo McNeil2000 com o algoritmo
# do Efron1993
ru_estest.tbl <- os_risk.tbl %>% 
  left_join(realized, by = "indice") %>% 
  transmute(indice = indice,
            id_name = id_name,
            model_type = model_type,
            coverage = coverage,
            EStest = pmap(., ~ESTest(..4, -coredata(..7), -coredata(..6), -coredata(..5), # troca o sinal!!
                                     boot = TRUE, n.boot = 1000)))

# ATENCAO!! Este teste so faz sentido ser realizados apos os teste de VaR
# aprovarem o modelo
# Chamando a funcao propria es_test
es.test <- teste_risk.tbl %>%
  left_join(teste_realized, by = "indice") %>%
  transmute(indice = indice,
            id_name = id_name,
            model_type = model_type,
            coverage = coverage,
            EStest = pmap(., ~es_test(..4, coredata(..7), coredata(..6), coredata(..5)))) # Nao troca o sinal

losses <- coredata(riskmeasures$loss_in[[2]])
VaR <- riskmeasures$VaR975[[2]]
ES <- riskmeasures$ES975[[2]]
t_es.test <- es_test(0.0275, losses, ES, VaR, n.boot = 100)


