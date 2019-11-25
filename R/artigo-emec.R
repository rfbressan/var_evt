#' # Revisão do artigo para a revista Empirical Economics
#' 
#' ## Ativos
#' 
#' O parecerista recomendou utilizar outros ativos ou índices que não somente
#' de ações. Podemos considerar as seguintes classes de ativos e buscar seus 
#' melhores benchmarks
#' 
#' * Stocks 
#'     - MSCI World
#' * Commodities
#'     - MSCI World Commodity Producers Price Return USD
#' * Infrastructure
#'     - MSCI World Infrastructure Price Return USD
#' * Currency pairs
#' * Cryptocurrencies (bitcoin)
#'     - Bitcoin USD (BTC-USD) a partir de 2010 somente
#' * Alternative assets
#'     - 
#' * Long-term interest rates
#'     - Treasury Yield 30 Years (^TYX)
#' * Real Estate
#'     - MSCI US REIT (^RMZ)
#'     - FTSE NAREIT (FNER)
#'     
#' ## Códigos BBG
#' 
#' MSCI   BBG       Nome
#' 990100	MSDUWI	  MSCI World USD
#' 128456	RMZ	   	  MSCI US REIT INDEX 
#' 126174	MXWO0CMP	World Commodity Producers Price Return USD
#' 129867	MXWO0INF	World Infrastructure Price Return USD

knitr::opts_knit$set(root.dir = "~/Documentos/Projetos R/mefca") 

#' ## Períodos in-sample e out-of-sample
#' 
#' Evitar terminar o período in-sample exatamente no meio da crise financeira
#' de 2008/09. Pode ser interessante fazer o estudo com 2 backtests diferentes.
#' Um primeiro com período in-sample terminando em 31/12/2007 e backtest indo 
#' até final de 2010. E outro com período in-sample a partir de 01/01/2011 até 
#' 31/12/2014 e out-of-sample até 31/12/2017. 
#+ echo = TRUE
# Inicio ------------------------------------------------------------------

packages <- c("fExtremes", 
              "evir",
              "CADFtest",
              "rugarch", 
              "PerformanceAnalytics", 
              "xtable", 
              "stargazer",
              "tidyverse", 
              "broom", 
              "purrr", 
              "kableExtra", 
              "gridExtra", 
              "WeightedPortTest")
new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if (length(new.packages)) install.packages(new.packages)

library(evir)
#library(xtable)
library(WeightedPortTest)
#library(CADFtest) # Teste Dickey-Fuller
library(PerformanceAnalytics)
library(fExtremes)
library(rugarch)
library(tidyverse)
#library(broom)
#library(gridExtra)
library(kableExtra)
library(stargazer)
source("./R/artigo_fun.R") # Carrega a funcao roll_fit para fazer o backtest

# Não utilizar notação científica
options(scipen = 999)
# Options for kableExtra
options(knitr.table.format = "latex")
options(knitr.kable.NA = "")
# Paths for tables and figures
t_path <- "./tables2/"
fig_path <- "./figs2/"

#' ## Amostra
#' AMOSTRA COM DADOS A PARTIR DE 01-01-2003
#' BTC a partir de 16-07-2010
#' 
#'  Dois periodos para testes
#'  1. High Volatility
#'      - In-Sample: 01-01-2003 a 31-12-2006
#'      - Out-Sample: 01-01-2007 a 31-12-2009
#'  2. Tranquil times
#'      - In-Sample: 01-01-2011 a 31-12-2014
#'      - Out-Sample: 01-01-2015 a 31-12-2017
#+ echo = TRUE    
start1 <- as.Date("2003-01-01")
end1 <- as.Date("2007-06-30")
backstart1 <- end1 + 1
backend1 <- as.Date("2010-12-31")
start2 <- as.Date("2011-01-01")
end2 <- as.Date("2015-06-30")
backstart2 <- end2 + 1
backend2 <- as.Date("2018-12-31")
indate1 <- paste0(start1, "/", end1)
indate2 <- paste0(start2, "/", end2)
outdate1 <- paste0(start1, "/", backend1) # Out-of-Sample TS needs whole period
outdate2 <- paste0(start2, "/", backend2) # Out-of-Sample TS needs whole period

u_quant <- 0.90 # quantile for treshold u
options(xtable.booktabs = TRUE) # utilizar o pacote booktabs no Latex
options(xtable.tabular.environment = "longtable") # Utiliza o pacote longtable no Latex
options(xtable.floating = FALSE) # Nao pode ser utilizado em conjunto com longtable
# Xtable Stargazer notas nas tabelas
xnota <- paste0("\\hline \n \\multicolumn{4}{l}",
                "{\\scriptsize{Note: * p < 0.1; ** p < 0.05; *** p < 0.001}} \n")
note <- "Note: * p < 0.1; ** p < 0.05; *** p < 0.001"
knote <- "* p < 0.1; ** p < 0.05; *** p < 0.001"

# Function to compute log-losses from prices
losses_fun <- function(asset, start) {
  tb <- read_csv(paste0("./input/emec-", asset, ".csv"), 
                 col_types = cols_only(Date = col_date(),
                                       `Adj Close` = col_double()))
  prices <- xts(tb$`Adj Close`, order.by = tb$Date)[paste0(start, "/")]
  colnames(prices) <- "close"
  return(-na.omit(Return.calculate(prices, method = "log")))
}

#' ## Tabela com ativos
#' Gera um tible com uma coluna com o codigo do ativo, a serie de PERDAS - ts e
#' o nome do indice - id_name
#+ echo = TRUE
# Ativos --------------------------------------------------------------

#assets <- c("BTC", "FNER", "MSDUWI", "MXWO0CMP", "MXWO0INF", "TYX")
assets <- c("EURUSD", "FNER", "MSDUWI", "MXWO0CMP", "MXWO0INF", "TYX")
assets_names <- c("Euro/Dollar", "NAREIT", "MSCI-W", "MSCI-C", "MSCI-I", "Treasuries")
lista <- lapply(assets, losses_fun, start1)
names(lista) <- assets
assets.tbl <- enframe(lista) %>% 
  rename(indice = name,
         ts = value) %>% 
  bind_cols(tibble(id_name = as.factor(assets_names))) %>% 
  crossing(period = c("period1", "period2")) %>% 
  mutate(indate = if_else(period == "period1", indate1, indate2),
         outdate = if_else(period == "period1", outdate1, outdate2),
         #' Time series of insample and out-of-sample respectively
         ints = map2(ts, indate, ~.x[.y]),
         outs = map2(ts, outdate, ~.x[.y]),
         #' Size of insample windows?
         insize = map_int(ints, ~as.integer(ndays(.x)))) %>% 
  arrange(period)

# Remove a variavel lista e assets que agora sao desnecessarios
rm(lista)

#' ## Modelos GARCH in-sample
#+ echo = TRUE
# Modelos GARCH in Sample-----------------------------------------------------------
# Ja as distribuicoes de zt a normal e t-Student nao apresentam bom fit
# A Johnson, GED, NIG, SkewStudent e a Ghyp sao melhores
# Lembrando, o modelo das perdas é AR(1) e a volatilidade é eGARCH(1,1)
# L_t=mu_t+e_t      mu_t=mu+ar1*mu_t-1+e_t
# e_t=sigma_t*z_t   ln(sigma^2_t)=omega+alpha1*z_t-1+gamma1(|z_t-1|-E[|zt-1|])+beta1*ln(sigma^2_t-1)
# LEMBRAR: ts contem as perdas!

# Multiple specifications
g_models <- c("sGARCH", "eGARCH", "gjrGARCH")
# Garch orders
orders <- tibble(order = list(c(1,1), c(1, 2), c(2, 1), c(2, 2)))
# Combination of all assets but BTC !!
# Couldn't find an ARMA specfication that eliminates serial autocorrelation
# on BTC
# **Variance targeting** is needed to insure unconditional variance is finite
# and equivalent to sample variance of mean equation residuals. See Introduction
# to rugarch package, section 2.2.1 The standard GARCH model (’sGARCH’)
comb <- tidyr::crossing(indice = assets, g_models, orders) %>% 
  mutate(spec = map2(g_models, order,
                     ~ugarchspec(mean.model = list(armaOrder = c(1,0)),
                                 variance.model = list(model = .x,
                                                       garchOrder = .y,
                                                       variance.targeting = TRUE),
                                 distribution.model = "norm"))) %>% 
  select(indice, spec)

# GARCH specs
garch.specs <- assets.tbl %>% 
  inner_join(comb, by = "indice") %>% 
  select(-ts) 

# GARCH multiple fits. Look at convergence and best BIC
# Modelando as PERDAS!! parametros (in|ou)ts são perdas!
# We need fixed.se = 1 to compute omega's standard error, since our specification
# is now variance targeting
# 
#             jbstat = map(ints, ~unlist(
# jarqueberaTest(as.timeSeries(.x))@test[c("statistic", "p.value")])),

garch.models <- garch.specs %>% 
  # filter(indice == "FNER") %>%
  select(indice, id_name, period, spec, ints, outs, insize) %>% 
  mutate(fit = map2(spec, ints, 
                    ~ugarchfit(.x, .y, solver = "hybrid")),
         convergence = map_dbl(fit, ~convergence(.x)),
         bic = map_dbl(fit, ~infocriteria(.x)[2]),
         model = garch_name(spec)) %>% 
  group_by(id_name, period) %>% 
  filter(bic == min(bic)) %>% 
  mutate(res = map(fit, ~rugarch::residuals(.x, standardize = TRUE)),
         q10 = map(res, ~unlist(
           Weighted.Box.test(as.timeSeries(.x), lag = 10, 
                             type = "Ljung-Box")[c("statistic", "p.value")])),
         q2_10 = map(res, ~unlist(
           Weighted.Box.test(as.timeSeries(.x), lag = 10, 
                             type = "Ljung-Box", sqrd.res = TRUE)[c("statistic", "p.value")])),
         q10_star = map_chr(q10, ~format_stars(.x, 2)),
         q2_10_star = map_chr(q2_10, ~format_stars(.x, 2)),
         jbstat = map(res, ~unlist(
           jarqueberaTest(as.timeSeries(.x))@test[c("statistic", "p.value")])),
         jb_star = map_chr(jbstat, ~format_stars(.x, 2)))

# Table informing models chosen and criteria
g_models.tbl <- garch.models %>% 
  select(period, id_name, bic, model, q10_star, q2_10_star, jb_star) %>% 
  mutate(bic = format(bic, digits = 3, nsmall = 2))


#' ## Modelos EVT in-sample
#+ echo = TRUE
# Modelo EVT para os residuos padronizados in Sample --------------------------

## Os residuos podem ser retirados atraves do metodo residuals com a opcao "standardize=T"
## os valores zq (quantile) e sq (shortfall)
evt.models <- garch.models %>% 
  mutate(mut = map(fit, ~coredata(fitted(.x))),
         sigmat = map(fit, ~coredata(sigma(.x))),
         gpdfit = map(res, ~gpd(.x, threshold = quantile(.x, u_quant))),
         Nobs = map_int(gpdfit, ~.x$n),
         Nu = map_int(gpdfit, ~.x$n.exceed),
         u = map_dbl(gpdfit, ~.x$threshold),
         xi = map_dbl(gpdfit, ~.x$par.ests[1]),
         beta = map_dbl(gpdfit, ~.x$par.ests[2]),
         xi_se = map_dbl(gpdfit, ~.x$par.ses[1]),
         beta_se = map_dbl(gpdfit, ~.x$par.ses[2]),
         risk = map(gpdfit, ~riskmeasures(.x, c(0.975, 0.990))),
         zq975 = map_dbl(risk, ~.x[1, "quantile"]),
         zq990 = map_dbl(risk, ~.x[2, "quantile"]),
         sq975 = map_dbl(risk, ~.x[1, "sfall"]),
         sq990 = map_dbl(risk, ~.x[2, "sfall"])) %>% 
  select(period, id_name, mut, sigmat, gpdfit, Nobs, Nu, u, xi, beta, zq975, zq990)

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
evt_names <- c("$N_{obs}$", "$N_u$", "$u$", "$\\xi$", "$\\psi$", "$z_{0.975}$", "$z_{0.990}$")
evt_note <- "The columns meanings are as follows: $N_{obs}$ and $N_u$ are number of observations and excesses respectively, $u$ is the threshold value, $\\\\xi$ and $\\\\psi$ are respectively the shape and scale parameters of a GPD, and finally $z_{0.975}$ and $z_{0.990}$ are quantiles from standardized residuals $\\\\hat Z$."
evtcoef <- evt.models %>% 
  select(-c(mut, sigmat, gpdfit))
colnames(evtcoef) <- c("", "", evt_names)


#' ## Backtesting com refit
#+ echo = TRUE
# Backtesting com refit ----------------------------------------------
# Monta o tibble para as estimacoes out of sample
assets_os.tbl <- garch.models %>% 
  ungroup() %>% 
  mutate(n.roll = map2_int(outs, insize, ~length(.x)-.y)) %>% 
  select(period, id_name, outs, insize, n.roll, spec)
# Periodo fora da amostra para fazer as comparacoes entre
# VaRt e realizado t+1
# realized tera n.roll observacoes todas dentro do periodo out-of-sample
# O VaR e ES devido ao deslocamento das previsoes, se iniciam no ultimo dia
# in-sample e terminam no penultimo dia out-of-sample
# As series de VaR e ES no tibble de risco tambem terao n.roll observacoes e estarao nas 
# datas corretas para COMPARACAO (a medida foi calculada no dia anterior)
realized <- assets_os.tbl %>% 
  # real = outs[(insize+1):(insize+n.roll)]
  mutate(real = pmap(., ~..3[(..4+1):(..4+..5)])) %>% 
  select(period, id_name, real)
saveRDS(realized,
        file = paste0("./output/", format(Sys.Date(), "%Y-%m-%d"), "realized.rds"))

######### TESTE PARA OS 6 INDICES ###
# roll_fit(data, insize, n.roll, spec, models)
# de assets_os.tbl temos os seguintes numeros
# 1: period
# 2: id_name
# 3: outs 
# 4: insize 
# 5: n.roll 
# 6: spec 
models <- c("cevt", "cnorm", "ct", "uevt", "unorm", "ut", "riskmetrics")

## ATENCAO! Aqui eh a rotina de estimacao do backtesting!! Rodar somente na AWS!!

f <- file("log.txt", open = "wt")
sink(f) # Inicia o log no arquivo
sink(f, type = "message") # Inclusive mensagens de erro e avisos
cat("\nInicio do map roll_fit:", as.character(Sys.time()))
os_roll.tbl <- assets_os.tbl %>%
  transmute(period = period,
            id_name = id_name,
            roll.fit = pmap(., ~roll_fit(..3, ..4, ..5, ..6, models)))
cat("\nFim do map roll_fit:", as.character(Sys.time()))
saveRDS(os_roll.tbl,
        file = paste0("./output/", format(Sys.Date(), "%Y-%m-%d"), "os_roll_tbl.rds"))
sink(type = "message")
sink() # Finaliza o log

os_roll_unnest <- os_roll.tbl %>%
  unnest() %>%
  mutate(risk.tbl = map(roll, ~.x$risk.tbl),
         param.xts = map(roll, ~.x$param.xts))

# Extrai a evolucao dos parametros de cada modelo no tempo
param.tbl <- os_roll_unnest %>%
  transmute(period = period,
            id_name = id_name,
            model_type = model_type,
            param.xts = param.xts)
saveRDS(param.tbl,
        file = paste0("./output/", format(Sys.Date(), "%Y-%m-%d"), "param_tbl.rds"))

# Extrai a evolucao das medidas de risco de cada modelo, para cada cobertura, no tempo
os_risk.tbl <- os_roll_unnest %>%
  transmute(period = period,
            id_name = id_name,
            model_type = model_type,
            risk.tbl = risk.tbl) %>%
  unnest()
format(object.size(os_risk.tbl), units = "Kb") # Verifica o tamanho do objeto
saveRDS(os_risk.tbl,
        file = paste0("./output/", format(Sys.Date(), "%Y-%m-%d"), "os_risk_tbl.rds"))

## Fim da AWS -----------------------------------------------------

# Load files from last fit -----------------------------------------
# os_roll.tbl <- readRDS(file = "./output/2018-01-02os_roll_tbl.rds")
# realized <- readRDS(file = "./output/2019-02-20realized.rds")
# os_risk.tbl <- readRDS(file = "./output/2019-02-20os_risk_tbl.rds")
# param.tbl <- readRDS(file = "./output/2019-02-15param_tbl.rds")

#' ## Code exclusively for testing purposes
#+ echo = TRUE
# Testing roll_fit -----------------------------------------------------
# teste_assets_os <- assets_os.tbl %>%
#   slice(12)
# teste_realized <- realized %>%
#   slice(12)
# models <- c("cevt")
# teste_os_roll <- teste_assets_os %>%
#   # Fits the 7 models for each asset and period
#   transmute(period = period,
#             id_name = id_name,
#             roll.fit = pmap(., ~roll_fit(..3, ..4, ..5, ..6, models)))
# 
# teste_unnest <- teste_os_roll %>%
#     unnest() %>%
#     mutate(risk.tbl = map(roll, ~.x$risk.tbl),
#            param.xts = map(roll, ~.x$param.xts))
# 
# teste_risk.tbl <- teste_unnest %>%
#   transmute(period = period,
#             id_name = id_name,
#             model_type = model_type,
#             risk.tbl = risk.tbl) %>%
#   unnest()
# teste_param.tbl <- teste_unnest %>%
#   transmute(period = period,
#             id_name = id_name,
#             model_type = model_type,
#             param.xts = param.xts)
# ## Filtra os_risk.tbl para teste 
# os_risk.tbl <- os_risk.tbl %>%
#   filter(period == "period2" & id_name == "Treasuries" & model_type == "cevt")
# plot(teste_realized$real[[1]], col = "darkgreen")
# lines(teste_risk.tbl$VaR.xts[[1]])
# lines(os_risk.tbl$VaR.xts[[1]], col = "red")

#' ## Tabelas do documento
#' ### Descriptive Statistics
#+ echo = TRUE
# Estatisticas descritivas -----------------------------------------
df.descritivas <- assets.tbl %>% 
  transmute(id_name = id_name,
            period = period,
            media = map_chr(ints, ~format(mean(.x), digits = 2, nsmall = 2)),
            mediana = map_chr(ints, ~format(median(.x), digits = 2, nsmall = 2)),
            maximo = map_chr(ints, ~format(max(.x), digits = 2, nsmall = 2)),
            minimo = map_chr(ints, ~format(min(.x), digits = 2, nsmall = 2)),
            desvp = map_chr(ints, ~format(sd(.x), digits = 2, nsmall = 2)),
            assim = map_chr(ints, ~format(skewness(.x), digits = 2, nsmall = 2)),
            curtose = map_chr(ints, ~format(kurtosis(.x), digits = 2, nsmall = 2)),
            jbstat = map(ints, ~unlist(
              jarqueberaTest(as.timeSeries(.x))@test[c("statistic", "p.value")])),
            q10stat = map(ints, ~unlist(
              Weighted.Box.test(as.timeSeries(.x), lag = 10, type = "Ljung-Box")[c("statistic", "p.value")])),
            q2_10stat = map(ints, ~unlist(
              Weighted.Box.test(as.timeSeries(.x), lag = 10, type = "Ljung-Box", 
                                sqrd.res = TRUE)[c("statistic", "p.value")])),
            nobs = map_chr(ints, ~as.character(length(.x)))) %>% 
  mutate(jbstat = map_chr(jbstat, ~format_stars(.x, 2)),
         q10stat = map_chr(q10stat, ~format_stars(.x, 2)),
         q2_10stat = map_chr(q2_10stat, ~format_stars(.x, 2))) %>% 
  gather(key = stat_name, value = stat_value, -c(id_name, period), factor_key = TRUE) %>% 
  spread(key = id_name, value = stat_value)

df.descritivas$stat_name <- rep(c("Mean", "Median", "Maximum", "Minimum", 
                              "Std.Dev.", "Skew", "Exc.Kurtosis",
                              "Jarque-Bera", "Q(10)", "$Q^2(10)$", "N.obs"), 2)
colnames(df.descritivas)[2] <- "Description"

# kableExtra for descriptive statistics
cap <- paste("Descriptive statistics for losses during in-sample periods. Period 1 spans from", 
             start1, "to",  end1, "and Period 2 from", start2, "to", end2, ".")
# Sets label for this table
knitr::opts_current$set(label = "descritivas")
desc_tbl <- df.descritivas[, -1] %>% 
  kable(caption = cap,
        # col.names = c("", "BIC", "Model", "Q(10)", "$Q^2(10)$"),
        escape = FALSE,
        booktabs = TRUE) %>% 
  kable_styling(latex_options = "HOLD_position", 
                font_size = 8,
                full_width = FALSE) %>% 
  group_rows("Period 1", 1, 11) %>%
  group_rows("Period 2", 12, 22) %>%
  footnote(general = knote, general_title = "Note: ", footnote_as_chunk = TRUE)
write(desc_tbl, file = paste0(t_path, "emec-tab-descritivas.tex"))

#' ### BTC stationarity
#+ echo = TRUE
# Estacionariedade BTC ----------------------------------------------------

adf_tbl <- assets.tbl %>%
#  filter(indice == "BTC") %>%
  mutate(adf = map(ints, ~table_adf(.x))) %>%
  dplyr::select(indice, adf) %>%
  unnest()
colnames(adf_tbl) <- c("Ativo", "Deriva", "Tendência")
# kableExtra for ADF statistics
cap <- "Teste aumentado de Dickey-Fuller para raiz unitária."
# Sets label for this table
knitr::opts_current$set(label = "dickeyfuller")
desc_tbl <- adf_tbl %>% 
  kable(caption = cap,
        # col.names = c("", "BIC", "Model", "Q(10)", "$Q^2(10)$"),
        escape = FALSE,
        booktabs = TRUE) %>% 
  kable_styling(latex_options = "HOLD_position", 
                font_size = 12,
                full_width = FALSE) %>% 
  group_rows("Period 1", 1, 6) %>%
  group_rows("Period 2", 7, 12) %>%
  footnote(general = knote, general_title = "Nota: ", footnote_as_chunk = TRUE)
write(desc_tbl, file = paste0(t_path, "emec-tab-adf.tex"))
# stargazer(adf_tbl,
#           type = "latex",
#           summary = FALSE,
#           notes = note,
#           header = FALSE,
#           title = "Augmented Dickey-Fuller unit root test.",
#           rownames = FALSE,
#           out = "./tables2/emec-tab-adf.tex",
#           label = "tab:dickeyfuller")




#' ### Garch Models
#+ echo = TRUE
# GARCH selected models ----------------------------------------------
# The best GARCH models have been selected for each asset based on BIC
# Sets label for this table
knitr::opts_current$set(label = "garchmodels")
m_tbl <- kable(g_models.tbl[, -1], 
               caption = "Model selection for each asset and in-sample period.",
               col.names = c("", "BIC", "Model", "Q(10)", "$Q^2(10)$", "J-B"),
               escape = FALSE,
               booktabs = TRUE) %>% 
  kable_styling(latex_options = "HOLD_position", 
                # font_size = 8,
                full_width = FALSE) %>% 
  group_rows("Period 1", 1, 6) %>% 
  group_rows("Period 2", 7, 12) %>% 
  footnote(general = knote, general_title = "Note: ", footnote_as_chunk = TRUE)
write(m_tbl, file = paste0(t_path, "emec-tab-garchmodels.tex"))

#' ### Parametros GARCH
#+ echo = TRUE
# Tabela de parametros GARCH ---------------------------------------
# Construindo a tabela com os parametros estimados do GARCH in sample
garch.par.names <- c("$\\mu$", "$\\phi_1$", "$\\omega$", "$\\alpha_1$", 
                     "$\\beta_1$", "$\\gamma_1$", "$\\alpha_2$", 
                     "$\\beta_2$", "$\\gamma_2$")
garch.pars <- garch.models %>% 
  mutate(values = map(fit, ~garch_pars(.x)),
         names = map(values, ~names(.x))) %>% 
  select(period, id_name, names, values) %>% 
  unnest() %>% 
  mutate(names = factor(names, levels = unique(names))) %>% 
  spread(key = id_name, value = values) %>% 
  mutate(names = garch.par.names) %>% # It's grouped period no need to replicate
  # Filters only parameters that is present in at least one asset
  filter_at(3:8, any_vars(. != ""))

colnames(garch.pars)[2] <- "Parameter"
cap <- "Garch parameters estimated for both in-sample periods."
# Last line of period 1
last_p1 <- dplyr::last(which(garch.pars$period == "period1"))
# Last line of period 2
last_p2 <- nrow(garch.pars)

# kableExtra
# Sets label for this table
knitr::opts_current$set(label = "garchpars")
pars_tbl <- garch.pars[, -1] %>% 
  kable(caption = cap,
        escape = FALSE,
        booktabs = TRUE) %>% 
  kable_styling(latex_options = "HOLD_position",
                font_size = 8,
                full_width = FALSE) %>% 
  group_rows("Period 1", 1, last_p1) %>% 
  group_rows("Period 2", last_p1 + 1L, last_p2) %>% 
  footnote(general = knote, general_title = "Note: ", footnote_as_chunk = TRUE)
write(pars_tbl, file = paste0(t_path, "emec-tab-garchpars.tex"))




#' ### Parametros EVT
#+ echo = TRUE
# Tabela EVT -----------------------------------------------------------------
# kableExtra
cap <- "EVT results and Generalized Pareto distribution parameters estimated from both in-sample periods."
# Sets label for this table
knitr::opts_current$set(label = "evtcoef")
evt_tbl <- kable(evtcoef[, -1],
                 caption = cap,
                 digits = 3,
                 escape = FALSE,
                 booktabs = TRUE) %>% 
  kable_styling(latex_options = "HOLD_position",
                # font_size = 10,
                full_width = FALSE) %>% 
  group_rows("Period 1", 1, 6) %>% 
  group_rows("Period 2", 7, 12) %>% 
  footnote(general = evt_note, footnote_as_chunk = TRUE, 
           escape = FALSE, threeparttable = TRUE)
write(evt_tbl, file = paste0(t_path, "emec-tab-evtcoef.tex"))


#' ### Violations percentage points
#+ echo = TRUE
# Tabela com os percentuais de violacoes -----------------------------------
varviol.tbl <- os_risk.tbl %>% 
  left_join(realized, by = c("period", "id_name")) %>% 
  mutate(violations = map2_dbl(VaR.xts, real, ~100*sum(.y > .x)/length(.y)),
         cov = 100*coverage) %>% 
  select(cov, period, id_name, model_type, violations) %>% 
  spread(key = id_name, value = violations) %>% 
  mutate(model_type = map_chr(model_type, ~switch(.x,
                                              cevt = "Cond. EVT",
                                              cnorm = "Cond. Normal",
                                              ct = "Cond. Student-t",
                                              uevt = "Uncond. EVT",
                                              unorm = "Uncond. Normal",
                                              ut = "Uncond. Student-t",
                                              riskmetrics = "RiskMetrics")),
         cov = case_when(cov == 1 ~ "Coverage of 1.0%",
                         cov == 2.5 ~ "Coverage of 2.5%"))
colnames(varviol.tbl)[3] <- "Model"

# kableExtra
cap <- paste("Percentage of violations for each coverage.")
# Sets label for this table
knitr::opts_current$set(label = "varviol")
varviol_tbl <- kable(varviol.tbl,
                 caption = cap,
                 digits = 3,
                 escape = FALSE,
                 booktabs = TRUE) %>% 
  kable_styling(latex_options = "HOLD_position",
                font_size = 8,
                full_width = FALSE) %>% 
  collapse_rows(1:3, row_group_label_position = "stack",
                latex_hline = "major")  
  # footnote(general = evt_note, footnote_as_chunk = TRUE, 
  #          escape = FALSE, threeparttable = TRUE)
write(varviol_tbl, file = paste0(t_path, "emec-tab-varviol.tex"))

#' ### Kupiec and Christoffersen tests
#+ echo = TRUE
# Kupiec and Christoffersen tests -------------------------------------
## Faz os testes de VaR para os 6 indices
## VaRTest(alpha = 0.05, actual, VaR, conf.level = 0.95)
## var_test(cover = 0.025, ret, var, conf.level = 0.95) Com bootstrap para cc
## vartest(cover = 0.025, ret, var, conf.level = 0.95) engloba uc e DurTest
## LRdur = 2*(uLL-rLL)
var_risk.tbl <- os_risk.tbl %>% 
  left_join(realized, by = c("period", "id_name")) %>% 
  transmute(period = period,
            id_name = id_name,
            model_type = model_type,
            coverage = coverage,
            VaRtest = pmap(., ~vartest(..4, -coredata(..7), -coredata(..5)))) %>% 
  unnest() %>% 
  mutate(uc_stat = map2_chr(uc.LRstat, uc.LRp, ~format_stars(c(.x, .y), 
                                                             digits = 2)),
         dur_stat = map2_chr(dur.LR, LRp, ~format_stars(c(.x, .y), 
                                                        digits = 2))) 
vartest.tbl <- var_risk.tbl %>% 
  select(coverage, period, id_name, model_type, uc_stat, dur_stat) %>%
  gather(key = stat_name, value = stat_value, 
         -c(period, id_name, coverage, model_type), 
         factor_key = TRUE) %>% 
  spread(key = id_name, value = stat_value) %>% 
  mutate(coverage = case_when(coverage == 0.01 ~ "Coverage of 1.0%",
                              coverage == 0.025 ~ "Coverage of 2.5%"))
# Changes names of stats
levels(vartest.tbl$stat_name) <- c("LRuc", "LRdur")
colnames(vartest.tbl)[3:4] <- c("Model", "Statistic")

# kableExtra
cap <- paste("Statistical tests for VaR. Unconditional test of Kupiec, \\emph{LRuc},
             and duration test of Christoffersen e Pelletier, \\emph{LRdur}. Tested models
             are: conditional EVT (cevt), conditional Normal (cnorm), conditional Student-t 
             (ct), Riskmetrics (riskmetrics), unconditional EVT (uevt), unconditional Normal 
             (unorm) and unconditional Student-t (ut).")
# Sets label for this table
knitr::opts_current$set(label = "vartest")
vartest_tbl <- kable(vartest.tbl,
                     caption = cap,
                     digits = 3,
                     escape = FALSE,
                     longtable = TRUE,
                     booktabs = TRUE) %>% 
  collapse_rows(1:3, row_group_label_position = "stack",
                latex_hline = "major") %>% 
  # longtable cannot be scaled down
  kable_styling(latex_options = c("HOLD_position", "repeat_header"),
                repeat_header_method = "replace",
                font_size = 8,
                full_width = FALSE) %>% 
  footnote(general = knote, footnote_as_chunk = TRUE,
           escape = FALSE, threeparttable = TRUE)
# Do not forget to edit raw tex file to
# \setlength{\tabcolsep}{3pt} after \begingroup command
write(vartest_tbl, file = paste0(t_path, "emec-tab-vartest.tex"))

#' ### Summary Kupiec and Christoffersen tests
#+ echo = TRUE
# Summary of Kupiec and Christoffersen tests -------------------------------------

vartest_suma <- var_risk.tbl %>% 
  select(coverage, period, id_name, model_type, uc.LRp, LRp) %>% 
  gather(key = stat_p, value = p_valor, -c(coverage, period, id_name, model_type)) %>% 
  group_by(coverage, period, model_type, stat_p) %>% 
  summarise(n = sum(p_valor <= 0.05, na.rm = TRUE)) %>% 
  ungroup() %>% 
  mutate(cov_test = paste0(as.character(coverage), stat_p)) %>% 
  select(-c(coverage, stat_p)) %>% 
  spread(key = cov_test, value = n) %>% 
  mutate(total = rowSums(.[, 3:6]))

# colnames(vartest_suma) <- c("period", "Model", "LRdur", "LRuc", "LRdur", "LRuc")
# kableExtra  
cap <- paste("Summary for the number of rejections of a well specified model. 
             Confidence level at 95\\%.")
# Sets label for this table
knitr::opts_current$set(label = "vartest_suma")
vartest_suma_tbl <- kable(vartest_suma[, -1], 
                         caption = cap, 
                         escape = FALSE,
                         col.names = c("Model", rep(c("LRdur", "LRuc"), 2), "Total"),
                         align = "c",  
                         booktabs = TRUE) %>% 
  add_header_above(c("", "Coverage 1.0%" = 2, "Coverage 2.5%" = 2)) %>% 
  group_rows("Period 1", 1, 7) %>% 
  group_rows("Period 2", 8, 14) %>% 
  column_spec(6, bold = TRUE) %>% 
  kable_styling(latex_options = "HOLD_position")
write(vartest_suma_tbl, paste0(t_path, "emec-tab-vartest_suma.tex"))

###### Teste da funcao var_test
# teste_risk <- os_risk.tbl %>% 
#   subset(subset = (indice == "IPSA" & model_type == "cevt" & coverage == 0.01)) %>% 
#   left_join(realized, by = "indice")
# var_test.tbl <- var_test(cover = teste_risk$coverage,
#                          loss = unlist(teste_risk$real), 
#                          var = unlist(teste_risk$VaR.xts))

#' ### MCS test
#+ echo = TRUE
# MCS para os VaR atraves da funcao VaRLoss ---------------------------------
# MCS significancia a 5% 
# a sequencia dos modelos eh igual a variavel "models"
mcs_test <- os_risk.tbl %>% 
  left_join(realized, by = c("period", "id_name")) %>% 
  transmute(period = period,
            id_name = id_name,
            model_type = model_type,
            coverage = coverage,
            VaRloss = pmap(., ~VaRloss(..4, -coredata(..7), -coredata(..5)))) %>% # Troca o sinal!!
  group_by(period, id_name, coverage) %>% 
  summarise(loss_matrix = list(do.call(cbind, VaRloss))) %>% 
  mutate(mcs_test = map(loss_matrix, ~mcsTest(.x, alpha = 0.05, nboot = 1000))) %>% 
  ungroup()

models_look <- tibble(num = as.numeric(1:7), models = models)
# pvals are ordered, lower to higher. Model number is ordered acording to
# excluded (lower pvals) and included (higher pvals)
mcs.tbl <- mcs_test %>%
  mutate(pvals = map(mcs_test, ~.x$pvalsR),
         num = map(mcs_test, ~c(.x$excludedR, .x$includedR))) %>% 
  select(-c(loss_matrix, mcs_test)) %>% 
  unnest() %>% 
  left_join(models_look, by = "num") %>% 
  select(coverage, period, id_name, models, pvals) %>% 
  spread(key = id_name, value = pvals) %>% 
  mutate(coverage = case_when(coverage == 0.01 ~ "Coverage of 1.0%",
                              coverage == 0.025 ~ "Coverage of 2.5%"))
colnames(mcs.tbl)[3] <- "Model"

cap <- paste("Model Confidence Set. For each model, period and coverage the 
p-values from range statistic are presented. A value below a specified significance
level excludes the model from the superior set.")
# Sets label for this table
knitr::opts_current$set(label = "mcs")
mcs_tbl <- kable(mcs.tbl,
                 caption = cap,
                 digits = 3,
                 escape = FALSE,
                 longtable = TRUE,
                 booktabs = TRUE) %>% 
  collapse_rows(1:3, row_group_label_position = "stack",
                latex_hline = "major") %>% 
  kable_styling(latex_options = c("HOLD_position", "repeat_header"),
                repeat_header_method = "replace",
                font_size = 8,
                full_width = FALSE)
write(mcs_tbl, file = paste0(t_path, "emec-tab-mcs.tex"))

# MCS sumario de exclusoes ------------------------------
mcs.suma <- mcs.tbl %>%
  gather(key = id_name, value = pvals, -c(coverage, period, Model)) %>% 
  group_by(coverage, period, Model) %>% 
  summarise(excl = sum(pvals < 0.05)) %>% 
  spread(key = coverage, value = excl) %>% 
  ungroup() %>% 
  mutate(total = rowSums(.[, 3:4]))
# kableExtra  
cap <- paste("Summary of exclusions from superior set.")
# Sets label for this table
knitr::opts_current$set(label = "mcs_suma")
mcs_suma_tbl <- kable(mcs.suma[, -1], 
                          caption = cap, 
                          escape = TRUE,
                          col.names = c("Model", "Coverage 1.0%", 
                                        "Coverage 2.5%", "Total"),
                          align = "c",  
                          booktabs = TRUE) %>% 
  group_rows("Period 1", 1, 7) %>% 
  group_rows("Period 2", 8, 14) %>% 
  column_spec(4, bold = TRUE) %>% 
  kable_styling(latex_options = "HOLD_position")
write(mcs_suma_tbl, paste0(t_path, "emec-tab-mcs_suma.tex"))

#' ## Gráficos do documento
#' ### Losses plot
#+ echo = TRUE
# Graficos perdas------------------------------------------------------
# Creates a data frame compatible to ggplot2 facets
ts <- do.call(xts::merge.xts, 
              dplyr::filter(assets.tbl, period == "period1")$ts)
ts <- ts[paste0(start1, "/", backend2)]
colnames(ts) <- assets_names
ts_df <- ggplot2::fortify(ts) %>% 
  gather(key = symbol, value = value, -Index) %>% 
  na.omit()
shades <- data.frame(xmin = c(backstart1, backstart2), 
                     xmax = c(backend1, backend2))
# Opens a pdf file to hold the plots
pdf(file = paste0(fig_path, "emec-Fig1.pdf"),
    width = 7,
    height = 8,
    colormodel = "grey")
ggplot() + 
  geom_rect(data = shades,
            aes(xmin = xmin, xmax = xmax, 
                ymin = -Inf, ymax = Inf),
            alpha = 0.4) +
  geom_line(data = ts_df, aes(x = Index, y = value)) +
  facet_wrap(~symbol, nrow = 3, scales = "free") +
  labs(title = "",
       x = "",
       y = "loss") +
  scale_y_continuous() +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme_bw()
dev.off() # fecha o arquivo pdf

#' ### QQ plots
#+ echo = TRUE
# QQ plots ------------------------------------
# # Teste para graficos QQ normal
# Opens a pdf file to hold the plots
pdf(file = paste0(fig_path, "emec-qqplot.pdf"),
    width = 7,
    height = 8,
    colormodel = "grey")
op <- par(mfrow = c(3,2),
          mar = c(4, 3, 3, 2))
for (i in 1:dim(assets.tbl)[1]) {
  qqnormPlot(assets.tbl$ts[[i]], labels = FALSE, title = FALSE, mtext = FALSE,
             main = assets.tbl$id_name[[i]],
             xlab = "")
             #ylab = "Amostra")
}
par(op)
dev.off()

#' ### ACF dos modelos GARCH
#+ echo = TRUE
# Gerar 6 figuras com estes 4 graficos ACF --------------------------------

# change name Euro/Dollar to Euro-Dollar in garch.models
garch.models1 <- garch.models %>% 
  ungroup() %>% 
  mutate(id_name = str_replace(id_name, "/", "-"))
for (i in 1:dim(garch.models)[1]) {
  pdf(file = paste0(fig_path, "emec-acf-", garch.models1$id_name[i], ".pdf"),
       width = 7, height = 7, colormodel = "grey")
  op <- par(mfrow = c(2,2))
  plot(garch.models$fit[[i]], which = 4)
  plot(garch.models$fit[[i]], which = 5)
  plot(garch.models$fit[[i]], which = 10)
  plot(garch.models$fit[[i]], which = 11)
  par(op)
  dev.off()
}
# file.rename(c("./figs/artigo-acf-S&P500.pdf", "./figs/artigo-acf-S&P TSE.pdf"), 
#             c("./figs/artigo-acf-SP500.pdf", "./figs/artigo-acf-SP-TSE.pdf"))

#' ### EVT goodness-of-fit plots
#+ echo = TRUE
# Gráficos para analisar a qualidade do gpdFit ------------------------------
# Selecionar primeiro 1 e em seguida 0
# repetir 6 vezes
evt.models1 <- evt.models %>% 
  filter(period == "period1")
pdf(file = paste0(fig_path, "emec-gpdfit1.pdf"),
    width = 7, height = 8,
    colormodel = "grey")
op <- par(mfrow = c(3,2))
for (i in seq_len(dim(evt.models1)[1])) {
  plot(evt.models1$gpdfit[[i]], main = evt.models1$id_name[i])
  # qplot(evt.models1$gpdfit[[i]]$data,
  #       xi = evt.models1$gpdfit[[i]]$par.ests["xi"],
  #       threshold = evt.models1$gpdfit[[i]]$threshold,
  #       main = evt.models1$id_name[i])
}
par(op)
dev.off()

evt.models2 <- evt.models %>% 
  filter(period == "period2")
pdf(file = paste0(fig_path, "emec-gpdfit2.pdf"),
    width = 7, height = 8,
    colormodel = "grey")
op <- par(mfrow = c(3,2))
for (i in seq_len(dim(evt.models2)[1])) {
  plot(evt.models2$gpdfit[[i]], main = evt.models2$id_name[i])
}
par(op)
dev.off()
# Remove evt.models[i]
rm(evt.models1, evt.models2)

#' ### VaR violations
#+ echo = TRUE
# Grafico VaR e violacoes -------------------------------------------------

# Plota um grafico da evolucao do VaR e das perdas realizadas
mtype <- c("cevt", "uevt")
cover <- 0.01
real.tbl <- realized %>% 
  mutate(real = map(real, ~fortify(.x))) %>% 
  unnest(real)
var.tbl <- os_risk.tbl %>% 
  filter(model_type %in% mtype & coverage == cover) %>% 
  select(-c(ES.xts, coverage)) %>% 
  mutate(var = map(VaR.xts, ~fortify(.x))) %>% 
  select(-VaR.xts) %>% 
  unnest(var) %>% 
  spread(key = model_type, value = Zq990) %>% 
  left_join(real.tbl, by = c("period", "id_name", "Index")) #%>% 
# tidyr::gather(key = name, value = varloss, -c(period, id_name, Index), factor_key = TRUE)
# Opens a pdf file to hold the plots
pdf(file = paste0(fig_path, "emec-var.pdf"),
    width = 7,
    height = 8,
    colormodel = "grey")
# Plot of all 6 assets for each period of selected model and coverage
ggplot(var.tbl, aes(x = Index)) +
  geom_line(aes(y = cevt), linetype = "solid") +
  geom_line(aes(y = uevt), linetype = "longdash") +
  geom_segment(aes(y = close, xend = Index, yend = 0), size = 0.5) +
  guides(linetype = FALSE) +
  facet_grid(id_name ~ period, scales = "free") +
  labs(title = "",
       x = "",
       y = "Loss / VaR") +
  scale_y_continuous() +
  scale_x_date(date_breaks = "1 years", date_labels = "%Y") +
  theme_bw()
dev.off() # Closes pdf file

## Plot do VaR e violacoes
asset <- "NAREIT"
viol.tbl <- var.tbl %>% 
  dplyr::filter(id_name %in% asset) %>% 
  mutate(viol1 = if_else(close > cevt, cevt, NA_real_),
         viol2 = if_else(close > uevt, uevt, NA_real_))
# Opens a pdf file to hold the plots
pdf(file = paste0(fig_path, "emec-viol.pdf"),
    width = 7,
    height = 8,
    colormodel = "grey")
# Plot of VaR violations on each period of selected asset, models and coverage
ggplot(viol.tbl, aes(x = Index)) +
  geom_line(aes(y = cevt), linetype = "solid") +
  geom_line(aes(y = uevt), linetype = "longdash") +
  geom_segment(aes(y = close, xend = Index, yend = 0), size = 0.4) +
  geom_point(aes(y = viol1), shape = 2, size = 2) +
  geom_point(aes(y = viol2), shape = 4, size = 2) +
  # guides(linetype = FALSE) +
  facet_wrap(~ period, nrow = 2, scales = "free") +
  labs(title = "",
       x = "",
       y = "Loss / VaR") +
  scale_y_continuous() +
  scale_x_date(date_breaks = "1 years", date_labels = "%Y") +
  theme_bw()
dev.off() # Closes pdf file

#' ### Looking closer at MSCI-W
#+ echo = TRUE
# Violations of MSCI-W ---------------------------------------------
msciw <- os_risk.tbl %>% 
  filter(id_name == "MSCI-W", 
         model_type %in% c("cevt", "riskmetrics")) %>% 
  select(-ES.xts) %>% 
  mutate(var = map(VaR.xts, ~fortify(.x))) %>% 
  select(-VaR.xts) %>% 
  unnest(var) %>% 
  mutate(var = select(., Zq975:Zq990) %>% rowSums(na.rm = TRUE)) %>% 
  select(-c(Zq975, Zq990)) %>% 
  left_join(real.tbl, by = c("period", "id_name", "Index")) %>% 
  filter(period == "period1") %>% 
  mutate(viol = if_else(close > var, var, NA_real_)) %>% 
  mutate(coverage = case_when(coverage == 0.01 ~ "Coverage of 1.0%",
                              coverage == 0.025 ~ "Coverage of 2.5%"))

# Opens a pdf file to hold the plots
pdf(file = paste0(fig_path, "emec-msciw.pdf"),
    width = 7,
    height = 8,
    colormodel = "grey")
ggplot(msciw, aes(x = Index, y = var)) +
  geom_line(aes(color = model_type)) +
  geom_segment(aes(y = close, xend = Index, yend = 0), size = 0.4) +
  geom_point(aes(y = viol), shape = 4, size = 2) +
  guides(color = FALSE) +
  facet_grid(model_type ~ coverage, scales = "free") +
  labs(title = "",
       x = "",
       y = "Loss / VaR") +
  scale_y_continuous() +
  scale_x_date(date_breaks = "1 years", date_labels = "%Y") +
  theme_bw()
dev.off() # Closes pdf file

# What happens to NAREIT at 2.5% coverage and unconditional models? ---------
mtype <- c("unorm", "ut")
cover <- 0.025
asset <- "NAREIT"
real.tbl <- realized %>% 
  mutate(real = map(real, ~fortify(.x))) %>% 
  unnest(real)
nareit.tbl <- os_risk.tbl %>% 
  filter(model_type %in% mtype & coverage == cover & id_name == asset) %>% 
  filter(period == "period1") %>% 
  select(-c(ES.xts, coverage)) %>% 
  mutate(var = map(VaR.xts, ~fortify(.x))) %>% 
  select(-VaR.xts) %>% 
  unnest(var) %>% 
  spread(key = model_type, value = Zq975) %>% 
  left_join(real.tbl, by = c("period", "id_name", "Index")) %>% 
  gather(key = risk, value = var, -c(period, id_name, Index, close)) %>% 
  mutate(viol = if_else(close > var, var, NA_real_)) 

# Plot of all 6 assets for each period of selected model and coverage
ggplot(nareit.tbl, aes(x = Index, y = var, color = risk)) +
  geom_line() +
  geom_segment(aes(y = close, xend = Index, yend = 0), size = 0.5) +
  geom_point(aes(y = viol), shape = 4, size = 2) +
  facet_grid(rows = vars(risk), scales = "free") +
  labs(title = "",
       x = "",
       y = "Loss / VaR") +
  scale_y_continuous() +
  scale_x_date(date_breaks = "1 years", date_labels = "%Y") +
  theme_bw()

# Table of Statistics
nareit_risk.tbl <- os_risk.tbl %>% 
  filter(model_type %in% mtype & coverage == cover & id_name == asset) %>% 
  filter(period == "period1") %>% 
  left_join(realized, by = c("period", "id_name")) %>% 
  transmute(period = period,
            id_name = id_name,
            model_type = model_type,
            coverage = coverage,
            VaRtest = pmap(., ~vartest(..4, -coredata(..7), -coredata(..5)))) %>% 
  unnest() %>% 
  mutate(uc_stat = map2_chr(uc.LRstat, uc.LRp, ~format_stars(c(.x, .y), 
                                                             digits = 2)),
         dur_stat = map2_chr(dur.LR, LRp, ~format_stars(c(.x, .y), 
                                                        digits = 2))) %>% 
  select(coverage, period, id_name, model_type, uc_stat, dur_stat) %>%
  gather(key = stat_name, value = stat_value, 
         -c(period, id_name, coverage, model_type), 
         factor_key = TRUE) %>% 
  spread(key = id_name, value = stat_value) %>% 
  mutate(coverage = case_when(coverage == 0.01 ~ "Coverage of 1.0%",
                              coverage == 0.025 ~ "Coverage of 2.5%"))


#' ### Statistical tests for ES
#+ echo = TRUE
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

#' # Appendix of other codes
#+ echo = TRUE
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

