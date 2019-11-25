# roll_fit_unorm ----------------------------------------------------------
# Ajusta os dados para um modelo Normal incondicional
mc_roll_fit_unorm <- function(data, n.roll, window.size) {
  tic <- Sys.time()
  tmp.list <- mclapply(1:n.roll, function(i){
    xts <- data[i:(window.size+i)]
    mean <- mean(xts)
    sd <- sd(xts)
    Zq975 <- qnorm(0.975, mean, sd)
    Zq990 <- qnorm(0.990, mean, sd)
    # Equacao do ES retirada de Pfaff2013 p. 36, eq. 4.5
    # ES_a = 1/(1-a) * int_a^1 q_l(x)dx
    Sq975 <- integrate(function(x){
      qnorm(x, mean = mean, sd = sd)},
      0.975,
      1)$value / (1-0.975)
    Sq990 <- integrate(function(x){
      qnorm(x, mean = mean, sd = sd)},
      0.990,
      1)$value / (1-0.990)
    return(cbind(mean, sd, Zq975, Zq990, Sq975, Sq990))
  },
  mc.cores = cores) # Fim do lapply
  # ao final da iteracao teremos uma lista com n.roll elementos, cada um correspondente
  # a uma data onde foi feito o ajuste dos dados. Nas colunas de cada elemento da lista estao 
  # os parametros e medidas de risco estimados para cada um dos dias fora da amostra
  # Junta-se tudo por rbind e joga fora o ultimo elemento, pois nao tem 
  # perda realizada para comparar com.
  # Depois forma um xts indexado pelos dias fora da amostra a partir do segundo dia
  ans <- do.call(rbind, tmp.list)
  ans <- ans[-dim(ans)[1],]           
  ans <- xts(ans, order.by = index(data[(window.size+2):(window.size+n.roll)]))
  colnames(ans) <- c("mu", "sigma", "Zq975", "Zq990", "Sq975", "Sq990")
  # Preenche os NA com a ultima observacao conhecida
  ans <- na.locf(ans)
  toc <- Sys.time()
  #print(toc-tic)
  return(ans)
}

# roll_fit_unorm ----------------------------------------------------------
# Ajusta os dados para um modelo Normal incondicional
or_roll_fit_unorm <- function(data, n.roll, window.size) {
  tic <- Sys.time()
  tmp.list <- mclapply(1:n.roll, function(i){
    xts <- data[i:(window.size+i)]
    mean <- mean(xts)
    sd <- sd(xts)
    Zq975 <- qnorm(0.975, mean, sd)
    Zq990 <- qnorm(0.990, mean, sd)
    # Equacao do ES retirada de Pfaff2013 p. 36, eq. 4.5
    # ES_a = 1/(1-a) * int_a^1 q_l(x)dx
    Sq975 <- integrate(function(x){
      qnorm(x, mean = mean, sd = sd)},
      0.975,
      1)$value / (1-0.975)
    Sq990 <- integrate(function(x){
      qnorm(x, mean = mean, sd = sd)},
      0.990,
      1)$value / (1-0.990)
    return(cbind(mean, sd, Zq975, Zq990, Sq975, Sq990))
  },
  mc.cores = cores) # Fim do lapply
  # ao final da iteracao teremos uma lista com n.roll elementos, cada um correspondente
  # a uma data onde foi feito o ajuste dos dados. Nas colunas de cada elemento da lista estao 
  # os parametros e medidas de risco estimados para cada um dos dias fora da amostra
  # Junta-se tudo por rbind e joga fora o ultimo elemento, pois nao tem 
  # perda realizada para comparar com.
  # Depois forma um xts indexado pelos dias fora da amostra a partir do segundo dia
  ans <- do.call(rbind, tmp.list)
  ans <- ans[-dim(ans)[1],]           
  ans <- xts(ans, order.by = index(data[(window.size+2):(window.size+n.roll)]))
  colnames(ans) <- c("mu", "sigma", "Zq975", "Zq990", "Sq975", "Sq990")
  # Preenche os NA com a ultima observacao conhecida
  ans <- na.locf(ans)
  toc <- Sys.time()
  #print(toc-tic)
  return(ans)
}


bench <- microbenchmark(mc_roll_fit_unorm(ibov.xts, n.roll, window.size),
               or_roll_fit_unorm(ibov.xts, n.roll, window.size))
autoplot(bench)
boxplot(bench)
