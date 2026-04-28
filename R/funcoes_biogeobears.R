
# Funções auxiliares para o BioGeoBEARS App

# Função para formatar resultados
format_results <- function(results) {
  # Implementação futura
  return(results)
}

# Função para calcular teste de razão de verossimilhança
calc_LRT <- function(res1, res2) {
  # res1 deve ser o modelo mais simples (aninhado)
  # res2 deve ser o modelo mais complexo
  
  LnL1 <- res1$total_loglikelihood
  LnL2 <- res2$total_loglikelihood
  
  # Diferença em número de parâmetros
  df <- length(res2$inputs$BioGeoBEARS_model_object@params_table[res2$inputs$BioGeoBEARS_model_object@params_table$type == "free", "type"]) - 
        length(res1$inputs$BioGeoBEARS_model_object@params_table[res1$inputs$BioGeoBEARS_model_object@params_table$type == "free", "type"])
  
  # Estatística de teste
  chisq <- 2 * (LnL2 - LnL1)
  
  # Valor p
  pval <- pchisq(chisq, df = df, lower.tail = FALSE)
  
  # Retornar resultados
  result <- list(
    model1 = res1$inputs$description,
    model2 = res2$inputs$description,
    LnL1 = LnL1,
    LnL2 = LnL2,
    df = df,
    chisq = chisq,
    pval = pval,
    significant = pval < 0.05
  )
  
  return(result)
}

# Outras funções auxiliares seriam adicionadas aqui

