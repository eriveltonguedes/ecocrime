# ==============================================================================
# MODELO LOGIT HIERÁRQUICO: DETERMINANTES DA ENTRADA NO CRIME ORGANIZADO
# Integração de Microdados Criminais (IPEA, 2023) e Indicadores Municipais (2022)
# ==============================================================================

# 1. PACOTES NECESSÁRIOS -------------------------------------------------------
pacotes <- c("readxl", "dplyr", "lme4", "performance", "broom.mixed", "marginaleffects")
novos_pacotes <- pacotes[!(pacotes %in% installed.packages()[, "Package"])]
if(length(novos_pacotes)) install.packages(novos_pacotes)

library(readxl)
library(dplyr)
library(lme4)
library(performance)
library(broom.mixed)
library(marginaleffects)

# 2. IMPORTAÇÃO DAS BASES DE DADOS ---------------------------------------------
# Insira o nome exato do arquivo Excel (deve estar no diretório de trabalho)
arquivo_excel <- "base_dados_monografia.xlsx"

# Planilha 1: Microdados individuais dos réus (Nível 1)
# Nota: A coluna "cod_mun" deve estar preenchida com o código IBGE (7 dígitos)
df_ipea <- read_excel(arquivo_excel, sheet = "Variáveis IPEA")

# Planilha 2: Indicadores socioeconômicos e de infraestrutura municipal (Nível 2)
# Nota: Caso o nome da aba no Excel seja diferente de "Municípios", ajuste abaixo
df_mun <- read_excel(arquivo_excel, sheet = "Municípios")

# 3. PREPARAÇÃO E PAREAMENTO (MERGE) -------------------------------------------
# Padronização da chave municipal para cruzamento seguro
df_ipea <- df_ipea %>% mutate(cod_mun = as.character(cod_mun))
df_mun  <- df_mun  %>% mutate(cod_mun = as.character(cod_mun))

# Tratamento de missing em "raça":
# Como há células vazias na variável "raça" (capturadas pela dummy "missing_raça"),
# atribuímos 0 aos NAs de "raça" para evitar que o R descarte essas observações via 
# listwise deletion. Assim, o efeito da não declaração fica isolado em "missing_raça".
df_ipea <- df_ipea %>%
  mutate(raça_tratada = ifelse(is.na(raça), 0, raça))

# Junção Nível 1 (indivíduos) + Nível 2 (municípios)
df_analise <- df_ipea %>%
  inner_join(df_mun, by = "cod_mun")

# Padronização (Z-score) das variáveis contínuas para garantir convergência numérica
# do algoritmo de máxima verossimilhança (Laplace) no modelo multinível
df_analise <- df_analise %>%
  mutate(across(c(idade_referência, pop, area_urb, area_equip, Gini_2010, 
                  pop_urb, grau_urb, tx_alfab, tx_hom_2019, educ_percapita, 
                  assist_percapita, dens, evasao, renda, bolsa_milhab),
                ~ as.numeric(scale(.)),
                .names = "{.col}_z"))

# 4. ESTIMAÇÃO DOS MODELOS MULTINÍVEL (GLMM) -----------------------------------

# 4.1 Modelo Nulo (Apenas intercepto aleatório por município para cálculo do ICC)
modelo_nulo <- glmer(
  facção ~ 1 + (1 | cod_mun),
  data = df_analise,
  family = binomial(link = "logit"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

# 4.2 Modelo Completo (Nível 1 + Nível 2)

modelo_completo <- glmer(
  facção ~ 
    # --- Nível 1: Características Individuais e Processuais ---
    sexo + idade_referência_z + raça_tratada + missing_raça +
    fund + méd + sup +
    desemp + estud + aut +
    D_cond + D_org + D_armas + D_drogas +
    D_rein + D_rein_org + D_rein_armas + D_rein_drogas +
    
    # --- Nível 2: Infraestrutura e Socioeconomia Municipal ---
    area_urb_z + area_equip_z + grau_urb_z + dens_z +
    Gini_2010_z + renda_z + bolsa_milhab_z +
    tx_alfab_z + evasao_z + educ_percapita_z + assist_percapita_z +
    tx_hom_2019_z +
    
    # --- Efeito Aleatório Municipal (Intercepto Aleatório) ---
    (1 | cod_mun),
  data = df_analise,
  family = binomial(link = "logit"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

# 5. RESULTADOS E MÉTRICAS DE AJUSTE -------------------------------------------
# Resumo das estimativas (Log-Odds)
summary(modelo_completo)

# Coeficiente de Correlação Intraclasse (ICC)
icc_nulo     <- performance::icc(modelo_nulo)
icc_completo <- performance::icc(modelo_completo)
print(icc_nulo)
print(icc_completo)

# Tabela de Odds Ratios (Razões de Chances) com Intervalo de Confiança de 95%
tabela_or <- broom.mixed::tidy(
  modelo_completo, 
  effects = "fixed", 
  conf.int = TRUE, 
  exponentiate = TRUE
)
print(tabela_or, n = Inf)

# Efeitos Marginais Médios (AME - Average Marginal Effects) em pontos percentuais
ame_resultados <- marginaleffects::avg_slopes(modelo_completo)
print(ame_resultados)

# 6. EXPORTAÇÃO DOS RESULTADOS AGREGADOS ---------------------------------------
# Exporta apenas as tabelas de parâmetros estimados (sem expor microdados individuais)
write.csv(tabela_or, "resultados_odds_ratios.csv", row.names = FALSE)
write.csv(as.data.frame(ame_resultados), "resultados_efeitos_marginais.csv", row.names = FALSE)
