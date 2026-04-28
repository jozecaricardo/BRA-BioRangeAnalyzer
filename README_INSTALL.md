# BRA (BioRangeAnalyzer Shiny) - Installation Guide

## Pré-requisitos

- R 4.0 ou superior
- RStudio (recomendado)

## Instalação Passo a Passo

### 1. Instalar Dependências

Abra o R ou RStudio e execute:

```r
# Copie e cole este código no console do R:
source("install_dependencies.R")
```

Ou manualmente:

```r
install.packages(c(
  "shiny", "shinydashboard", "shinyjs", "shinyWidgets", "shinythemes",
  "leaflet", "DT", "ape", "phytools", "dplyr", "tidyr", "ggplot2", 
  "plotly", "RColorBrewer", "sp", "sf", "raster", "terra", "viridis", 
  "geosphere", "vegan", "spdep", "devtools"
))
```

### 2. Instalar o pacote BRA (BioRangeAnalyzer Shiny)

```r
devtools::install_local("biogeoshiny")
```

Ou se estiver no diretório do pacote:

```r
devtools::install_local(".")
```

### 3. Executar o App

```r
biogeoshiny::run_biogeoshiny()
```

## Formato de Dados Esperado

### Arquivo CSV/TXT de ocorrencias

Seu arquivo CSV/TXT **DEVE** ter exatamente 3 colunas:

| spp | long | lat |
|-----|------|-----|
| Homo_sapiens | -50.5 | -25.3 |
| Homo_sapiens | -51.2 | -26.1 |
| Pan_troglodytes | -48.3 | -24.5 |

**Importante:**
- Coluna `spp`: Nome da espécie (texto)
- Coluna `long`: Longitude em graus decimais (número)
- Coluna `lat`: Latitude em graus decimais (número)
- As colunas podem estar em qualquer ordem
- Nomes das colunas são case-insensitive (SPP, spp, Spp funcionam)

### Arquivo de Árvore Filogenética

Formatos aceitos:
- Newick (.nwk, .newick, .txt, .tre)
- Nexus (.nex, .nexus)

Exemplo Newick:
```
((Homo_sapiens:0.1,Pan_troglodytes:0.1):0.2,Gorilla_gorilla:0.3);
```

## Fluxo de Trabalho

1. **Data Input** - Carregue seus dados (CSV/TXT) e arvore
2. **Tree Validation** - Valide a árvore filogenética
3. **Data Preprocessing** - Remova taxas problemáticas (singletons)
4. **Range Extrapolation** - Escolha método de extrapolação:
   - Buffers circulares
   - Convex Hull (polígonos convexos mínimos)
   - MST (Minimum Spanning Tree)
5. **Visualizations** - Visualize os resultados
6. **Export Files** - Exporte em múltiplos formatos:
   - BioGeoBEARS (.data)
   - NEXUS (.nex)
   - TNT (.tnt)
   - NDM (.xyd)

## Solução de Problemas

### Erro: "Package 'X' not found"

Execute `install_dependencies.R` novamente ou instale o pacote manualmente:

```r
install.packages("nome_do_pacote")
```

### Erro: "Colunas indefinidas selecionadas"

Verifique se seu CSV tem as 3 colunas obrigatórias: `spp`, `long`, `lat`

### Erro: "Não foi possível encontrar a função"

Certifique-se de que o pacote foi instalado corretamente:

```r
devtools::install_local("biogeoshiny")
```

## Projeto e autores

Projeto: https://jozecaricardo.github.io/biogeografiaAULAS/

- Jose Ricardo Inacio Ribeiro (Universidade Federal do Pampa (UNIPAMPA), campus Sao Gabriel, Rio Grande do Sul State, Brazil)
- Augusto Ferrari (Universidade Federal do Rio Grande (FURG), Rio Grande State, Brazil)
