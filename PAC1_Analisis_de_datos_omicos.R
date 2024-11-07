
### Ruta de gestión en archivo  
setwd("C:/Users/ASUS/OneDrive/Escritorio/Maestria en Bioinformatica y Bioestadistica/6. Analisis de Datos Omicos/PAC1")

### Con gert podemos obtener datos o clona e repositorios en Github

if (!require("gert")) {
  install.packages("gert")
  library(gert)
}

### Script para obtención de datos o repositorios

git_clone("https://github.com/nutrimetabolomics/metaboData.git",
          path = "C:/Users/ASUS/OneDrive/Escritorio/Maestria en Bioinformatica y Bioestadistica/6. Analisis de Datos Omicos/PAC1")


### Instalación de BiocManager para el uso de SummarizedExperiment


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SummarizedExperiment")
BiocManager::install("Biobase")

### Configuración de SummarizedExperiment con los datos de 2018-MetabotypingPaper/DataValues_S013.csv


library(SummarizedExperiment)

# Cargar como objeto archivo CSV con datos de expresión, demogragráfico antropomórfico y metabolomica
data <- read.csv("C:/Users/ASUS/OneDrive/Escritorio/Maestria en Bioinformatica y Bioestadistica/6. Analisis de Datos Omicos/PAC1/Datasets/2018-MetabotypingPaper/DataValues_S013.csv", row.names = 1)

# Crear un objeto SummarizedExperiment
# Separar los datos en 'assay' y 'colData' según el formato de tu dataset
assay_data <- as.matrix(data[, -1])  # Matriz con datos de expresión, demogragráfico, de medicación y metabolitos paraclinicos

sample_info <- data.frame(sampleID = colnames(assay_data))  # Información de muestras

# Crear el contenedor SummarizedExperiment
se <- SummarizedExperiment(assays = list(counts = assay_data),
                           colData = sample_info)


# resumen de datos
se

# Se lista las variables y promeros datos
head(assay(se))

# Resumen de ensayo
summary(assay(se))

# Caracteristica de las variables
colData(se)

##################################################

# Analisis de datos: a continuación se realizará una comparación de datos antropomórficos y de metabolitos 
# de los dos grupos uno metabolicamente sano (Grupo 1) y otro metabolicamente no sano (Grupo 2) y se compara
# con t student antes de inicio de la intervención bariátrica y seis meses despues con el fin de evaluar si hay
# cambios significativos despues de la intervención durante ese periodo 


# Para este analisis se filyro datos por grupo
group1 <- subset(data, Group == 1)
group2 <- subset(data, Group == 2)

# Variables T0 y T5 a comparar (solo se tomó antropomorficos y de metabolitos para este ejemplo donde TO es la
# linea base y T5 son las mediciiones despues de seis meses)

variables_T0 <- c("MEDDM_T0", "MEDCOL_T0", "MEDINF_T0", "MEDHTA_T0",
                  "GLU_T0", "INS_T0", "HOMA_T0", "HBA1C_T0", "HBA1C.mmol.mol_T0",
                  "PESO_T0", "bmi_T0", "CC_T0", "CINT_T0", "CAD_T0", "TAD_T0",
                  "TAS_T0", "TG_T0", "COL_T0", "LDL_T0", "HDL_T0", "VLDL_T0",
                  "PCR_T0", "LEP_T0", "ADIPO_T0", "GOT_T0", "GPT_T0", "GGT_T0",
                  "URICO_T0", "CREAT_T0", "UREA_T0", "HIERRO_T0", "TRANSF_T0", "FERR_T0")

variables_T5 <- c("MEDDM_T5", "MEDCOL_T5", "MEDINF_T5", "MEDHTA_T5",
                  "GLU_T5", "INS_T5", "HOMA_T5", "HBA1C_T5", "HBA1C.mmol.mol_T5",
                  "PESO_T5", "bmi_T5", "CC_T5", "CINT_T5", "CAD_T5", "TAD_T5",
                  "TAS_T5", "TG_T5", "COL_T5", "LDL_T5", "HDL_T5", "VLDL_T5",
                  "PCR_T5", "LEP_T5", "ADIPO_T5", "GOT_T5", "GPT_T5", "GGT_T5",
                  "URICO_T5", "CREAT_T5", "UREA_T5", "HIERRO_T5", "TRANSF_T5", "FERR_T5")

# Crear data frame para los resultados
results <- data.frame(Variable = variables_T0,
                      Group1_Mean_T0 = numeric(length(variables_T0)),
                      Group1_SD_T0 = numeric(length(variables_T0)),
                      Group2_Mean_T0 = numeric(length(variables_T0)),
                      Group2_SD_T0 = numeric(length(variables_T0)),
                      Group1_Mean_T5 = numeric(length(variables_T5)),
                      Group1_SD_T5 = numeric(length(variables_T5)),
                      Group2_Mean_T5 = numeric(length(variables_T5)),
                      Group2_SD_T5 = numeric(length(variables_T5)),
                      p_value_group1 = numeric(length(variables_T0)),
                      p_value_group2 = numeric(length(variables_T0)))

# Calcular estadísticas y pruebas t
for (i in 1:length(variables_T0)) {
  var_T0 <- variables_T0[i]
  var_T5 <- variables_T5[i]
  
  # Calcular media y desviación estándar para el grupo 1 en T0 y T5
  results$Group1_Mean_T0[i] <- mean(group1[[var_T0]], na.rm = TRUE)
  results$Group1_SD_T0[i] <- sd(group1[[var_T0]], na.rm = TRUE)
  results$Group1_Mean_T5[i] <- mean(group1[[var_T5]], na.rm = TRUE)
  results$Group1_SD_T5[i] <- sd(group1[[var_T5]], na.rm = TRUE)
  
  # Calcular media y desviación estándar para el grupo 2 en T0 y T5
  results$Group2_Mean_T0[i] <- mean(group2[[var_T0]], na.rm = TRUE)
  results$Group2_SD_T0[i] <- sd(group2[[var_T0]], na.rm = TRUE)
  results$Group2_Mean_T5[i] <- mean(group2[[var_T5]], na.rm = TRUE)
  results$Group2_SD_T5[i] <- sd(group2[[var_T5]], na.rm = TRUE)
  
  # Realizar prueba t entre T0 y T5 para cada grupo
  results$p_value_group1[i] <- t.test(group1[[var_T0]], group1[[var_T5]], paired = TRUE)$p.value
  results$p_value_group2[i] <- t.test(group2[[var_T0]], group2[[var_T5]], paired = TRUE)$p.value
}

# Resultados
print(results)



# Como resultados en el grupo de metabolicamente sanos disminuyó la medicación, antropomorfica y los metabolitos asociados a:

#MEDCOL: Medición de colesterol o diagnóstico relacionado con colesterol.
#MEDINF: Medición médica de infección o inflamación.
#GLU: Glucosa en sangre (nivel de glucosa en ayunas).
#INS: Insulina (nivel de insulina en sangre).
#HOMA: HOMA (Modelo de evaluación de homeostasis), fórmula que se usa para estimar la resistencia a la insulina.
#HBA1C: Hemoglobina A1c (indicador de los niveles promedio de glucosa en sangre durante los últimos 2-3 meses).
#HBA1C, mmol/mol_T0: Hemoglobina A1c expresada en unidades de mmol/mol.
#PESO: Peso corporal.
#BMI: Índice de Masa Corporal (IMC), medida del peso corporal en relación con la altura.
#CINT: Circunferencia de la cintura, utilizada para evaluar la distribución de la grasa corporal.
#CAD: Enfermedad arterial coronaria (CAD, por sus siglas en inglés).
#TG: Triglicéridos (tipo de grasa en la sangre).
#COL: Colesterol (total o alguna fracción, como LDL o HDL).
#PCR: Proteína C reactiva (marcador de inflamación en el cuerpo).
#LEP: Leptina (hormona relacionada con el control del apetito y la energía).
#ADIPO: Adiponectina (hormona que influye en el metabolismo de la glucosa y los lípidos).
#GGT: Gamma-glutamil transferasa (enzima hepática, utilizada como indicador de daño hepático).
#URICO: Ácido úrico (compuesto relacionado con el metabolismo de las purinas).
#UREA: Urea (producto de desecho del metabolismo de las proteínas que se elimina en la orina).
#TRANSF: Transferrina (proteína que transporta hierro en la sangre).


# Como resultados en el grupo de metabolicamente no sanos disminuyó la medicación, antropomorfica y los metabolitos asociados a:


#MEDDM: Medicación diabetes mellitus (medicación para diabetes tipo 2) al inicio.
#MEDCOL: Medicación para colesterol al inicio.
#GLU: Glucosa en sangre al inicio.
#INS: Insulina al inicio.
#HOMA: Índice de Homeostasis para la Resistencia a la Insulina al inicio.
#PESO: Peso corporal al inicio.
#BMI: Índice de masa corporal (Body Mass Index) al inicio.
#CINT: Circunferencia de la cintura al inicio.
#CAD: Coronary Artery Disease (Enfermedad Arterial Coronaria) al inicio.
#TG: Triglicéridos en sangre al inicio.
#COL: Colesterol total al inicio.
#LDL: Colesterol LDL (lipoproteínas de baja densidad) al inicio.
#PCR: Proteína C reactiva (indicador de inflamación) al inicio.
#LEP: Leptina (hormona relacionada con la regulación del peso corporal) al inicio.
#GPT: Alanina aminotransferasa (enzima hepática) al inicio.
#GGT: Gamma-glutamil transferasa (enzima hepática) al inicio.
#URICO: Ácido úrico en sangre al inicio.
#CREAT: Creatinina en sangre, indicador de función renal al inicio.
#TRANS: Transferrina (proteína que transporta hierro) al inicio.

# Conclusión: Existe efectos positivos del procedimiento de la cirujia bariatrica tanto en pacientes metabolicamente
# estable como los que no , lo que sugiere que este proceso es util tanto para disminucion de pesa y de masa como 
# mejoramoento de los paraclinicos y enfermedades cardiovasculares