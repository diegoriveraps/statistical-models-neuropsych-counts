## Modelos lineales y modelos lineales generalizados con distribución Poisson aplicados a datos de conteo en neuropsicología: un análisis comparativo

### Objetivo

El objetivo de este estudio fue evaluar el comportamiento estadístico de las puntuaciones de tipo conteo, obtenidas en un tiempo definido, en un conjunto de pruebas neuropsicológicas de fluidez verbal y velocidad de procesamiento. Se compararon modelos lineales (LM) y modelos lineales generalizados con distribución de Poisson (GLM) para determinar qué metodología ofrece estimaciones más adecuadas a la naturaleza discreta y no negativa de los datos y una interpretación clínica coherente.

### Metodología

1. Datos: 11 puntuaciones neuropsicológicas de participantes mexicanos.
2. Ajustes del modelo:

   * Edad: polinomio de orden 2.
   * Escolaridad: transformación logarítmica.
   * Sexo e interacciones.

3. Selección de variables: inferencia bayesiana y selección exhaustiva.
4. Estimación: MCMC con JAGS.
5. Validación: validación cruzada
6. Generación de normas: cálculo de percentiles normativos para cada prueba.
7. Curvas ROC (controles sanos vs. Alzheimer).

   Estructura del repositorio

### Estructura del repositorio

* scripts/: códigos R organizados por etapas:

  * 01\_libraries\_preprocess.R
  * 02\_data\_analysis.R
  * 03\_cross\_validation.R
  * 04\_generate\_norms.R
  * 05\_plots\_tables.R
  
* Material suplementario (Trace & Density plots)

  ### Requisitos

* R 4.0 o superior.
* Paquetes R: tidyverse, rjags, coda, future.apply, caret, patchwork, cowplot, pROC, readxl.
* JAGS instalado en el sistema.

  ### 

