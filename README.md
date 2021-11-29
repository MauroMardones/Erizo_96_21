------------------------------------------------------------------------

------------------------------------------------------------------------

# PBRs\_Erizo 2021 serie 1996-2021

Este repositorio contiene los análisis realizados para el cálculo de PBRs para las unidades de stock del Erizo.

**Descripción del contenido de la carpeta "codigos\_admb"**

El código base corresponde a un Modelo de Captura a la talla ("Modbento") creado por cristian.canales\@ifop.cl. Los códigos utilizados para este estudio corresponde a una modificación del "Modbento". Esta modificación es realizada por mariajose.zuniga\@ifop.cl y mauricio.mardones\@ifop.cl se detallará a continuación.

Contiene los códigos de ADMB principalmente MAET.tpl y MAET.dat del modelo de evaluación de stock de erizo de la zona X Norte (MAETXN) , X Sur (MAETXS) y XI (MAETXI).

Los modelos corresponden a modelos anuales (año calendario) con dinámica a la edad e información en tallas (utiliza matriz de probabilidad edad\_talla).

Los principales datos que ingresan al modelo corresponden a los desembarques, CPUE y frecuencia de tallas de la flota. Utiliza un vector de pesos medios a la talla, madurez a la talla, y parámetros de crecimiento y mortalidad natural.

## 

**CARPETAS:**

1.  **codigos\_admb**

-   código fuente (.tpl)

-   datos y controles (.dat)

2.  **funciones**

-   Fn\_Retrospectivo.R = Retrospectivo (genera la carpeta con archivos admb)

-   functions.R = lisread (función que lee archivo ".dat") y reptoRlist (función que lee archivo ".rep")

**Códigos en Rmarkdown**

-   Informe\_Estatus\_Erizo\_word.Rmd = corresponde al informe de estatus del las tres unidades de stock del erizo que debe ser entregado entre diciembre y enero.

-   Erizo.Rmd = contendrá los análisis exploratorios y cálculo de PBRs

-   

**Archivos:**

-   Informe\_Estatus\_Erizo\_word.docx = es la salida del código .Rmd con el mismo nombre

-   TemplateWord\_Erizo.docx = es utilizado en el código "Informe\_Estatus\_Erizo\_word.Rmd" para darle el formato deseado (tipo de letra, títulos, etc).

-   Readme\_codADMB.pdf = se describen los archivos de datos y controles (.dat), reporte (.rep), el código fuente (.tpl), los parámetros estimados por el modelo (.par) y parámetros y variables estimadas con error (.std).
