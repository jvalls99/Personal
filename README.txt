Esta carpeta contiene los programas necesarios para la construcción de un core genome y un pangenoma.

TFM:

1. construir_core.sh: Recibe una lista de identificadores de NCBI (GCF) como argumento. El usuario elige mediante parámetros si desea construir el core o pangenoma, si desea elaborar una base de datos plasmídica y en ese caso, puede introducir el fichero FASTA con las secuencias. Se automatizan los procesos de
descarga, anotación, eliminación de plásmidos (únicamente para core) y obtención de los alineamientos con Panaroo. También se obtiene los SNPs y la filogenia
core.

  Entrada:
        - lista_accesiones: Un ID de NCBI por linea.
  Salidas:
        - Carpeta con los genomas descargados y sus ficheros de anotación e identificación de plásmidos-
        - Carpeta con los alineamientos individuales de genes y su concatenación (core/pangenoma)
        - Árbol filogenético con los SNPs del alineamiento.
  
            
EXPRESIÓN DIFERENCIAL:

1. pipeline.nf: Está escrito con lenguaje Nextflow para facilitar la construcción de la pipeline. Permite obtener a partir de identificadores SRA y una referencia los ficheros BAM, a partir de los cuales se va a construir la matriz de conteo y finalmente estudiar la expresión Diferencial con edgeR y DeSeq2.

2. build_matrix.R: Construye la matriz y ejecuta los análisis de Expresión Diferencial.


BÚSQUEDA DE PATRONES Y MANIPULACIÓN DE SEUCENCIAS VÍRICAS.

1. search_pattern_virus.py
