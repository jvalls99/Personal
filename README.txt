Esta carpeta contiene los programas necesarios para la construcción de un core genome y un pangenoma.

1. construir_core.sh: Recibe una lista de identificadores de NCBI (GCF) como argumento. El usuario elige mediante parámetros si desea construir el core o pangenoma, si desea elaborar una base de datos plasmídica y en ese caso, puede introducir el fichero FASTA con las secuencias. Se automatizan los procesos de
descarga, anotación, eliminación de plásmidos (únicamente para core) y obtención de los alineamientos con Panaroo. También se obtiene los SNPs y la filogenia
core.

  Entrada:
        - lista_accesiones: Un ID de NCBI por linea.
  Salidas:
        - Carpeta con los genomas descargados y sus ficheros de anotación e identificación de plásmidos-
        - Carpeta con los alineamientos individuales de genes y su concatenación (core/pangenoma)
        - Árbol filogenético con los SNPs del alineamiento.
  
            
