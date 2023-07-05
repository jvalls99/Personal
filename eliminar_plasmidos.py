import sys

fichero_entrada = sys.argv[1]
fichero_salida = sys.argv[2]
cabeceras = sys.argv[3:]

contenido_fasta = open(fichero_entrada,"r")
salida = open(fichero_salida,"w")
for line in contenido_fasta:
    line = line.strip("\n")
    if line.startswith(">"):
        header = line[1:].split(" ")[0]
        if header not in cabeceras:
            REMOVE = False
            salida.write(line + "\n")
        else:
            REMOVE = True
    else:
        if REMOVE is False:
            salida.write(line + "\n")
salida.close()

