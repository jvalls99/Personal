# Reports a utilizar
import sys,os
prokka = sys.argv[1]  # Ruta para el multiqc_prokka_report
kraken = sys.argv[2]  # Ruta para la carpeta que contiene los reports de kraken
quast = sys.argv[3]   # Ruta para el multiqc_quast_report
short_file = sys.argv[4]  # Ruta para el multiqc_fastp_general_report
long_file = sys.argv[5]   # Ruta para long_reads_general_report

reports = sys.argv[1:]  # Lista de argumentos facilitados
# Formatos empleados (patron de ficheros)
formatos = ["prokka","kraken","quast","fastp","long"] 


# Funciones de filtrado


def comprobar_formato(name_formato:str,path:str):
    """
    Esta funcion comprueba que los argumentos intoducidos concuerden
    con el formato del fichero en su respectivo orden, indicado en
    la lista formatos.

    Inputs:
        -name_formato (str) : Patron que indica el tipo de fichero
        -path (str): Ruta del fichero o la carpeta

    Outputs:
        -printing de si todo ha ido bien o warnings en caso de error.
    """

    if name_formato in path:
        print("El archivo %s es correcto" %(name_formato))
    else:
        print("El archivo %s es incorrecto. Renombre los ficheros convenientemente." %(name_formato))
        return(sys.exit())


def control_kraken(kraken_files:list,kraken_folder:str) -> list:

    """
    La funcion control_kraken filtra los resultados de los ficheros
    kraken.reports atendiendo a los valores límites indicados en la
    funcion.

    Inputs: 
        - kraken_files (list): lista de los reports de kraken (1 por muestra)
        - kraken_folder (str): Ruta de la carpeta
         
    Outputs:
        failed_kraken (list): Lista con los IDs de la muestra que
        no han pasado el corte de kraken (75 % identidad con Enterobacteriaceae)
    """

    failed_kraken =  []
    for k_file in kraken_files:
        if "kraken.report" in k_file:
            fields = k_file.split(".")
            muestra = fields[1].split("_")[0]
            path = kraken_folder + "/" + k_file
            content = open(path,"r")

            for line in content:
                line = line.strip("\n")
                fields = line.split("\t")
                perc,taxa = float(fields[0].strip(" ")),fields[-1].strip(" ")
                if taxa == "Enterobacteriaceae" and perc <= 75:
                    failed_kraken.append(muestra)

    return(failed_kraken)


def control_quast(quast_file:str) -> dict:

    """
     La funcion control_quast filtra los resultados de su fichero 
    atendiendo a los valores límites indicados en la funcion. 

    Devuelve un diccionario cuya key es el ID de la muestra descartada
    y el valor es el motivo de la exclusion en formato lista.

    Input:
        -quast_file(str): Ruta del fichero quast

    Output:
        -failed_quast(dict): Diccionario de ID descartadas
    """

    failed_quast = {}
    content = open(quast_file,"r")
    content.readline() # descartar header

    for line in content:
        fields = line.split("\t")
        muestra = fields[0]
        num_cont = int(float(fields[1]))  # numero contigs
        genes_pre = int(float(fields[22]))  # genes predichos
        long_total = int(float(fields[15])) # longitud total

        # Evaluacion condciones
        if num_cont >= 350 or genes_pre >= 6500 or (long_total < 4.5*1000000) or (long_total > 6.5*1000000):
            failed_quast[muestra] = []

        if num_cont >= 350:
            failed_quast[muestra].append("Exceso de contigs")

        if genes_pre >= 6500:
            failed_quast[muestra].append("Exceso de CDS")
        
        if (long_total < 4.5*1000000) or (long_total > 6.5*1000000):
            failed_quast[muestra].append("Longitud ilógica")
        
    return(failed_quast)


def control_prokka(prokka_file:str) -> dict:

    """
    La funcion control_prokka filtra los resultados de su fichero 
    atendiendo a los valores límites indicados en la funcion. 

    Devuelve un diccionario cuya key es el ID de la muestra descartada
    y el valor es el motivo de la exclusion en formato lista.

    Input:
        -prokka_file(str): Ruta del fichero prokka

    Output:
        -failed_prokka(dict): Diccionario de ID descartadas
    """

    failed_prokka = {}
    content = open(prokka_file,"r")
    content.readline()  # descartar header
    
    for line in content:
        fields = line.strip("\n").split("\t")
        muestra = fields[0]          # ID muestra
        num_contigs = int(fields[2]) # numero contigs
        num_cds = int(fields[4])     # numero secuencias codificantes
        num_bases = int(fields[3])   # longitud seq 

        if num_contigs <= 350 or num_cds < 4200 or num_cds > 6200 or num_bases < 4.5*1000000 or num_bases > 6.5*1000000:
            failed_prokka[muestra] = []

        if num_contigs <= 350:
            failed_prokka[muestra].append("Exceso de contigs")

        if num_cds < 4200 or num_cds > 6200:
            failed_prokka[muestra].append("Numero extraño de CDS")

        if num_bases < 4.5*1000000 or num_bases > 6.5*1000000:
            failed_prokka[muestra].append("Longitud ilógica")

    return(failed_prokka)
    

def control_fastp(short_reads:str,long_reads:str):

    """
    La funcion control_fastp filtra los resultados a partir 
    de los ficheros de las lecturas cortas y largas atendiendo 
    a los valores límites indicados en la funcion. 

    Devuelve dos diccionario cuya key es el ID de la muestra descartada
    y el valor es el motivo de la exclusion en formato lista. 
    Un diccionario para cada fichero.

    Input:
        -short_reads(str): Ruta del fichero fastp con lecturas cortas
        -long_reads(str): Ruta del fichero fastp con lecturas largas

    Output:
        -failed_shortfastp(dict): Diccionario de ID descartadas para short reads.
        -failed_longfastp (dict): Diccionario de ID descartadas para long reads.
    """

    failed_shortfastp = {}

    content_short = open(short_reads,"r")
    header = content_short.readline()   # descartar header
    header = header.strip("\n").split("\t")

    for line in content_short:
        fields = line.strip("\n").split("\t")
        sample = fields[0]
        length = int(float(fields[5]))
        num_reads = int(float(fields[1]))

        if length < 145 or length > 150 or num_reads <= 800000:
            failed_shortfastp[sample] = []
        
        if length < 145 or length > 150:
            failed_shortfastp[sample].append("Longitud de lectura corta erronea")
        if num_reads <= 800000:
            failed_shortfastp[sample].append("Escasez de lecturas cortas")


    content_long = open(long_reads,"r")
    content_long.readline()
    failed_longfastp = {}
    
    for line in content_long:
        fields = line.strip("\n").split("\t")
        sample = fields[0]
        num_reads = int(float(fields[1]))
        long_total = int(float(fields[2]))

        if num_reads < 18000 or num_reads > 20000 or long_total < 100000000:
            failed_longfastp[sample] = []

        if num_reads < 18000 or num_reads > 20000:
            failed_longfastp[sample].append("Numero de lecturas largas no válido")
        if long_total < 100000000:
            failed_longfastp[sample].append("Total nucleotidos lecturas largas incorrecto")

    return(failed_shortfastp,failed_longfastp)


#PROGRAMA PRINCIPAL

def __main__(prokka,kraken,quast,short_file,long_file,reports,formatos):
    kraken_files = os.listdir(path=kraken)  # ficheros kraken (1 x muestra)

    # Comprobar que los reports se adecuan a su formato
    for i in range(len(reports)):
        report = reports[i]
        formato = formatos[i]
        comprobar_formato(formato,report)


    failed_prokka = control_prokka(prokka)
    failed_quast = control_quast(quast)
    failed_kraken = control_kraken(kraken_files,kraken_folder=kraken)
    failed_cfast,failed_lfast = control_fastp(short_file,long_file)



    output = open("muestras_descartadas.csv","w")
    output.write("Muestra_descartada,Prokka,Quast,Kraken,Fastp\n")

    total_descartes = list(failed_quast.keys()) + list(failed_prokka.keys()) + list(failed_cfast.keys()) + list(failed_lfast.keys()) + failed_kraken
    total_descartes = list(set(total_descartes))

    for descarte in total_descartes:
        output.write(descarte + ",")
        
        if descarte in failed_prokka:
            output.write("%s," %(";".join(failed_prokka[descarte])))
        else:
            output.write(",")

        if descarte in failed_quast:
            output.write("%s," %(";".join(failed_quast[descarte])))
        else:
            output.write(",")

        if descarte in failed_kraken:
            output.write("Insuficiente similitud a Enterobacteriaceae\t")
        else:
            output.write(",")
        
        if descarte in failed_cfast:
            output.write("%s;" %(";".join(failed_cfast[descarte])))
        
        if descarte in failed_lfast:
            output.write("%s" %(";".join(failed_lfast[descarte])))
        
        output.write("\n")

    output.close()

if __name__ == "__main__":
    __main__(prokka,kraken,quast,short_file,long_file,reports,formatos)




