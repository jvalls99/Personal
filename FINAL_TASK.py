###############################################################################
## Este programa recopila la clase fasta y sus métodos asociados como parte  ##
## de la tarea final de la asignatura de Programación.                       ##
##                                                                           ##
## Versión: 2.0                                                              ##
##                                                                           ##
## Alumno: Francisco Javier Valls                                            ##
##                                                                           ##
############################################################################### 

# Módulos a importar
import sys
import re
import numpy as np
from matplotlib import pyplot as plt


fasta = sys.argv[1]          #fichero fasta: input

# 1. Generación de la clase Fasta

class Fasta:

    def __init__(self, header, secuencia):
        self.header = header    # la header del fichero fasta
        self.seq = secuencia    # las lineas de secuencia del fichero FASTA
    
    
    def evaluar_secuencia(self) -> str:
        """
        Esta función devuelve un mensaje donde se indica si la secuencia input
        (str) es DNA o no.

        Inputs: Objeto tipo FASTA

        output: Str
        """
        secuencia = self.seq
        # conjunto de letras que se permiten en la secuencia FASTA
        bases_permitidas = ['N','n','a','c','g','t','A','C','G','T', '\n']      
        numero = 0
        for caracter in bases_permitidas:
            numero += secuencia.count(caracter)

        if numero != len(secuencia):
            return(False)
        else:
            return(True)
    

    def constructor(file) -> "Fasta":
        """
        Esta función a partir de un fichero FASTA, devuelve un objeto de 
        la clase Fasta.

        Precondicion: No se aceptan ficheros MULTIFASTA, es decir, solo puede
        haber una secuencia en el fichero.

        Input: Fichero FASTA

        Output: objeto Fasta con la header y secuencia
        """
        #capturar la linea header
        file = open(file,"r")  # archivo FASTA
        header = file.readline() 
        header = header.strip("\n")
        # capturar la secuencia
        seq = []
        for line in file:
            line = line.strip("\n")
            seq.append(line)
        # convertir la secuencia a string
        seq = "".join(seq)
        fasta_obj = Fasta(header,seq)     # objeto tipo Fasta del fichero
        return(fasta_obj)             

    
    def graba_fa(self,nombre_out):
        """
        Esta función se encarga de recibir una header y una secuencia
        y obtener un fichero FASTA escrito.

        Inputs:
            - nombre_out (str) : El nombre del fichero de salida
        """
        # abrir fichero de salida en formato escritura
        nombre_out = open(nombre_out + ".fa", "w")
        nombre_out.write(self.header + "\n")

        try:
            delimitador = int(input("Cuantos caracteres por linea de \
secuencia desea: "))             # Delimitador de carácteres por linea

        except:
            print("Seleccione un entero como delimitador")
        
        else:
            secuencia = self.seq
            for i in range(0,len(secuencia),delimitador):
                linea = secuencia[i: i + delimitador]
                nombre_out.write(linea + "\n")
            
            # cerrar el fichero
            nombre_out.close()


    def virus(self,motivo) -> list:
        """
        Esta función determina las posiciones en las que se encuentra un 
        determinado patrón del virus en la secuencia.

        Input:
            - Motivo (str)
        Output:
            - Lista de posiciones de match. 
        """
        secuencia = self.seq
        indices = re.finditer(motivo,secuencia)
        # la variable indices almacena toda la informacion de cada uno de los 
        # matches del método finditer.

        posiciones = []
        for indice in indices:
            posicion = indice.start()
            posiciones.append(posicion)

        if posiciones == []:
            return("No se ha encontrado el motivo")
        else:
            return(posiciones)


    def tripletesint(self, bases = ["A","T","G","C"]) -> str:
        """
        La función sec2int convierte la secuencia de DNA a secuencia de 
        números. Hay un número por cada triplete.

        Input: 
            - Lista de bases
        Output:
            - Str con la secuencia numérica
        """
        secuencia = self.seq
        
        dic_ints = {}         
        numero_tr = 0
        for base1 in bases:
            for base2 in bases:
                for base3 in bases:
                    numero_tr += 1
                    triplete = base1 + base2 + base3
                    # se le asigna un número a cada codon
                    dic_ints[triplete] = numero_tr
        
        seq_int = []
        for i in range(0,len(secuencia),3):
            triplete = secuencia[i:i + 3]
            if (len(triplete) != 3):
                continue
            elif triplete not in dic_ints.keys():
                seq_int.append("N")   #en caso de que haya caracteres extraños se añade una N.
                continue    

            num = dic_ints[triplete]
            seq_int.append(str(num))
        
        # dado que hay codones de más de 1 digito, se separaran los codigos 
        # por ;
        seq_int =  ";".join(seq_int)
        return(seq_int)


    def calcular_GC(self) -> float:
        """
        Esta función es capaz de extraer la secuencia del objeto FASTA
        y calcular su porcentaje GC.

        Inputs:
            - Objeto Fasta
        Outputs
            - Porcentaje GC (float)
        """
        secuencia = self.seq
        numG,numC = secuencia.count("G"),secuencia.count("C")
        longitud = len(secuencia)
        perc_GC = (numG + numC)*100/longitud
        perc_GC = round(perc_GC,3)
        return(perc_GC)


    def buscar_genes(self) -> list:
        """
        Esta función busca genes en una secuencia mediante una 
        expresión regular. 

        Inputs:
            - Objeto Fasta
        Outputs:
            - Secuencia de los genes en formato lista
        """
        secuencia = self.seq
        motivo = "(ATG)([ATGC]{3})*?(TGA|TAG|TAA)"
        indices = re.finditer(motivo,secuencia)
        lon_signific = int(input("Longitud mínima del gen: "))
        num_gen = 0
        for indice in indices:
            start,stop = indice.span()
            gen = secuencia[start:stop]
            frame = start % 3
            # solo mostrar genes con más de x nucleotidos
            if len(gen) >= lon_signific:
                num_gen += 1
                print("@gen%d_frame%d" %(num_gen,frame + 1))
                print(gen + "\n")


    def mutar_secuencia(self) -> str:
        """
        La función mutar_secuencia permite realizar sustituciones de bases, 
        inserciones y delecciones sobre la secuencia problema. Es recomendable
        concoer muy bien la posición en la que se desean realizar los cambios.

        Input:
            - Objeto FASTA
        Outputs:
            - Secuencia FASTA mutada (str)
        """
        secuencia = self.seq
        secuencia = list(secuencia)

        num_cambios = int(input("Cuantos cambios quiere: "))
        for i in range(num_cambios):
            print("Qué tipo de mutación quiere para el cambio numero %d:\n " %(i + 1))
            print("1. Sustitución: Cambia una base por otra")
            print("2. Inserción: Añade 1 o más bases")
            print("3. Delección: Elimina 1 o más bases de la secuencia")

            decision = input("Escoja su opción: ")
            posibilidades = ["1","2","3"]

            if decision not in posibilidades:
                print("\nNo está en las opciones")
                exit()
            posicion = int(input("Qué posición quiere atacar: ")) - 1

            if decision == "1":
                secuencia[posicion] = input("Nueva base: ")
            elif decision == "2":
                posicion += 1
                secuencia.insert(posicion,input("Nueva base o bases: "))
            elif decision == "3":
                del secuencia[posicion: posicion + int(input("Cuantas posiciones eliminar: "))]


        secuencia = "".join(secuencia)
        return(secuencia)
    
    #Para aumentar el rango de posibilidades voy a importar a la clase
    #algunas funciones interesantes del boletín Funciones, más concretamente
    # del fichero funciones_extra.py.
    from funciones_extra import cuenta_nt_tipo
    from funciones_extra import cuenta_codones
    from funciones_extra import contador_mutaciones
    from funciones_extra import traduccion_seq
    from funciones_extra import reversa_complementaria


    def dot_plot(self):
        """
        Función que representa un dotplot de DNA. La función pedirá por teclado
        la secuencia con la que comparar al ejecutar la función.

        Prerequisitos: La secuencia con la que comparar ha de ser de la misma
        longitud que la secuencia problema.

        Inputs:
            - Objeto fasta
            - Secuencia a comparar
        Outputs:
            - dotplot
        """
        secuencia1 = self.seq
        secuencia2 = input("Secuencia de %d bases: " %(len(secuencia1)))

        if len(secuencia1) != len(secuencia2):
            print("No se puede mostrar el dotplot por diferencias de longitud")
        else:
            secuencia1, secuencia2 = np.array(list(secuencia1)), np.array(list(secuencia2))
            plt.imshow(secuencia1==secuencia2[:,None])        
            plt.xticks(np.arange(len(secuencia1)), secuencia1)
            plt.yticks(np.arange(len(secuencia2)), secuencia2)
            plt.show()    



# Programa principal/Menu de opciones

def __main__():
    fasta_obj = Fasta.constructor(fasta)
    secuencia = fasta_obj.seq

    if Fasta.evaluar_secuencia(fasta_obj) is False:
        print("El fasta no contiene una secuencia de DNA ") 
        return()

    lista_funciones = ["evaluar_secuencia", "graba_fa","virus","tripletesint"
,"calcular_GC","buscar_genes","mutar_secuencia","dot_plot" ,"cuenta_nt_tipo"
,"cuenta_codones", "contador_mutaciones", "reversa_complementaria",
"traduccion_seq"]
    
    print("\nFunciones posibles:\n")
    for j,funcion in enumerate(lista_funciones):
        opcion = j + 1
        print(str(opcion) + ". " + funcion)
    
    try:
        opcion = int(input("\nQué función quiere?: "))
    except ValueError:
        print("Eso no es siquiera un número")
        return()
    
    rango_funciones = list(range(1,len(lista_funciones)+1))
    if opcion not in rango_funciones:
        print("Esta función no está disponible")

    elif opcion == 1:
        resultado = Fasta.evaluar_secuencia(fasta_obj)
        if resultado is True:
            print("Es una secuencia de ADN")
    elif opcion == 2:
        nom_salida = input("Nombre del fichero de salida: ")
        Fasta.graba_fa(fasta_obj,nom_salida)
        
    elif opcion == 3:
        motivo = "(AG)+(C)?AGATA([ACGT])+?(GAT){2,3}"
        print(Fasta.virus(fasta_obj,motivo))
        
    elif opcion == 4:
        print(Fasta.tripletesint(fasta_obj, bases = ["A","T","G","C"] ))
        
    elif opcion == 5:
        print(Fasta.calcular_GC(fasta_obj))
        
    elif opcion == 6:
        Fasta.buscar_genes(fasta_obj)
        
    elif opcion == 7:
        print(Fasta.mutar_secuencia(fasta_obj))
        
    elif opcion == 8:
        Fasta.dot_plot(fasta_obj)
        
    elif opcion == 9:
        print(Fasta.cuenta_nt_tipo(secuencia))
        
    elif opcion == 10:
        frame = int(input("Qué frame desea(1,2,3): "))
        print(Fasta.cuenta_codones(secuencia,frame))
        
    elif opcion == 11:
        seq2 = input("Introduzca la secuencia con la que comparar: ")
        print(Fasta.contador_mutaciones(secuencia,seq2))
        
    elif opcion == 12:
        print(Fasta.reversa_complementaria(secuencia))

    elif opcion == 13:
        frame = int(input("Elija el frame de lectura (1,2,3): "))
        if frame in [1,2,3]:
            secuencia = secuencia[frame - 1:]
            print(Fasta.traduccion_seq(secuencia))
        else:
            print("\nEse frame no es coherente")

if __name__ == __main__():
    __main__()



