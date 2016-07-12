c Programa master - practica 1 programación paralela
c Alumno Ivan Gomez
       	program master
       	implicit none
      	include 'fpvm3.h'
c Declaración de variables
       	integer nhost, numt, tids, info, var

c Lee la variable de entrada

c Inicializa PVM

c Expande los procesos hijos
        call expande('slave', tids, nhost, numt)

c Inicializa el buffer para el intercambio y envía la información
c Envia:
c   	- tid del proceso padre
c	- tids de los procesos disponibles (esclavos)
c  	- valor del input a incrementar
c	- 
	call pvmfinitsend(PVMDATADRAW, 1)
	call pvmfpack(INTEGER2, var, 1, 1, info)

	
c Espera el resultado final del último exclavo

	call pvmfmcast(nhost, tids, 1, info)



c Sale de pvm y exit del programa
	call pvmfexit( info)
        print*, 'He salido de PVM; codigo de estado: ',info
       	end
