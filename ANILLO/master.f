c Programa master - practica 1 programación paralela
c Alumno Ivan Gomez
      program master
      implicit none
      include 'fpvm3.h'

c Declaración de variables
      integer nhost, numt, tids, info, var, mytid

c Inicializa PVM
	  call pvmfmytid(mytid)

c Lee la variable de entrada
      print *,'Introduce un número entero:'
      read *, var

c Expande los procesos hijos
      call expande('slave', tids, nhost, numt)

c Comprobamos que tenemos más se ha distribuido a más de un host
	  if (nhost.lt.0) then

c       Inicializa el buffer para el intercambio y envía la información
c         Envia:
c   	  - tid del proceso padre (al que hay que devolver al valor)
c	      - tids de los procesos disponibles (esclavos)
c         - número de tids
c  	      - valor del input a incrementar
	  	call pvmfinitsend(PVMDATADRAW, 1)
	  	call pvmfpack(INTEGER4, mytid, 1, 1, info)
		call pvmfpack(INTEGER4, nhost, 1, 1, info)
		call pvmfpack(INTEGER4, tids, nhost, 1, info)
		call pvmfpack(INTEGER4, var, 1, 1, info)
		call pvmfsend(tids(i), 1, info)
	  endif
c Espera el resultado final del último exclavo
	  
	  call pvmfrecv(-1, 3, info)
      call pvmfunpack(INTEGER4, var, 1, 1, info)

	  print*, 'La suma del programa es: ', var

c Sale de pvm y exit del programa
	  call pvmfexit( info)
      print*, 'He salido de PVM; codigo de estado: ',info
      end
