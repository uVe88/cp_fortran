c Programa master - practica 1 programación paralela
c Alumno Ivan Gomez
      program master
      implicit none
      include 'fpvm3.h'

c Declaración de variables
	integer maxtids
      parameter(maxtids=8)
      integer nhost, numt, tids(maxtids), info, var, mytid, i

c Inicializa PVM
	call pvmfmytid(mytid)
	call pvmfcatchout(1,info)

	print*, 'Maestro tid:', mytid

c Lee la variable de entrada
      print *,'Introduce un número entero:'
      read *, var
	print *,'Introduce un número procesos:'
      read *, nhost
	print *,'Leido'
c Expande los procesos hijos
      call expande('slave', tids, nhost, numt)

c	call pvmfspawn('slave',PVMTASKDEFAULT,'*',nhost,tids,numt)
	print *,'Propagado a:', nhost, numt, tids

c Comprobamos que se ha distribuido a más de un host
	if (numt.gt.0) then

c       Inicializa el buffer para el intercambio y envía la información
c         Envia:
c   	  - tid del proceso padre (al que hay que devolver al valor)
c	      - tids de los procesos disponibles (esclavos)
c         - número de tids
c  	      - valor del input a incrementar
	  print*,'Maestro va a enviar a esclavo con tid: ',tids(1)
	  print*,'Info: ',mytid, numt, tids, var
	  call pvmfinitsend(PVMDATARAW, info)
	  call pvmfpack(INTEGER4, mytid, 1, 1, info)
	  call pvmfpack(INTEGER4, numt, 1, 1, info)
	  call pvmfpack(INTEGER4, tids, numt, 1, info)
	  call pvmfpack(INTEGER4, var, 1, 1, info)
	  call pvmfsend(tids(1), 1, info)
	  print*,'Maestro envió a esclavo con tid: ',tids(1)
	
c Espera el resultado final del último exclavo
	  print*, 'Espero los datos del ultimo nodo del anillo'
	  call pvmfrecv(-1, 3, info)
        call pvmfunpack(INTEGER4, var, 1, 1, info)
	  print*, 'Datos recibidos'
	  print*, 'La suma del programa es: ', var

c Sale de pvm y exit del programa
	  call pvmfexit( info)
        print*, 'He salido de PVM; codigo de estado: ',info
	endif
      end
