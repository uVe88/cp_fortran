c Programa master - practica 1 programación paralela
c Alumno Ivan Gomez
      program master
      implicit none
      include 'fpvm3.h'

c Declaración de variables
	integer maxtids, clusterHost
      parameter(maxtids=100)
      integer nhost, numt, tids(maxtids), info, var, mytid, i

	clusterHost=1

c Inicializa PVM
	call pvmfmytid(mytid)
	call pvmfcatchout(1,info)

	print*, 'Maestro tid:', mytid

c Lee la variable de entrada
      print *,'Introduce un número entero:'
      read *, var

c Expande los procesos hijos

	if (clusterHost.eq.0) then
	  print *,'Introduce un número procesos:'
        read *, nhost
	  call pvmfspawn('slave',PVMTASKDEFAULT,'*',nhost,tids,numt)
	  print *,'Propagado a:', nhost, numt, (tids(i), i=1, numt)
	else
        call expande('slave', tids, nhost, numt)
	  print *,'Propagado a:', nhost, numt, (tids(i), i=1, numt)
	  
	  if (nhost-1.le.0) then
	  	print*,'Warning: No hay hosts'
	  	call pvmfexit(info)
		stop
	  endif
 	endif

c Comprobamos que se ha distribuido a más de un host
	if (numt.gt.0) then

c       Inicializa el buffer para el intercambio y envía la información
c         Envia:
c   	  - tid del proceso padre (al que hay que devolver al valor)
c	      - tids de los procesos disponibles (esclavos)
c         - número de tids
c  	      - valor del input a incrementar
	  print*,'Maestro va a enviar a esclavo con tid: ',tids(1)
	  print*,'Info: ',mytid, numt, var, (tids(i), i=1, numt)
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
	else
	  print*, 'Error, no pids. Value:', numt
	endif
      end
