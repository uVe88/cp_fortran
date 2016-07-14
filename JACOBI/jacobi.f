  	  program jacobi 
      include 'fpvm3.h'
c ---------------------------------------------------------
c Ejemplo de programa fortran para ilustrar el uso de PVM 3
c ---------------------------------------------------------
      
      parameter( maxm=100, maxntids=2)
      integer i, info, nproc, nhost, msgtype
      integer mytid, iptid,numt,filasproc,resto,numfilas,gnum,gnump
	  real*8 A(maxm,maxm)
	  real*8 b(maxm)
	  real*8 diag(maxm)
	  integer m,n, filaini, filafin, k,z
	  integer tids(maxntids)
	  real*8 bk,diagk
	  real*8 sum,sumx
	  real*8 xsol(maxm),x(maxm)

	  
c     Incluimos este programa en PVM
c     ----------------------------------------	  
	  call pvmfmytid( mytid )
c     preguntamos por el proceso padre
	  call pvmfparent( iptid )
	  call pvmfjoingroup('group1',gnum)
	  print *, 'Mi tid es ', mytid,' y mi numero en el grupo es ',gnum, 'mi padre es', iptid
c	  CALL pvmfcatchout( 1, info )
c 	  Inicializamos resto para los procesos hijos
	  resto = 0
      print *,'Tengo padre:', iptid, pvmnoparent


      if (iptid.eq.pvmnoparent) then
c	  Si no tiene padre hacemos el proceso principal 
c	  Leemos los inputs
	    print *,'Numero de filas'
        read *, n		
 	    call construirDatos(n,A,b,diag,maxm)
c  TODO. Hacer 	que si n < mxntids -> maxntids = n	 
c	TODO. Machacar todo si hay un numero mayor de procesos que maxntids
        call pvmfspawn('jacobi',PVMTASKDEFAULT,'*',maxntids,tids,numt)
        
		filasproc = n/numt
		resto = n -(filasproc*numt)
		print *,'filasproc: ',filasproc
c     inicializamos al resultado
		 filaini = 1
		 print *,'Comenzando envio de info',numt  
c	  mandamos los dato iniciales a los procesos
	    do i=1,numt
	       filafin = filaini+filasproc -1 
c		   Mandamos filaini y filafin
		   call escribe(A,maxm,n,n)	
           print *,'Comenzando envio de info',numt 
	       call pvmfinitsend( PVMDATARAW, info)
           call pvmfpack( INTEGER4, filaini, 1, 1, info)
		   call pvmfpack( INTEGER4, filafin, 1, 1, info)
		   call pvmfpack( INTEGER4,n,1,1,info)
c Enviamos el numero de tareas para hacer los barriers
		   call pvmfpack( INTEGER4,numt,1,1,info)
c  Enviamos el gnum del parent pava el pvmfreduce
		   call pvmfpack(INTEGER4,gnum,1,1,info)
		   gnump = gnum
		   print *,'filaini: ',filaini,'filafin: ',filafin
c          Mandamos las filas de a
           do j=filaini,filafin
				call pvmfpack (REAL8,A(j,1),n,maxm,info)
c			
c				
		   enddo
c	Enviamos el vector b
		   do j=filaini,filafin
                call pvmfpack (REAL8,b(j),1,1,info)
		   enddo
c   Enviamos el vector con las diagonales
		  do j=filaini,filafin
			  call pvmfpack (REAL8,diag(j),1,1,info)
		  enddo

           call pvmfsend( tids(i), 2, info)
		   print *,'Comenzando envio de info',numt
c		Esto lo hacemos para empezar a trabajar con el resto 
	       filaini = filafin + 1;
		   filafin = filaini + filasproc -1
	    enddo
		filasproc = resto

      else
	  	print *,'iniciando recepcion de datos'
c 	    Recibimos los datos
		call pvmfrecv(iptid,2,info)
		call pvmfunpack(INTEGER4,filaini,1,1,info)
		call pvmfunpack(INTEGER4,filafin,1,1,info)
		call pvmfunpack(INTEGER4,n,1,1,info)
		call pvmfunpack( INTEGER4,numt,1,1,info)
		call pvmfunpack(INTEGER4,gnump,1,1,info)
		filasproc = filafin - filaini+1
		print *,'filaini: ',filaini,'filafin: ',filafin
		print *,'n: ',n,'.maxm: ',maxm
		do j=filaini,filafin     
			call pvmfunpack (REAL8,A(j,1),n,maxm,info)
		enddo
c	Recibimos el vector b
		do j=filaini,filafin
            call pvmfunpack (REAL8,b(j),1,1,info)
		enddo
c   Recibimos el vector con las diagonales
		  do j=filaini,filafin
			  call pvmfunpack (REAL8,diag(j),1,1,info)
		  enddo
		call escribe(A,maxm,filafin,n)
		print *,'datos recibidos'
	   
c     Si tiene padre hacemos el proceso de los hijos  

	  endif 

c	Inicializamos las x a cero como ponia el enunciado, si se quisiera leer de un fichero, habria que mandarlo por send
	  do i=filaini,filafin
		  x(i) = 0
	  enddo
	  print *,'x inicializadas'
	  fin=0
	  k=1
	  print *,'fin: ',fin
c      do while ( fin .eq. 0)
      do z=1,4
c   Comprobacion si hay resto en el proceso principal
      print *,'Entrando en bucle'
	  print *,'filasproc:',filasproc
		  if (filasproc .gt. 0) then
c si hay mas de una fila hacemos el calculo de jacobi para las filas
			do i=0,filasproc-1
			  k = filaini+i
              sumx = 0
	  			do j=1,n
          			sumx = sumx  + (A(k,j)*x(k))
	  			enddo
	 print *,'k: ',k,'bk: ',b(k),'sum x: ',sumx, 'diagk: ',diag(k)	  
			  xsol(i) = (b(k)-sumx)/diag(k)
			  print *,'x ',filaini,'= ',xsol(i)
			  sum = sum + abs((x(i)- xsol(i)))
			enddo
		  endif
          print *,'sum:', sum,gnump
          print *,'ntask: ',numt
c   Llamamos a pvmreduce para calcular la norma y ver la salida, TODO sustituir pvmSum por norma
		  call pvmfreduce(PvmSum,sum,1,REAL8,2,'group1',gnump,info)
		 print *,'barrera'
c		  call pvmfbarrier('group1',numt+1,info)
c 	Necesario llamar a pvmfbarrier ya que pvmfreduce no es bloqueante
          print *,'Final de todas las llamadas de: ',gnum, gnump
		
		 
		  if (gnum .eq. gnump) then
		     print *,'suma',sum
		  	 fin = norma(sum)  
		  endif
		  
	  enddo
	  
	  call pvmflvgroup('group1',info)
	  call pvmfexit(info)
	  
      end


	  

	
c     */Construimos la funcion principal en función de los parametros dados	 /* 
	  subroutine construirDatos (n,A,b,diag,maxm)
	  implicit none
	  real*8 A(maxm,maxm)
	  real*8 b(maxm)
	  real*8 diag(maxm)
	  integer n,i,j,maxm
	  
	  do j=1,n
        do i=1,n
		  if ( i .eq. j ) then
			A(i,j)=n
			diag(i)=n
		  else
			A(i,j)=1
		  endif
        enddo
      enddo
	  do i=1,n
		b(i) = 2*n-1 	
	  enddo
      return
      end	  
      
	

	  