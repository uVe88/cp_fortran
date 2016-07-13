  	  program jacobi 
      include 'fpvm3.h'
c ---------------------------------------------------------
c Ejemplo de programa fortran para ilustrar el uso de PVM 3
c ---------------------------------------------------------
      
      parameter( maxm=100, maxntids=8)
      integer i, info, nproc, nhost, msgtype
      integer mytid, iptid,numt,filasproc,resto,numfilas,gnum
	  real*8 A(maxm,maxm)
	  real*8 b(maxm)
	  real*8 diag(maxm)
	  integer m,n, filaini, filafin, k
	  integer tids(maxntids)
	  real*8 bk,diagk
	  real*8 sum
	  real*8 xsol(maxm),x(maxm)

	  
c     Incluimos este programa en PVM
c     ----------------------------------------	  
	  call pvmfmytid( mytid )
c     preguntamos por el proceso padre
	  call pvmfparent( iptid )
	  call pvmfjoingroup('test',gnum)
	  print *, 'Mi tid es ', mytid,' y mi numero en el grupo es ',gnum, 'mi padre es', iptid
	  CALL pvmfcatchout( 1, info )
c 	  Inicializamos resto para los procesos hijos
	  resto = 0
      print *,'Tengo padre:', iptid, pvmnoparent


      if (iptid.eq.pvmnoparent) then
c	  Si no tiene padre hacemos el proceso principal 
c	  Leemos los inputs
	    print *,'Numero de filas'
        read *, n		
 	    call construirDatos(n,A,b,diag)
		call pvmfspawn( 'jacobi', PVMDEFAULT, '*', maxntids, tids, numt)
c	TODO. Machacar todo si hay un numero mayor de procesos que maxntids
		filasproc = int(n/numt)
		resto = n -(filasproc*numt)
c     inicializamos al resultado
		 filaini = 1
		 print *,'Comenzando envio de info',numt  
c	  mandamos los dato iniciales a los procesos
	    do i=1,numt
	       filafin = filaini+filasproc  
c		   Mandamos filaini y filafin
           print *,'Comenzando envio de info',numt 
	       call pvmfinitsend( PVMDATARAW, info)
           call pvmfpack( INTEGER4, filaini, 1, 1, info)
		   call pvmfpack( INTEGER4, filafin, 1, 1, info)
c          El numero total de datos es necesario para recorrer las matrices
		   call pvmfpack(INTEGER4,n,1,1,info)
c          Mandamos las filas de a
           do j=filaini,filafin
				call pvmfpack (REAL8,A(filaini,1),n,maxm,info)
				call pvmfpack (REAL8,b(filaini),1,1,info)
				call pvmfpack (REAL8,diag(filaini),1,1,info)
		   enddo
           call pvmfsend( tids(i), 2, info)
		   print *,'Comenzando envio de info',numt 
	       filaini = filafin + 1;
	    enddo
		filasproc = resto

      else
	  	print *,'iniciando recepcion de datos'
c 	    Recibimos los datos
		call pvmfrecv(iptid,2,info)
		call pvmfunpack(INTEGER4,filaini,1,1,info)
		call pvmfunpack(INTEGER4,filafin,1,1,info)
		filasproc = filafin - filaini
		call pvmfpack(INTEGER4,n,1,1,info)
		do j=filaini,filafin
				call pvmfunpack (REAL8,A(filaini,1),n,maxm,info)
				call pvmfunpack (REAL8,b(filaini),1,1,info)
				call pvmfunpack (REAL8,diag(filaini),1,1,info)
		enddo
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
      do while ( fin .eq. 1)
c   Comprobacion si hay resto en el proceso principal
		  if (filasproc .gt. 0) then
			do i=1,filasproc
			  k = filaini+i-1
			  diagk=diag(k)
			  xsol(i) = call calcularJacobi(n,A,b(k),diagk,x(k),k)
			  sum = sum + abs((x(i)- xsol(i)))
			enddo
		  endif

c   Llamamos a pvmreduce para calcular la norma y ver la salida, TODO sustituir pvmSum por norma
		  call pvmfreduce(PvmSum,sum,1,REAL8,2,'jacobi',0,info)
		  print *,'suma',sum
c 	Necesario llamar a pvmfbarrier ya que pvmfreduce no es bloqueate
          print *,'Final de todas las llamadas'
		  call pvmfbarrier('jacobi',ntask+1,info)
		  if (gnum .eq. 0) then
		  	 fin = norma(sum)  
		  endif
		  
	  enddo
	  
	  call pvmflvgroup('jacobi',info)
	  call pvmfexit(info)
	  
      end


	  
c     subrutina para calcular jacobi a partir de un determinado numero de valores
	  real*8 function  calcularJacobi(n,A,bk,diag,x,k)
	  implicit none
	  integer n    
	  real*8 A(n,n)
	  real*8 bk,diag,sum,x
	  integer i,k
	  sum = 0
	  do i=1,n
          sum = sum + (A(k,i)*x)
	  enddo
	  calcularJacobi = (bk-sum)/diag
	  return
	  end
	  
	
c     */Construimos la funcion principal en función de los parametros dados	 /* 
	  subroutine construirDatos (n,A,b,diag)
	  implicit none
	  real*8 A(n,n)
	  real*8 b(n)
	  real*8 diag(n)
	  integer n,i,j
	  
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

	  