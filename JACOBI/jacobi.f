  	  program jacobi 
      include 'fpvm3.h'
c ---------------------------------------------------------
c Ejemplo de programa fortran para ilustrar el uso de PVM 3
c ---------------------------------------------------------
      
      parameter( maxm=100, maxntids=2)
      integer i, info, nproc, nhost, msgtype
      integer mytid, iptid,numt,filasproc,resto,numfilas,gnum,root
	  real*8 A(maxm,maxm)
	  real*8 b(maxm)
	  real*8 btemp(maxm)
	  integer m,n, filaini, filafin, filapadre,filaact, k,z,fin
	  integer tids(maxntids)
	  real*8 sum,sumx
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
 	    call construirDatos(n,A,b,d,maxm)
c  TODO. Hacer 	que si n < mxntids -> maxntids = n	 
c	TODO. Machacar todo si hay un numero mayor de procesos que maxntids
        call pvmfspawn('jacobi',PVMTASKDEFAULT,'*',maxntids,tids,numt)
        
		filasproc = n/numt
		resto = n -(filasproc*numt)
c     inicializamos al resultado
		 filaini = 1
c	  mandamos los dato iniciales a los procesos
        call escribe(A,maxm,n,n)	
	    do i=1,numt
	       filafin = filaini+filasproc -1 
c		   Mandamos filaini y filafin
		  
           
	       call pvmfinitsend( PVMDATARAW, info)
		   call pvmfpack(INTEGER4,filaini,1,1,info)
           call pvmfpack( INTEGER4, filasproc, 1, 1, info)
		   call pvmfpack( INTEGER4,n,1,1,info)
		   
c Enviamos el numero de tareas para hacer los barriers
		   call pvmfpack( INTEGER4,numt,1,1,info)
c  Enviamos el gnum del parent pava el pvmfreduce
           root = gnum  
		   
		   call pvmfpack(INTEGER4,root,1,1,info)
		
		
c          Mandamos las filas de a
           do j=filaini,filafin
				call pvmfpack (REAL8,A(j,1),n,maxm,info)				
		   enddo
		  
		   
		   k=1
		   do j=filaini,filafin
		      btemp(k) = b(j)
              k=k+1	
		   enddo
		   call pvmfpack(REAL8,btemp,filasproc,1,info)
		   filaini=filafin+1
		   filafin=filafin+resto-1
           call pvmfsend( tids(i), 2, info)
		 
c		Inicializamos el resultado a cero
	      do j=1,n
		       x(j) = 0
	      enddo
	     
	    enddo
        filafin=filafin+resto
		filapadre = 1

      else
c     Si tiene padre hacemos el proceso de los hijos  
	  	print *,'iniciando recepcion de datos'
c 	    Recibimos los datos
		call pvmfrecv(iptid,2,info)

        call pvmfunpack(INTEGER4,filapadre,1,1,info)  
		print *,'filapadre',filapadre   
        call pvmfunpack(INTEGER4, filasproc, 1, 1, info)
		call pvmfunpack(INTEGER4,n,1,1,info)
		call pvmfunpack( INTEGER4,numt,1,1,info)
		call pvmfunpack(INTEGER4,root,1,1,info)
		
		do j=1,filasproc     
			call pvmfunpack (REAL8,A(j,1),n,maxm,info)
		enddo
		
        call pvmfunpack (REAL8,b,filasproc,1,info)

c		call escribe(A,maxm,filafin,n)
		print *,'datos recibidos'
	    filaini=1
		filafin=filasproc
       

	  endif 

      call pvmfbarrier('test',numt+1,info)
c 	Hacemos un scatter de b y de la diagonal


      

	   fin = 0
      do while ( fin .eq. 0)
c Hacemos un broadcast de toda la solucion x
      print *,'Entra en bucle'
	  print *,gnum,root
      if (gnum .eq. root) then
	       call pvmfinitsend(PVMDATARAW,info)
		   call pvmfpack(REAL8,x,n,1,info)
		   call pvmfbcast('test',1,info) 
           print *,'xpadre:',x(1),x(2),x(3),x(4)		
	  else
	 	   call pvmfrecv(iptid,1,info)
		   call pvmfunpack(REAL8,x,n,1,info)
		   print *,'xhijo:',x(1),x(2),x(3),x(4)			  
	  endif
     
c Calculamos la solucion parcial con los datos
      	
      print *,'filaini: ',filaini,'filafin: ',filafin
	  print *,'filapadre: ',filapadre
	  print *,filasproc
c	  call escribe(A,maxm,filaini,filafin) 
      do i=filaini,filafin
	    sumx=0
		filaact = filaini+filapadre-1
		do j=1,n
		  if (j .ne. filaact) then
		     sumx = sumx+(A(i,j)*x(j))
			 print *,'parcial j',j,sumx,A(i,j),x(j)
		  endif
		enddo
		print *,'parcial i',i,b(i),sumx,A(i,i)
		x(i) = (b(i)-sumx/A(i,filaact))
		print *,'x parcial',filapadre,x(i)
	  enddo	
c Recogemos las soluciones parciales en un vector

	   
       call pvmfgather(x,x,filasproc,REAL8,2,'test',root,info)

c Realizamos el calculo de Jacobi para las filas que quedan 
	   if (gnum .eq. root) then
c   Calculo la norma
         fin = norma(xsol,x,maxm,n)
         print *,'Calculando norma'
		 print *,'resultado',x(1),x(2),x(3),x(4)

	

        
	   endif

      enddo
	  			

c Matamos todos los procesos hijos		  
	  call shutdown(numt, tids)
	  print *,'Salgo'
	  call pvmflvgroup('test',info)
	  call pvmfexit(info)
	  
      end

      subroutine shutdown (numt,tids)
      integer numt,i, tids(numt)

	  do i=1,numt
         call pvmfkill(tids(i),info) 
	  enddo
      return   
	  end 
	  
	
c     */Construimos la funcion principal en función de los parametros dados	 /* 
	  subroutine construirDatos (n,A,b,d,maxm)
	  implicit none
	  real*8 A(maxm,maxm)
	  real*8 b(maxm)
	  real*8 d(maxm)
	  integer n,i,j,maxm
	  
	  do j=1,n
        do i=1,n
		  if ( i .eq. j ) then
			A(i,j)=n
			d(i)=n
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
      
	

	  