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
	  real*8 d(maxm)
	  integer m,n, filaini, filafin, k,z,fin
	  integer tids(maxntids)
	  real*8 sum,sumx
	  real*8 xsol(maxm),x(maxm)

	  
c     Incluimos este programa en PVM
c     ----------------------------------------	  
	  call pvmfmytid( mytid )
c     preguntamos por el proceso padre
	  call pvmfparent( iptid )
	  call pvmfjoingroup('group1',gnum)
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
           call pvmfpack( INTEGER4, filasproc, 1, 1, info)
		   call pvmfpack( INTEGER4,n,1,1,info)
c Enviamos el numero de tareas para hacer los barriers
		   call pvmfpack( INTEGER4,numt,1,1,info)
c  Enviamos el gnum del parent pava el pvmfreduce
           root = gnum  
		   
		   call pvmfpack(INTEGER4,root,1,1,info)
		
		   print *,'filaini: ',filaini,'filafin: ',filafin
c          Mandamos las filas de a
           do j=filaini,filafin
				call pvmfpack (REAL8,A(j,1),n,maxm,info)				
		   enddo
		   filaini=filafin+1
		   filafin=filafin+resto-1
           call pvmfsend( tids(i), 2, info)
		   print *,'Comenzando envio de info',numt
c		Inicializamos el resultado a cero
	      do j=1,n
		       x(j) = 0
	      enddo
		  print *,'datos padre',b(1),d(1),d(2)  
	    enddo


      else
c     Si tiene padre hacemos el proceso de los hijos  
	  	print *,'iniciando recepcion de datos'
c 	    Recibimos los datos
		call pvmfrecv(iptid,2,info)
c       Recibimos
        call pvmfunpack(INTEGER4, filasproc, 1, 1, info)
		call pvmfunpack(INTEGER4,n,1,1,info)
		call pvmfunpack( INTEGER4,numt,1,1,info)
		call pvmfunpack(INTEGER4,root,1,1,info)
		
		do j=1,filasproc     
			call pvmfunpack (REAL8,A(j,1),n,maxm,info)
		enddo

		call escribe(A,maxm,filafin,n)
		print *,'datos recibidos'
	    filaini=1
		filafin=filasproc
        print *,'filaini: ',filaini,'filafin: ',filafin
		print *,'n: ',n,'.maxm: ',maxm

	  endif 
      call pvmfbarrier('group1',numt+1,info)
c 	Hacemos un scatter de b y de la diagonal
      print *,filasproc
	  call pvmfscatter(b,b,filasproc,REAL8,1,'group1',root,info)
	  call pvmfscatter(d,d,filasproc,REAL8,2,'group1',root,info)
	  print *,'datos inicializados',b(1),d(1),d(2)  
      

	  fin = 0
c      do while ( fin .eq. 0)
	   do z=1,3
c Hacemos un broadcast de toda la solucion x
       print *,'Entra en bucle'
      if (gnum .eq. gnump) then
           call pvmfrecv(iptid,1,info)
		   call pvmfunpack(REAL8,x,n,maxm,info)
	  else
	 	   call pvmfinitsend(PVMDATARAW,info)
		   call pvmfpack(REAL8,x,n,maxm,info)
		   call pvmfbcast('group1',1,info)  		
	  endif
     
c Calculamos la solucion parcial con los datos
      	
      print *,'filaini: ',filaini,'filafin: ',filafin
	  call escribe(A,maxm,filaini,filafin) 
      do i=filaini,filafin
	    sumx=0
		do j=1,n
		  if (j .ne. i) then
		     sumx = sumx+A(i,j)*x(j)
		  endif
		enddo
		print *,'sumx : ',sumx,'di:',d(i)
		x(i) = (b(i)-sumx/d(i))
	  enddo	
c Recogemos las soluciones parciales en un vector
       print *,'antes del gather',filasproc,root
	   
       call pvmfgather(x,x,filasproc,REAL8,2,'group1',root,info)
	   print *,'en el grupo soy: ',gnum
       print *,('x(',i,')=',x(i),'; ',i=1,n,filasproc)
c Realizamos el calculo de Jacobi para las filas que quedan 
	   if (gnum .eq. root) then
c   Calculo la norma
            	

	   endif
c      Sincronizo procesos	  


      enddo
	  			
	 	
            
		  
		  
c 	Necesario llamar a pvmfbarrier ya que pvmfreduce no es bloqueante
      print *,'Final de todas las llamadas de: ',gnum, gnump

		  
	 
	  
	  call pvmflvgroup('group1',info)
	  call pvmfexit(info)
	  
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
      
	

	  