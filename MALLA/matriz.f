      program matriz
      implicit none
      include 'fpvm3.h'
        
      integer ntask, info, mytid, mygid, i, blqtam, row, col
      integer uptid, downtid, dchatid, izqtid
      integer ATAG, BTAG, AATAG, BBTAG
      integer MAXNTIDS, DIMTAG, j, test
      integer MAXR, MAXBLOQTAM, mygtid, nfilas, tgid, ttag

      parameter(MAXNTIDS=16)
      parameter(ATAG=2,BTAG=3,AATAG=4,BBTAG=5,DIMTAG=6)
      parameter(MAXR = 4, MAXBLOQTAM = 10)
	    integer tids(MAXNTIDS-1)
      real*8 a(MAXBLOQTAM*MAXBLOQTAM), b(MAXBLOQTAM*MAXBLOQTAM)
      real*8 c(MAXBLOQTAM*MAXBLOQTAM), cprima(MAXBLOQTAM*MAXBLOQTAM)
	 
c     Obtengo el id de mi tarea
      call pvmfmytid(mytid)
      call pvmfcatchout(1,info)

c     Verificamos que pvm ha iniciado bien.
      if (mytid.lt.0) then
        print* , 'Error, tid erroneo'
        stop
      endif

c     Nos unimos al grupo mmult
      call pvmfjoingroup('mmult', mygid)

      if (mygid.lt.0) then
        print*,'Error al obtener el id de grupo'
        call pvmfexit()
        stop
      endif
        
c       /* if my group id is 0 then I must spawn the other tasks */
      if (mygid.eq.0) then

c       Soy el root del grupo
        print *,'Introduce el valor de nfilas r (entero):'
        read *, nfilas

        print *,'Introduce el tamaño de bloque (entero):'
        read *, blqtam

c       Asigno el numero de procesos de la malla y compruebo que las tareas están en los límites de máximos tids
        ntask = nfilas*nfilas
        if ((ntask.lt.1).or.(ntask.ge.MAXNTIDS)) then
          print*,'ntask = ',ntask,' not valid'
          call pvmflvgroup('mmult',info)
          call pvmfexit(info)
          stop
        endif

c       Comprueba que el tamaña del bloque es menos que el máximo
        if ((blqtam.lt.1).or.(blqtam.ge.MAXBLOQTAM)) then
          print*,'ERROR bloqtam =',blqtam,'not valid;max:',MAXBLOQTAM
          call pvmflvgroup('mmult',info)
          call pvmfexit(info)
          stop
        endif

c       Comprueba que el tamaña del bloque es menos que el máximo
        if ((nfilas.lt.1).or.(nfilas.ge.MAXR)) then
          print*,'ERROR r = ',nfilas,' not valid; max:', MAXR 
          call pvmflvgroup('mmult',info)
          call pvmfexit(info)
          stop
        endif

c       Saltamos la propagación en caso de ser una única tarea.
        if (ntask.eq.1) go to 100

c       Propagamos todas las tareas menos unas que es la mia (root)
        call pvmfspawn('matriz',PVMTASKDEFAULT,'*',ntask-1,tids,info)

c       Comprobamos que se han lanzados todos los procesos de la malla
        if (info.ne.ntask-1) then
          print*,'ERROR!!!', info
          call pvmflvgroup('mmult', info)
          call pvmfexit(info)
          stop
        endif
           
c       Mandadamos la dimension de la matriz a todas las tareas
        
        call pvmfinitsend(PVMDATADEFAULT,info)
        call pvmfpack(INTEGER4, nfilas, 1, 1, info)
        call pvmfpack(INTEGER4, blqtam, 1, 1, info)
        call pvmfmcast(ntask-1, tids, DIMTAG, info)
        
      else
c       Soy miembro de la malla para pero no root
c       Obtengo mi id de grupo

        call pvmfgettid('mmult', 0, mygtid)
        
        
c       Recibo los parámetros de la malla enviados por el root
        call pvmfrecv(mygtid, DIMTAG, info)
        call pvmfunpack(INTEGER4, nfilas, 1, 1, info)
        call pvmfunpack(INTEGER4, blqtam, 1, 1, info)
        
c       Reinicio el número de tareas.       
        ntask = nfilas*nfilas
      endif

c     Sincronizamos con barrier
100   call pvmfbarrier('mmult',ntask, info)

c     Busco mis bloques
      row = mygid/nfilas+1
      col = mod(mygid,nfilas)+1
      

c     Calculo mis vecinos de arriba y abajo 
	    if (row.gt.1) then
        call pvmfgettid('mmult', mygid-nfilas, uptid)
      else
       	call pvmfgettid('mmult', mygid+(nfilas*(nfilas-1)), uptid)
        
      endif
      
      if (row.lt.nfilas) then
        call pvmfgettid('mmult', mygid+nfilas, downtid)
      else
       	call pvmfgettid('mmult', mygid-(nfilas*(nfilas-1)), downtid)   
      endif

c     Calculo mis vecinos de izq y dcha
c     Calculo mi vecinos de izq
	    if (col.gt.1) then
        call pvmfgettid('mmult', mygid-1, izqtid)
      else
       	call pvmfgettid('mmult', mygid+nfilas-1, izqtid) 
      endif
c     Calculo mi vecino dcha
	    if (col.eq.nfilas) then
	   	  call pvmfgettid('mmult', mygid-nfilas+1, dchatid)
	    else
	   	  call pvmfgettid('mmult', mygid+1, dchatid)
      endif

c      print*,'P:',mygid,'F:',row,'C:',col,'tid',mytid
c      print*,'UP:',uptid,'DOWN:',downtid,'IZQ:',izqtid,'DCHA:',dchatid
      
      call pvmfbarrier('mmult',ntask, info)
      
c     Inicializo los bloques
      call iniciarbloques(a, b, c, blqtam, row, col)

c     Calculos la multiplicación de la matriz
      do i=1, nfilas-1
      
c     Alineación inicial y cálculo local y guardar en temporal
        if (row.ne.1) then
c     Los datos del bloque de la matriz a de las filas mayores que 1 debo mandarlo al procesador de mi izq y esperar los mios para computar
          call pvmfinitsend(PVMDATADEFAULT,info)
          call pvmfpack(REAL8, a, blqtam*blqtam, 1, info)
          print*,'A - F:',row,'C:',col,'Envío A a:',izqtid
          call pvmfsend(izqtid, i*ATAG, info)
c     Espero los mios de mi derecha

          call pvmfrecv(dchatid, i*ATAG, info)
          call pvmfunpack(REAL8, a, blqtam*blqtam, 1, info)
        endif

        if (col.ne.1) then
c     Los datos del bloque de la matriz b de las columnas mayores que uno las muevo hacia arriba n veces en funcion de su col n = col-1
          call pvmfinitsend(PVMDATADEFAULT,info)
          call pvmfpack(REAL8, b, blqtam*blqtam, 1, info)
          print*,'B - F:',row,'C:',col,'Envío B a:',uptid, mytid
          call pvmfsend(uptid, i*BTAG, info)
c     Espero los datos del procesador de abajo
          print*,'B - F:',row,'C:',col,'Espero datos de:',dchatid, mytid
          call pvmfrecv(downtid, i*BTAG, info)
          call pvmfunpack(REAL8, b, blqtam*blqtam, 1, info)
          print*,'B - F:',row,'C:',col,'Recibo B de:',dchatid
        endif

        print*, 'Computo AB'
c     Calculamos con los datos ordenados en cada procesador
        call multiplicarbloque(c,a,b,blqtam)

c     Alineación posterior (Movel bloques A y B)

c     Movemos todos los bloques de data A al proc de la izq
        call pvmfinitsend(PVMDATADEFAULT,info)
        call pvmfpack(REAL8, b, blqtam*blqtam, 1, info)
        print*,'F:',row,'C:',col,'Envío a:',izqtid
        call pvmfsend(izqtid, i*AATAG, info)

c     Movemos todos los bloques de datos B al proc de arriba        
        call pvmfinitsend(PVMDATADEFAULT,info)
        call pvmfpack(REAL8, b, blqtam*blqtam, 1, info)
        print*,'F:',row,'C:',col,'Envío a:',uptid
        call pvmfsend(uptid, i*BBTAG, info)
c     Esperamos los bloques de A y B
        call pvmfrecv(dchatid, i*AATAG, info)
        call pvmfunpack(REAL8, a, blqtam*blqtam, 1, info)
        call pvmfrecv(downtid, i*BBTAG, info)
        call pvmfunpack(REAL8, b, blqtam*blqtam, 1, info)

        print*, 'Computo AABB'
        call multiplicarbloque(c,a,b,blqtam)
      enddo

      call pvmflvgroup('mmult',info)
      call pvmfexit(info)
    
      print*, 'Fin programa'
      stop
      end
        
      subroutine iniciarbloques (a, b, c, blk, row, col)
      implicit none
	    real*8 a(blk*blk), b(blk*blk), c(blk*blk)
      integer blk, row, col
      integer len
      integer i,j
        
      len = blk*blk
      do i=1, len 
        a(i)=i*(row*col+1)**2/blk**2
        c(i) = 0.0
      enddo
      do i=1,blk
        do j=1, blk
          if (row.eq.col) then
            if (i.eq.j) then
              b(j*blk+i)=1.0
            else
              b(j*blk+i)=0.0
            endif
          else
            b(j*blk+i) = 0.0
          endif
        enddo
      enddo
      end

      subroutine multiplicarbloque(c, a, b, blk) 
      implicit none
      integer blk
      integer i,j,k
      real*8 a(blk*blk), b(blk*blk), c(blk*blk)

      do i=1, blk
        do j=1,blk
          do k=1,blk
            c(i*blk+j) =c(i*blk+j) + (a(i*blk+k) * b(k*blk+j))
          enddo
        enddo
      enddo
      end
