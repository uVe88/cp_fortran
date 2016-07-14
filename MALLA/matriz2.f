      program matriz
      implicit none
      include 'fpvm3.h'
        
      integer ntask, info, mytid, mygid, i, blqtam, row, col, up
      integer down, MAXNTIDS, MAXM, ATAG, BTAG, DIMTAG
      integer MAXR, MAXBLOQTAM, mygtid, filas, tgid, ttag

      parameter(MAXNTIDS=16, MAXM=10, ATAG=2, BTAG=3, DIMTAG=5)
      parameter(MAXR = 4, MAXBLOQTAM = 100)
	  integer tids(MAXNTIDS-1)
	  integer mifila(MAXM)
      real*8 a(MAXBLOQTAM*MAXBLOQTAM), b(MAXBLOQTAM*MAXBLOQTAM)
      real*8 c(MAXBLOQTAM*MAXBLOQTAM), atmp(MAXBLOQTAM*MAXBLOQTAM)
	 
c     Obtengo el id de mi tarea
      call pvmfmytid(mytid)
      call pvmfcatchout(1,info)
      print*,'mi tid es:', mytid

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
        print *,'Introduce el valor de filas r (entero):'
        read *, filas

        print *,'Introduce el tamaño de bloque (entero):'
        read *, blqtam

c       Asigno el numero de procesos de la malla y compruebo que las tareas están en los límites de máximos tids
        ntask = filas*filas
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
        if ((filas.lt.1).or.(filas.ge.MAXR)) then
          print*,'ERROR r = ',filas,' not valid; max:', MAXR 
          call pvmflvgroup('mmult',info)
          call pvmfexit(info)
          stop
        endif

c       Saltamos la propagación en caso de ser una única tarea.
        if (ntask.eq.1) go to 100

c       Propagamos todas las tareas menos unas que es la mia (root)
        print*,'ntask:', ntask
        call pvmfspawn('matriz',PVMTASKDEFAULT,'*',ntask-1,tids,info)
        print*, 'ntids:', info

c       Comprobamos que se han lanzados todos los procesos de la malla
        if (info.ne.ntask-1) then
          print*,'ERROR!!!', info
          call pvmflvgroup('mmult', info)
          call pvmfexit(info)
          stop
        endif
           
c       Mandadamos la dimension de la matriz a todas las tareas
        print*, 'Envia root la dimensión de las matrices'
        call pvmfinitsend(PVMDATADEFAULT,info)
        call pvmfpack(INTEGER4, filas, 1, 1, info)
        call pvmfpack(INTEGER4, blqtam, 1, 1, info)
        call pvmfmcast(ntask-1, tids, DIMTAG, info)
        print*, 'Enviada la info de root'
      else
c       Soy miembro de la malla para pero no root
c       Obtengo mi id de grupo

        call pvmfgettid('mmult', 0, mygtid)
        print*, 'group id:', mygtid
        
c       Recibo los parámetros de la malla enviados por el root
        print*, 'Espero datos'
        call pvmfrecv(mygtid, DIMTAG, info)
        call pvmfunpack(INTEGER4, filas, 1, 1, info)
        call pvmfunpack(INTEGER4, blqtam, 1, 1, info)
        print*, 'Recibo los datos', filas, blqtam
c       Reinicio el número de tareas.       
        ntask = filas*filas
      endif

c     Sincronizamos con barrier
100   call pvmfbarrier('mmult',ntask, info)
      print*, 'paso el barrier'
      print*, 'filas', filas
      print*, 'mygid', mygid

c     Ordenos y busco los tids
      do i=1, filas
        tgid = (mygid/filas)*filas + i
        print*,'tgid:',tgid 
        call pvmfgettid('mmult', tgid, mygtid)
        mifila(i) = mygtid 
      enddo

      print*, 'Ordenos y busco los tids'
c     Busco mis posición
      row = mygid/filas
      col = mod(mygid,filas)
      print*,'IdGrupo, Fila, Columna',mygid,row,col

c     Calculo mis vecinos de arriba y abajo 
	  if (row.gt.0) then
        call pvmfgettid('mmult', row-1, up)
      else
       	call pvmfgettid('mmult', (filas-1)*filas+col, up)	   
      endif

	  if (row.eq.filas-1) then
	   	call pvmfgettid('mmult', col, down)
	  else
	   	call pvmfgettid('mmult', (row+1)*filas+col, down)
      endif
      print*, 'Calculo mis vecinos de arriba y abajo'

c     Inicializo los bloques
      call iniciarbloques(a, b, c, blqtam, row, col)
      print*, 'Inicializo los bloques'
c     Calculos la multiplicación de la matriz
      do i=1, filas
c     Envio los bloques de la matriz A
        if (col.eq.mod((row + i), filas)) then
          print*,'Entro en if'
          call pvmfinitsend(PVMDATADEFAULT,info)
          call pvmfpack(REAL8, a, blqtam*blqtam, 1, info)
          ttag = ATAG
          call pvmfmcast(filas, mifila, ttag, info)
          print*,'Entro en if - datos enviados'
          call multiplicarbloque(c,a,b,blqtam)
          print*,'Entro en if - multiplico datos'
        else
          print*,'Entro en else'
          tgid = row*filas+mod((row +i),filas)
          call pvmfgettid('mmult',tgid, mygtid)
          ttag = ATAG
          call pvmfrecv(mygtid, ttag, info)
          call pvmfunpack(REAL8, a, blqtam*blqtam, 1, info)
          print*,'Entro en else - datos recibidos'
          call multiplicarbloque(c,atmp,b,blqtam)
          print*,'Entro en else - multiplico datos'
        endif
        print*,'Salgo de if-else'

c     Roto las columnas de la matriz B
        call pvmfinitsend(PVMDATADEFAULT,info)
        call pvmfpack(REAL8, b, blqtam*blqtam, 1, info)
        ttag = BTAG
        call pvmfsend(up, (i+1)*BTAG, info)
        ttag = BTAG
        call pvmfrecv(down, (i+1)*BTAG, info)
        call pvmfunpack(REAL8, b, blqtam*blqtam, 1, info)
      enddo

c     Compruebo los resultados
      do i=0, blqtam*blqtam 
        if (a(i).ne.c(i)) then
          print*,'Error a(', i,') ', a(i),' != c(', i,') ', c(i)
        endif
      enddo

      call pvmflvgroup('mmult', info)
      call pvmfexit(info)

c     Imprimo el resultado
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
            c((i-1)*blk+j) =c((i-1)*blk+j) + (a(i-1)*blk+k) * b((k-1)*blk+j))
          enddo
        enddo
      enddo
      end
