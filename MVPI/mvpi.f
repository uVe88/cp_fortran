      program mvpi
      implicit none
      include 'fpvm3.h'

      integer maxm, maxn, maxntids, info, ibuf, nproc, mtask
      Parameter(maxm=100, maxn=100, maxntids=8)
      integer mytid, ptid, ctid
      real*8 a(maxm, maxn), x(maxn), res(maxm*maxn)
      integer tids(1:maxntids)
      integer resto, ucom, i , j, c
      integer m, n, pindex, iup, idown
c     Inicializo pvm y obtengo mi tid y el tid de mi padre
      
      call pvmfmytid(mytid)
      call pvmfparent(ptid)
      call pvmfcatchout(1,info)

      m = maxm
      n = maxn
      mtask = maxntids

      m = 10
      n = 10
      mtask = 1

c     Compruebo si soy el padre
      if(ptid.eq.pvmnoparent) then
c       Soy el padre
        call defmatrizvector(a, x, m, n)
        call pvmfspawn('mvpi',PVMTASKDEFAULT, '*', mtask, tids, 
     *nproc)

c     Calculo las unidades de computacion por proceso y el restante de filas para el hilo principal
        ucom = m / nproc
        resto = mod(m, nproc)

        print*, 'Número de procesos: ', nproc
        print*, 'Unidades de computación: ', ucom
        print*, 'Resto para el padre: ', resto
        print*, x
        print*, 'Valores de la matriz'
        call escribe(a, m, n)

        if (ucom.gt.0) then
c       Envio los datos a computar.
            iup = ucom
            idown = 1
            do c=1, nproc
                call pvmfinitsend(PVMDATARAW, info)
                call pvmfpack(INTEGER4, idown, 1, 1, info)
                call pvmfpack(INTEGER4, iup, 1, 1, info)
                call pvmfpack(INTEGER4, ucom, 1, 1, info)
                call pvmfpack(REAL8,x,n,1,info)

                do i=1, m
                    call pvmfpack(REAL8, a(i,1), n, maxn, info)
                enddo
                
                call pvmfmcast(1, tids(c), 2, info)

                idown = iup +1
                iup = iup + ucom
            enddo
        endif
        print*, 'Padre-Valores Enviados'

c     Recibo los datos de los procesos hijo y reconstruyo los resultados
        do i=1, nproc
            print*, 'Padre-Espero valores'
            call pvmfrecv(-1, 3, info)
            call pvmfunpack(INTEGER4, idown, 1, 1, info)
            call pvmfunpack(INTEGER4, iup, 1, 1, info)
            call pvmfunpack(INTEGER4, ucom, 1, 1, info)
            call pvmfunpack(INTEGER4, ctid, 1, 1, info)
            do j=idown, ucom
                call pvmfunpack(REAL8,res(j),maxm*maxn, 1, info)
            enddo
            print*, 'Padre-Valores recibidos de:', ctid
        enddo

        if (resto.gt.0) then
          idown = ucom * nproc
          iup = m
          call matrizvector(a, x, res, iup, idown, m, n, info)
        endif

        print*, ('y(',i,')=',res(i),'; ',i=1,m*n)

      else
c     Soy hijo
c     Recibo los datos
        call pvmfrecv(ptid, 2, info)
        call pvmfunpack(INTEGER4, idown, 1, 1, info)
        call pvmfunpack(INTEGER4, iup, 1, 1, info)
        call pvmfunpack(INTEGER4, ucom, 1, 1, info)
        call pvmfunpack(REAL8,x,n,1, info)
        print*, 'Hijo - valor de m y n: ', m, n
        do i=1, m
            call pvmfunpack(REAL8, a(i,1), n, maxn, info)
            print*, 'Hijo - valor matrix i: ',i
            print*, ('a(',i,',',j,')=',a(i, j),'; ',j=1,n)
        enddo

        print*, 'Hijo - idown: ', idown
        print*, 'Hijo - iup: ', iup
        print*, 'Hijo - ucom: ', ucom
        print*, x
        print*, 'Hijo - Valores de la matriz'
        
        call escribe(a, m, n)
        
        call matrizvector(a, x, res, iup, idown, m, n, info)

c     Envio datos al padre ya computados
        call pvmfinitsend(PVMDATAINPLACE, ibuf)
        call pvmfpack(INTEGER4, idown, 1, 1, info)
        call pvmfpack(INTEGER4, iup, 1, 1, info)
        call pvmfpack(INTEGER4, ucom, 1, 1, info)
        call pvmfpack(INTEGER4, mytid, 1, 1, info)
        do i=idown, ucom
            call pvmfpack(REAL8,res(i),maxm*maxn, 1, info)
        enddo

        call pvmfsend(ptid, 3, info)
      endif
      
      call pvmfexit(info)
      end   
      
      subroutine defmatrizvector(a, x, maxm, maxn)
      implicit none
      integer info, i, j, maxm, maxn
      real*8 a(maxm, maxn), x(maxn)
      
      do j=1, maxn
        do i=1, maxm
              a(i,j) = (i+j)**2
        enddo
        x(j) = j**2
      enddo     
      end

      subroutine matrizvector(a, x, res, iup, idown, m, n, info)
      implicit none
      real*8 a(m, n), x(n)
      integer m, n, info
      integer iup, idown
      integer r, i, j
      real*8 res(*)
      info = 0

      if (idown.gt.n) then
        print*, 'Error: indice > n'
        info = -1
      else
        
        print*, 'proceso:', idown, 'indice:', idown
        do i=idown, iup-idown
            do j=1, n
                res(j + (i-1)*n) = a(i, j) * x(j)
            enddo
        enddo
      endif
      end

      subroutine shutdown(nproc, tids)
      integer nproc, tids(*), i
      do i=1, nproc
         call pvmfkill( tids(i), info )
      enddo
      call pvmfexit(info)
      end

       subroutine escribe(a,m,n)
c
c  ****  Escribe la matriz a con m filas y n columnas. 
c  ****  La primera dimension de a es maxm
c
       integer m,n
       real*8 a(m,n)
       do i=0,n
          write(*,100) i
       enddo
       write(*,'(:)')
       do i=1,m
          write(*,100) i
          do j=1,n
              write(*,200) a(i,j)
          enddo
          write(*,'(:)')
       enddo
 200   format($,F7.3)
 100   format($,I7)
       return
       end