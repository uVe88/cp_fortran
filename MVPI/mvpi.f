      program mvpi
      implicit none
      integer maxm, maxn, maxntids, info, ibuf, nproc
      Parameter(maxm=100, maxn=100, maxntids=8)
      integer mytid, ptid
      real*8 a(maxm, maxn), x(maxn), res(maxm, maxn)
      integer tids(1:maxntids)
      integer resto, ucom
      integer m, n, pindex, iup, idown
c     Inicializo pvm y obtengo mi tid y el tid de mi padre
      
      call pvmfmytid(mytid)
      call pvmfparent(ptid)

      m = maxm
      n = maxn

c     Compruebo si soy el padre
      if(ptid.eq.pvmnoparent) then
c       Soy el padre
        call defmatrizvector(a, x, maxm, maxn)
        call pvmfspawn('mvpi',PVMTASKDEFAULT, '*', maxntids, tids, nproc)

        ucom = m / nproc
        resto = mod(m, nproc)

        if (ucom.gt.0) then
c       Envio los datos a computar.
            iup = ucom
            idown = 1
            do i=1, nproc
                
                call pvmfinitsend(PVMDATAINPLACE, ibuf)
                call pvmfpack(INTEGER4, idown, 1, 1, info)
                call pvmfpack(INTEGER4, iup, 1, 1, info)
                call pvmfpack(INTEGER4, ucom, 1, 1, info)
                call pvmfpack(REAL8, x(1), maxn, 1, info)

                do i=idown, ucom
                    call pvmfpack(REAL8, a(i,1), maxn, 1, info)
                enddo
                
                call pvmfsend(tids(i), 2, info)

                idown = iup +1
                iup = iup + ucom
            enddo

        endif

        if (resto.gt.0) then
          idown = ucom * mproc
          iup = m
          call matrizvector(a, x, res, iup, idown, m, n, info)
        endif

c     Recibo los datos de los procesos hijo y reconstruyo los resultados
      

      else
c     Soy hijo
c     Recibo los datos
        call pvmfrecv()
        call pvmfpack(INTEGER4, idown, 1, 1, info)
        call pvmfpack(INTEGER4, iup, 1, 1, info)
        call pvmfpack(INTEGER4, ucom, 1, 1, info)
        call pvmfpack(REAL8, x(1), n, 1, info)

        do i=idown, ucom
            call pvmfpack(REAL8, a(i,1), n, 1, info)
        enddo
        
        call pvmfsend(tids(i), tids(i), info)
      endif

      end
      
      subroutine defmatrizvector(a, x, m, n)
      implicit none
      real*8 a(m, n), x(n)
      integer m, n, info
      
      do j=1, m
        do i=1, n
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
      real*8 res(iup-idown, n)
      info = 0

      if (index.gt.n) then
        print*, 'Error: indice > n'
        info = -1
      else
        print*, 'Calculo de la matrix por vector del índice,', index
        do i=idown, iup-idown
            do j=1, n
                res(i, j) = a(i, j) * x(j)
            enddo
      endif
      end

      subroutine shutdown(nproc, tids)
      integer nproc, tids(*)

      do i=1, nproc
         call pvmfkill( tids(i), info )
      enddo
      call pvmfexit( info )
      return
      end