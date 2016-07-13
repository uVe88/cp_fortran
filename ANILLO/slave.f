c Programa slave - practica 1 programación paralela
c Alumno Ivan Gomez
      program slave
	implicit none
	intenger mytid, rtid, tids, ntids, var, mpos, ntid
        
c Inicio el pvm y obtengo mi tid
      call pvmfmytid( mytid )
      print*, 'Esclavo: ', mytid

c Obtengo la información
c     - el tip del nodo master - rtid
c     - número total de tids - ntids
c     - el array de tids - tids
c     - la variable a incrementar - var
      call pvmfrecv(mytid, 1, info)
      call pvmfunpack(INTEGER4, rtid, 1, 1, info)
      call pvmfunpack(INTEGER4, ntids, 1, 1, info)
      call pvmfunpack(INTEGER4, tids, ntids, 1, info)
      call pvmfunpack(INTEGER4, var, 1, 1, info)

c Incrementamos el valor de la variable de paso
      var = var + 1

      do i=1, ntids
        if (mytid.eq.tids(i)) then
          mpos=i
        endif
      enddo

c Compruebo si soy el último y envío al master el resultado o envío al siguiente tid
      if (mpos.eq.ntids) then
c Soy el último
        call pvmfinitsend(PVMDATADRAW, 1)
        call pvmfpack(INTEGER4, var, 1, 1, info)
        call pvmfsend(rtid, 1, info)
      else
c Envio al siguiente nodo del anillo
        call pvmfinitsend(PVMDATADRAW, 1)
        call pvmfpack(INTEGER4, mytid, 1, 1, info)
        call pvmfpack(INTEGER4, tids, nhost, 1, info)
        call pvmfpack(INTEGER4, nhost, 1, 1, info)
        call pvmfpack(INTEGER4, var, 1, 1, info)
        call pvmfsend(tids(mpos+1), 1, info)
      endif
      end
