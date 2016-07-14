c Programa slave - practica 1 programación paralela
c Alumno Ivan Gomez
      program slave
	implicit none
      include 'fpvm3.h'
      integer maxtids, mpos, ntid, info, i, numt
      parameter(maxtids=8)
	integer mytid, rtid, tids(maxtids), ntids, var
        
c Inicio el pvm y obtengo mi tid
      call pvmfmytid(mytid)
      print*, 'Esclavo: ', mytid

c Obtengo la información
c     - el tip del nodo master - rtid
c     - número total de tids - ntids
c     - el array de tids - tids
c     - la variable a incrementar - var
      print*, 'Espero datos'
      call pvmfrecv(-1, 1, info)
      call pvmfunpack(INTEGER4, rtid, 1, 1, info)
      call pvmfunpack(INTEGER4, ntids, 1, 1, info)
      call pvmfunpack(INTEGER4, tids, ntids, 1, info)
      call pvmfunpack(INTEGER4, var, 1, 1, info)
      print*, "recibo parametros", rtid, ntids, tids, var

c Incrementamos el valor de la variable de paso
      var = var + 1

      do i=1, ntids
        if (mytid.eq.tids(i)) then
          mpos=i
        endif
      enddo

      print*, 'Mi posicion es:',mpos
c Compruebo si soy el último y envío al master el resultado o envío al siguiente tid
      if (mpos.eq.ntids) then
c Soy el último
        print*, 'Soy el último'
        call pvmfinitsend(PVMDATARAW, info)
        call pvmfpack(INTEGER4, var, 1, 1, info)
        call pvmfsend(rtid, 3, info)
        print*, 'Enviado a: ',rtid
      else
c Envio al siguiente nodo del anillo
        print*, 'No soy el último'
        print*, 'Enviando a: ',tids(mpos+1)
        call pvmfinitsend(PVMDATARAW, info)
        call pvmfpack(INTEGER4, rtid, 1, 1, info)
        call pvmfpack(INTEGER4, ntids, 1, 1, info)
        call pvmfpack(INTEGER4, tids, ntids, 1, info)
        call pvmfpack(INTEGER4, var, 1, 1, info)
        call pvmfsend(tids(mpos+1), 1, info)
        print*, 'Enviado a: ',tids(mpos+1)
      endif
      end
