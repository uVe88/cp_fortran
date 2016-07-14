c Subrutine expande - practica 1 programación paralela
c Alumno Ivan Gomez
      subroutine expande(ejec, tids, nhost, numt)
      include 'fpvm3.h'
      integer nhost, numt, info, narch, dtid, speed, mytid
      character*30 ejec, arch, name
      integer tids(*), ttids(100), i, j, ntask
c Una tarea por host
      ntask=1
 100  format(1x,A20,3x,i7,3x,A7,3x,i7)
      call pvmfmytid(mytid)
      print*, 'Expande tid:', mytid
      call pvmfconfig(nhost, narch, dtid, name, arch, speed, info)
      
      print*,'------ Master host -------'
      write(*,100) name, dtid, arch, speed
      print*,'--------------------------'

      print*, 'Numero de hosts en PVM: ', nhost - 1

      if (nhost-1.gt.0) then
        print*, 'Propaga a los host'
        do i=1, nhost
          call pvmfconfig(nhost,narch,dtid,name,arch,speed,info)

          write(*,100) name, dtid, arch, speed
          call pvmfspawn(ejec,PVMTASKHOST,dtid,ntask,ttids,numt)
          
          do j=1, numt
            tids(i)=ttids(j)
          enddo
        enddo
      endif
      nhost = nhost-1
	end
