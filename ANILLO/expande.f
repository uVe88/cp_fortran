c Subrutine expande - practica 1 programación paralela
c Alumno Ivan Gomez
      subroutine expande(ejec, tids, nhost, numt)
      include 'fpvm3.h'
      integer nhost, numt, info, narch, dtid, speed
      character*30 ejec, arch, name
      call pvmfconfig(nhost, narch, dtid, name, arch, speed, info)
	  print*, 'Numero de host: ', nhost
      call pvmfspawn(ejec, PVMTASKDEAFULT, '*', nhost, tids, numt)
	  end
