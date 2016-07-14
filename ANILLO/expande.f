c Subrutine expande - practica 1 programación paralela
c Alumno Ivan Gomez
      subroutine expande(ejec, tids, nhost, numt)
      include 'fpvm3.h'
      integer nhost, numt, info, narch, dtid, speed, mytid
      character*30 ejec, arch, name
      integer tids(*)

      call pvmfmytid(mytid)
c      call pvmfconfig(nhost, narch, dtid, name, arch, speed, info)
	print*, 'Numero de host: ', nhost
      call pvmfspawn(ejec,PVMTASKDEFAULT,'*',nhost,tids,numt)
      return
	end
