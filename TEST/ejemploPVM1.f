	program ejemploPVM1
	implicit none
	include 'fpvm3.h'
	integer mytid, info
	call pvmfmytid(mytid)
	print*, 'Mi número de identificación de proceso es: ', mytid

	if (mytid .lt. 0) then
		write(*,*) 'Error al entrar en PVM'
		call pvmfexit(info)
		stop
	endif
	call pvmfexit(info)
	print*, 'Ha salido de PVM; codigo de estado', info
	end

