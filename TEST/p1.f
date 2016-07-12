c Programa fortran que utiliza sumador
	programa p1
	implicit none
	call suma(5, 7)
	end

	subroutine suma(a,b)
	integer a, b, res
	res = a + b
	print*, 'Sumador: ', res
	end
	
