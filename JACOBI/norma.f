
       integer function norma (xsol,x,mmax,n)
       real*8 xsol(mmax), x(mmax),suma
       integer mmax,n,i 
       suma = 0
       do i=1,n
         suma = suma +abs(xsol(i)-x(i)) 
       end do 
       norma = 0
       if (sum .lt. 1e-8 ) then
         norma = 1
       endif

       return
       end 