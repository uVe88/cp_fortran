
       integer function norma (sum)
       norma = 0
       if (sum .lt. 1e-8 ) then
         norma = 1
       endif

       return
       end 