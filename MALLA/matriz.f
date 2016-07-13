        program matriz
        implicit none
        include 'fpvm3.h'
        
        integer ntask, info, mytid, mygid, i, m, blksize, row, col, up
        integer down, MAXNTIDS, MAXROW, ATAG, BTAG, DIMTAG
        integer r, blktam, mygtid

c        /* Maximum number of children this program will spawn */
        parameter(MAXNTIDS=100, MAXROW=10, ATAG=2, BTAG=3, DIMTAG=5)
        parameter(r = 4, blktam = 4)
	 integer child(MAXNTIDS-1)
	 integer myrow(MAXROW)
        real*8 a(blktam*blktam), b(blktam*blktam)
        real*8 c(blktam*blktam), atmp(blktam*blktam)
	 
c       number of tasks to spawn, use 3 as the default
        ntask = 2

c       children task id array

        
c       array of the tids in my row
        

c       Obtengo el id de mi tarea
        call pvmfmytid(mytid)
c       pvm_advise(PvmRouteDirect)

c       Verificamos que pvm ha iniciado bien.
        if (mytid.lt.0) then
c           print out the error
            print* , 'Error, tid erroneo' 
c           exit the program
            stop
        endif

c       Nos unimos al grupo mmult
        call pvmfjoingroup('mmult', mygid)

        if (mygid.lt.0) then
            print*,'Error al obtener el id de grupo'
            call pvmfexit()
            stop
        endif
        
c       /* if my group id is 0 then I must spawn the other tasks */
        if (mygid.eq.0) then
c           find out how many tasks to spawn
            
            m = r
            blksize = blktam
            
c           make sure ntask is legal
            ntask = m*m
            if ((ntask.lt.1).or.(ntask.ge.MAXNTIDS)) then
                print*,'ntask = ',ntask,' not valid'
                call pvmflvgroup('mmult',info)
                call pvmfexit(info)
                stop
            endif

c           no need to spawn if there is only one task
            if (ntask.eq.1) go to 100

c           /* spawn the child tasks */
            call pvmfspawn('mmult',PVMTASKDEFAULT, '*', ntask-1, child
     *,info)
            
c           /* make sure spawn succeeded */
            if (info.ne.ntask-1) then
                call pvmflvgroup('mmult', info)
                call pvmfexit(info)
                stop 
            endif
            
c           send the matrix dimension
            call pvmfinitsend(PVMDATADEFAULT,info)
            call pvmfpack(INTEGER4, m, 1, 1, info)
            call pvmfpack(INTEGER4, blksize, 1, 1, info)
            call pvmfmcast(child, ntask-1, DIMTAG, info)
        else
c           recv the matrix dimension
            call pvmfgettid('mmult', 0, mygtid)
            call pvmfrecv(mygtid, DIMTAG, info)
            call pvmfunpack(INTEGER4, m, 1, 1, info)
            call pvmfunpack(INTEGER4, blksize, 1, 1, info)
            ntask = m*m
        endif
c       make sure all tasks have joined the group */

c       barrier

100     call pvmfbarrier('mmult',ntask, info)

c        if (info.lt.0) then
c            pvm_perror(argv[0])

c        /* find the tids in my row */
        do i=0, m
            call pvmfgettid('mmult', (mygid/m)*m + i, mygtid)
            myrow(i) = mygtid 
        enddo
c       find my block's row and column
        row = mygid/m 
        col = mod(mygid,m)
c       calculate the neighbor's above and below
	 
	    if (row.gt.0) then
        	call pvmfgettid('mmult', row-1, up)
    	else
        	call pvmfgettid('mmult', (m-1)*m+col, up)	   
        endif

	    if (row.eq.m-1) then
	    	call pvmfgettid('mmult', col, down)
	    else
	    	call pvmfgettid('mmult', (row+1)*m+col, down)
	    endif

c       initialize the blocks
        call InitBlock(a, b, c, blksize, row, col)
c        do the matrix multiply

        do i=0, m
c           mcast the block of matrix A
            if (col.eq.mod((row + i), m)) then
                call pvmfinitsend(PVMDATADEFAULT,info)
                call pvmfpack(REAL8, a, blksize*blksize, 1, info)
                call pvmfmcast(myrow, m, (i+1)*ATAG, info)
                call BlockMult(c,a,b,blksize)
            else
                call pvmfgettid('mmult',row*m+mod((row +i),m), mygtid)
                call pvmfrecv(mygtid, (i+1)*ATAG, info)
                call pvmfunpack(REAL8, a, blksize*blksize, 1, info)
                call BlockMult(c,atmp,b,blksize)
            endif
c           rotate the columns of B
            call pvmfinitsend(PVMDATADEFAULT,info)
            call pvmfpack(REAL8, b, blksize*blksize, 1, info)
            call pvmfsend(up, (i+1)*BTAG, info)
            call pvmfrecv(down, (i+1)*BTAG, info)
            call pvmfunpack(REAL8, b, blksize*blksize, 1, info)
        enddo

c       check it
        do i=0, blksize*blksize 
            if (a(i).ne.c(i)) then
                print*,'Error a(', i,') ', a(i),' != c(', i,') ', c(i)
            endif
        enddo

        call pvmflvgroup('mmult', info)
        call pvmfexit(info)
        end
        
        subroutine InitBlock (a, b, c, blk, row, col)
        implicit none
	    real*8 a(blk*blk), b(blk*blk), c(blk*blk)
        integer blk, row, col
        integer len, ind
        integer i,j
        
        len = blk*blk
        do ind=1, len 
            a(ind)=ind*(row*col+1)**2/blk**2
            c(ind) = 0.0
        enddo
        do i=1,blk
            do j=1, blk
                if (row.eq.col) then
                    if (i.eq.j) then
                        b(j*blk+i)=1.0
                    else
                        b(j*blk+i)=0.0
                    endif
                else
                    b(j*blk+i) = 0.0
                endif
            enddo
        enddo
        end

        subroutine BlockMult(c, a, b, blk) 
        implicit none
        integer blk
        integer i,j,k
        real*8 a(blk*blk), b(blk*blk), c(blk*blk)

        do i=1, blk
            do j=1,blk
                do k=1,blk
                    c(i*blk+j) =c(i*blk+j) + (a(i*blk+k) * b(k*blk+j))
                enddo
            enddo
        enddo
        end
        
