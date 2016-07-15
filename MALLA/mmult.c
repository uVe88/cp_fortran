/* 
    Matrix Multiply
*/

/* defines and prototypes for the PVM library */
#include <pvm3.h>
#include <stdio.h>

/* Maximum number of children this program will spawn */
#define MAXNTIDS    100
#define MAXROW      10

/* Message tags */
#define ATAG        2
#define BTAG        3
#define DIMTAG      5

void
InitBlock(float *a, float *b, float *c, int blk, int row, int col)
{
    int len, ind;
    int i,j;

    srand(pvm_mytid());
    len = blk*blk;
    for (ind = 0; ind < len; ind++) 
        { a[ind] = (float)(rand()%1000)/100.0; c[ind] = 0.0; }
    for (i = 0; i < blk; i++) {
        for (j = 0; j < blk; j++) {
            if (row == col)
                b[j*blk+i] = (i==j)? 1.0 : 0.0;
            else
                b[j*blk+i] = 0.0;
            }
        }
}

void
BlockMult(float* c, float* a, float* b, int blk) 
{
    int i,j,k;

    for (i = 0; i < blk; i++)
        for (j = 0; j < blk; j ++)
            for (k = 0; k < blk; k++)
                c[i*blk+j] += (a[i*blk+k] * b[k*blk+j]);
}

int
main(int argc, char* argv[])
{

    /* number of tasks to spawn, use 3 as the default */
    int ntask = 2;
    /* return code from pvm calls */
    int info;
    /* my task and group id */
    int mytid, mygid;
    /* children task id array */
    int child[MAXNTIDS-1];
    int i, m, blksize;
    /* array of the tids in my row */
    int myrow[MAXROW];
    float *a, *b, *c, *atmp;
    int row, col, up, down;

    /* find out my task id number */
    mytid = pvm_mytid();
    pvm_advise(PvmRouteDirect);

    /* check for error */
    if (mytid < 0) { 
        /* print out the error */
        pvm_perror(argv[0]); 
        /* exit the program */ 
        return -1;
        }

    /* join the mmult group */
    mygid = pvm_joingroup("mmult");
    if (mygid < 0) { 
        pvm_perror(argv[0]); pvm_exit(); return -1; 
        }

    /* if my group id is 0 then I must spawn the other tasks */
    if (mygid == 0) {
        /* find out how many tasks to spawn */
        if (argc == 3) {
            m = atoi(argv[1]);
            blksize = atoi(argv[2]);
            }
        if (argc < 3) {
            fprintf(stderr, "usage: mmult m blk\n");
            pvm_lvgroup("mmult"); pvm_exit(); return -1; 
            }

        /* make sure ntask is legal */
        ntask = m*m;
        if ((ntask < 1) || (ntask >= MAXNTIDS)) { 
            fprintf(stderr, "ntask = %d not valid.\n", ntask);
            pvm_lvgroup("mmult"); pvm_exit(); return -1; 
            }
        /* no need to spawn if there is only one task */
        if (ntask == 1) goto barrier;

        /* spawn the child tasks */
        info = pvm_spawn("mmult", (char**)0, PvmTaskDefault, (char*)0,
            ntask-1, child);

        /* make sure spawn succeeded */
        if (info != ntask-1) { 
            pvm_lvgroup("mmult"); pvm_exit(); return -1; 
            }

        /* send the matrix dimension */
        pvm_initsend(PvmDataDefault);
        pvm_pkint(&m, 1, 1);
        pvm_pkint(&blksize, 1, 1);
        pvm_mcast(child, ntask-1, DIMTAG);
        }
    else {
        /* recv the matrix dimension */
        pvm_recv(pvm_gettid("mmult", 0), DIMTAG);
        pvm_upkint(&m, 1, 1);
        pvm_upkint(&blksize, 1, 1);
        ntask = m*m;
        }

    /* make sure all tasks have joined the group */
barrier:
    info = pvm_barrier("mmult",ntask);
    if (info < 0) pvm_perror(argv[0]);

    /* find the tids in my row */
    for (i = 0; i < m; i++) 
        myrow[i] = pvm_gettid("mmult", (mygid/m)*m + i);

    /* allocate the memory for the local blocks */
    a = (float*)malloc(sizeof(float)*blksize*blksize);
    b = (float*)malloc(sizeof(float)*blksize*blksize);
    c = (float*)malloc(sizeof(float)*blksize*blksize);
    atmp = (float*)malloc(sizeof(float)*blksize*blksize);
    /* check for valid pointers */
    if (!(a && b && c && atmp)) { 
        fprintf(stderr, "%s: out of memory!\n", argv[0]);
        free(a); free(b); free(c); free(atmp);
        pvm_lvgroup("mmult"); pvm_exit(); return -1; 
        }

    /* find my block's row and column */
    row = mygid/m; col = mygid % m;
    /* calculate the neighbor's above and below */
    up = pvm_gettid("mmult", ((row)?(row-1):(m-1))*m+col);
    down = pvm_gettid("mmult", ((row == (m-1))?col:(row+1)*m+col));

    /* initialize the blocks */
    InitBlock(a, b, c, blksize, row, col);

    /* do the matrix multiply */
    for (i = 0; i < m; i++) {
        /* mcast the block of matrix A */
        if (col == (row + i)%m) {
            pvm_initsend(PvmDataDefault);
            pvm_pkfloat(a, blksize*blksize, 1);
            pvm_mcast(myrow, m, (i+1)*ATAG);
            BlockMult(c,a,b,blksize);
            }
        else {
            pvm_recv(pvm_gettid("mmult", row*m + (row +i)%m), (i+1)*ATAG);
            pvm_upkfloat(atmp, blksize*blksize, 1);
            BlockMult(c,atmp,b,blksize);
            }
        /* rotate the columns of B */
        pvm_initsend(PvmDataDefault);
        pvm_pkfloat(b, blksize*blksize, 1);
        pvm_send(up, (i+1)*BTAG);
        pvm_recv(down, (i+1)*BTAG);
        pvm_upkfloat(b, blksize*blksize, 1);
        }

    /* check it */
    for (i = 0 ; i < blksize*blksize; i++) 
        if (a[i] != c[i]) 
            printf("Error a[%d] (%g) != c[%d] (%g) \n", i, a[i],i,c[i]);

    printf("Done.\n");
    free(a); free(b); free(c); free(atmp);
    pvm_lvgroup("mmult");
    pvm_exit();
    return 0;
}