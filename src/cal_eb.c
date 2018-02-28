#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double identify_tst(double*, int, double);

void cal_eb(double *eb, double *ef_tr, int n1, int n2, int crystal1, int crystal2, double *ef_atoms, int n_crystal_type, double eb0) {
    
    int i, ii;
    double tmp;
    // cluster2 -> cluster1
    for (i=0;i<n1+n2+1;i++)
        *(ef_tr+i) = 0.0;
    *(ef_tr + 0) = *(ef_atoms + n_crystal_type*(n1-1) + crystal1) + *(ef_atoms + n_crystal_type*(n2-1) + crystal2);
//    printf("ef_tr[%d] = %f\n",0,*(ef_tr+0));
    for (i=1;i<n2;i++)
        *(ef_tr + i) = *(ef_atoms + n_crystal_type*(n1+i-1) + crystal1) + *(ef_atoms + n_crystal_type*(n2-i-1) + crystal2);
    *(ef_tr + n2) = *(ef_atoms + n_crystal_type*(n1+n2-1) + crystal1);
//    for (i=0;i<n2+1;i++)
//        printf("ef_tr[%d] = %f\n",i,*(ef_tr+i));
    tmp = identify_tst(ef_tr,n2+1,eb0);
//    printf("tmp = %f\n",tmp);
    *(eb+crystal1) = tmp;
//    printf("crystal1->crystal2: eb[%d] = %f\n",crystal1,*(eb+crystal1));

    // crystal1 -> crystal2
//    ef_tr = (double*) realloc(ef_tr, (n1+1));
    *(ef_tr + 0) = *(ef_atoms + n_crystal_type*(n1-1) + crystal1) + *(ef_atoms + n_crystal_type*(n2-1) + crystal2);
    for (i=1;i<n1;i++)
        *(ef_tr + i) = *(ef_atoms + n_crystal_type*(n1-i-1) + crystal1) + *(ef_atoms + n_crystal_type*(n2+i-1) + crystal2);
    *(ef_tr + n1) = *(ef_atoms + n_crystal_type*(n1+n2-1) + crystal2);
//    for (i=0;i<n1+1;i++)
//        printf("ef_tr[%d] = %f\n",i,*(ef_tr+i));
    *(eb+crystal2) = identify_tst(ef_tr,n1+1,eb0);
//    printf("crystal2->crystal1: eb[%d] = %f\n",crystal2,*(eb+crystal2));
    
    // crystal 1&2 -> crystal3
    for (ii=0;ii<n_crystal_type;ii++) {
        if (ii==crystal1 || ii==crystal2) {
//            printf("%d crystal is existing: cry1=%d, cry2=%d\n",ii,crystal1,crystal2);
            continue;
        }
        else {
//            printf("%d crystal is new\n",ii);
//            ef_tr = (double*) realloc(ef_tr, (n1+n2+1));
            *(ef_tr + 0) = *(ef_atoms + n_crystal_type*(n1-1) + crystal1) + *(ef_atoms + n_crystal_type*(n2-1) + crystal2);
            for (i=1;i<n1;i++)
                *(ef_tr + i) = *(ef_atoms + n_crystal_type*(i-1) + ii) + *(ef_atoms + n_crystal_type*(n1-i-1) + crystal1) + *(ef_atoms + n_crystal_type*(n2-1) + crystal2);
            *(ef_tr + n1) = *(ef_atoms + n_crystal_type*(n1-1) + ii) + *(ef_atoms + n_crystal_type*(n2-1) + crystal2);
            for (i=1;i<n2;i++)
                *(ef_tr + n1+i) = *(ef_atoms + n_crystal_type*(n1+i-1) + ii) + *(ef_atoms + n_crystal_type*(n2-i-1) + crystal2);
            *(ef_tr + n1+n2) = *(ef_atoms + n_crystal_type*(n1+n2-1) + ii);
//            for (i=0;i<n1+n2+1;i++)
//                printf("ef_tr[%d] = %f\n",i,*(ef_tr+i));
            *(eb+ii) = identify_tst(ef_tr,n1+n2+1,eb0);
        }
    }
}



















