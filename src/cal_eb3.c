#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double identify_tst(double*, int, double);
double identify_tst2(double*, int, double*, int, int, int, int);

void cal_eb3(double *eb, double *ef_tr, int n1, int n2, int crystal1, int crystal2, double *ef_atoms, int n_crystal_type, double *eb0, int natoms_c) {
    
    int i, ii;
    double tmp;
    
    int ef_n1[n1+2], ef_n2[n2+2];
    int len_ef_n1, len_ef_n2;
    double eb12, eb21;
    
    if (n1<=natoms_c) {
        len_ef_n1 = 2;
        ef_n1[0] = n1;
        ef_n1[1] = 0;
    }
    else {
        len_ef_n1 = n1-natoms_c+2;
        for (i=0;i<=n1-natoms_c;i++)
            ef_n1[i] = n1-i;
        ef_n1[n1-natoms_c+1] = 0;
    }
    if (n2<=natoms_c) {
        len_ef_n2 = 2;
        ef_n2[0] = n2;
        ef_n2[1] = 0;
    }
    else {
        len_ef_n2 = n2-natoms_c+2;
        for (i=0;i<=n2-natoms_c;i++)
            ef_n2[i] = n2-i;
        ef_n2[n2-natoms_c+1] = 0;
    }
    
    // cluster2 -> cluster1
    for (i=0;i<n1+n2+1;i++)
        *(ef_tr+i) = 0.0;
    *(ef_tr + 0) = *(ef_atoms + n_crystal_type*(n1-1) + crystal1) + *(ef_atoms + n_crystal_type*(n2-1) + crystal2);
//    printf("ef_tr[%d] = %f\n",0,*(ef_tr+0));
    for (i=1;i<len_ef_n2-1;i++)
        *(ef_tr + i) = *(ef_atoms + n_crystal_type*(n1+n2-ef_n2[i]-1) + crystal1) + *(ef_atoms + n_crystal_type*(ef_n2[i]-1) + crystal2);
    *(ef_tr + len_ef_n2-1) = *(ef_atoms + n_crystal_type*(n1+n2-1) + crystal1);
//    for (i=0;i<len_ef_n2;i++)
//        printf("ef_tr[%d] = %f\n",i,*(ef_tr+i));
    eb21 = identify_tst(ef_tr,len_ef_n2,eb0[crystal2]);
//    printf("cluster2(%d)->cluster1(%d): eb = %f\n",crystal2,crystal1,eb21);

    // crystal1 -> crystal2
//    ef_tr = (double*) realloc(ef_tr, (n1+1));
    for (i=0;i<n1+n2+1;i++)
        *(ef_tr+i) = 0.0;
    *(ef_tr + 0) = *(ef_atoms + n_crystal_type*(n1-1) + crystal1) + *(ef_atoms + n_crystal_type*(n2-1) + crystal2);
    for (i=1;i<len_ef_n1-1;i++)
        *(ef_tr + i) = *(ef_atoms + n_crystal_type*(n2+n1-ef_n1[i]-1) + crystal2) + *(ef_atoms + n_crystal_type*(ef_n1[i]-1) + crystal1);
    *(ef_tr + len_ef_n1-1) = *(ef_atoms + n_crystal_type*(n1+n2-1) + crystal2);
//    for (i=0;i<len_ef_n1;i++)
//        printf("ef_tr[%d] = %f\n",i,*(ef_tr+i));
    eb12 = identify_tst(ef_tr,len_ef_n1,eb0[crystal1]);
//    printf("cluster1(%d)->cluster2(%d): eb = %f\n",crystal1,crystal2,eb12);
    
    if (crystal1==crystal2)
        *(eb + crystal1) = fmin(eb21,eb12);
    else {
        *(eb + crystal1) = eb21;
        *(eb + crystal2) = eb12;
    }
    
    // crystal 1&2 -> crystal3
    for (ii=0;ii<n_crystal_type;ii++) {
        if (ii==crystal1 || ii==crystal2) {
//            printf("%d crystal is existing: cry1=%d, cry2=%d\n",ii,crystal1,crystal2);
            continue;
        }
        else {
//            printf("%d crystal is new\n",ii);
//            ef_tr = (double*) realloc(ef_tr, (n1+n2+1));
            // cluster 1 and then 2 transforms
            *(ef_tr + 0) = *(ef_atoms + n_crystal_type*(n1-1) + crystal1) + *(ef_atoms + n_crystal_type*(n2-1) + crystal2);
            for (i=1;i<len_ef_n1-1;i++)
                *(ef_tr + i) = *(ef_atoms + n_crystal_type*(n1-ef_n1[i]-1) + ii) + *(ef_atoms + n_crystal_type*(ef_n1[i]-1) + crystal1) + *(ef_atoms + n_crystal_type*(n2-1) + crystal2);
            *(ef_tr + len_ef_n1-1) = *(ef_atoms + n_crystal_type*(n1-1) + ii) + *(ef_atoms + n_crystal_type*(n2-1) + crystal2);
            for (i=1;i<len_ef_n2-1;i++)
                *(ef_tr + len_ef_n1-1+i) = *(ef_atoms + n_crystal_type*(n1+n2-ef_n2[i]-1) + ii) + *(ef_atoms + n_crystal_type*(ef_n2[i]-1) + crystal2);
            *(ef_tr + len_ef_n1+len_ef_n2-2) = *(ef_atoms + n_crystal_type*(n1+n2-1) + ii);
//            for (i=0;i<len_ef_n1+len_ef_n2-1;i++)
//                printf("ef_tr[%d] = %f\n",i,*(ef_tr+i));
            eb12 = identify_tst2(ef_tr,len_ef_n1+len_ef_n2-1,eb0,crystal1,len_ef_n1,crystal2,len_ef_n2);
            
            // cluster 2 and then 1 transforms
            *(ef_tr + 0) = *(ef_atoms + n_crystal_type*(n1-1) + crystal1) + *(ef_atoms + n_crystal_type*(n2-1) + crystal2);
            for (i=1;i<len_ef_n2-1;i++)
                *(ef_tr + i) = *(ef_atoms + n_crystal_type*(n2-ef_n2[i]-1) + ii) + *(ef_atoms + n_crystal_type*(ef_n2[i]-1) + crystal2) + *(ef_atoms + n_crystal_type*(n1-1) + crystal1);
            *(ef_tr + len_ef_n2-1) = *(ef_atoms + n_crystal_type*(n2-1) + ii) + *(ef_atoms + n_crystal_type*(n1-1) + crystal1);
            for (i=1;i<len_ef_n1-1;i++)
                *(ef_tr + len_ef_n2-1+i) = *(ef_atoms + n_crystal_type*(n2+n1-ef_n1[i]-1) + ii) + *(ef_atoms + n_crystal_type*(ef_n1[i]-1) + crystal1);
            *(ef_tr + len_ef_n1+len_ef_n2-2) = *(ef_atoms + n_crystal_type*(n1+n2-1) + ii);
            //            for (i=0;i<len_ef_n1+len_ef_n2-1;i++)
            //                printf("ef_tr[%d] = %f\n",i,*(ef_tr+i));
            eb21 = identify_tst2(ef_tr,len_ef_n1+len_ef_n2-1,eb0,crystal2,len_ef_n2,crystal1,len_ef_n1);
            
            *(eb + ii) = fmin(eb12,eb21);
        }
    }
}






















