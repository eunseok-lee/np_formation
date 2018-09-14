#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double identify_tst2(double *ef_tr, int ef_tr_size, int c3, double *eb0, int c1, int n1, int c2, int n2) {
    
    int i, max_id;
    double max_val, eb;
    
    max_id = 0;
    max_val = *(ef_tr+0);
    for (i=1;i<ef_tr_size;i++)
        if (*(ef_tr+i) > max_val) {
            max_id = i;
            max_val = *(ef_tr+i);
        }
//    printf("max_id = %d, max_val = %f\n",max_id,max_val);
    if (max_id == 0)
        eb = fmax((*(ef_tr+0) + *(ef_tr+1))/2+fmax(eb0[c1],eb0[c3]),*(ef_tr+0)) - *(ef_tr+0);
    else if (max_id == ef_tr_size-1)
        eb = fmax((*(ef_tr+ef_tr_size-2) + *(ef_tr+ef_tr_size-1))/2+fmax(eb0[c2],eb0[c3]),*(ef_tr+0)) - *(ef_tr+0);
    else if (max_id < n1)
        eb = fmax(fmax((*(ef_tr+max_id-1) + *(ef_tr+max_id))/2,(*(ef_tr+max_id) + *(ef_tr+max_id+1))/2)+fmax(eb0[c1],eb0[c3]),*(ef_tr+0)) - *(ef_tr+0);
    else
        eb = fmax(fmax((*(ef_tr+max_id-1) + *(ef_tr+max_id))/2,(*(ef_tr+max_id) + *(ef_tr+max_id+1))/2)+fmax(eb0[c2],eb0[c3]),*(ef_tr+0)) - *(ef_tr+0);
//    printf("calcualted eb=%f\n",eb);
    return(eb);
}

