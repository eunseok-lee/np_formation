#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double identify_tst(double *ef_tr, int ef_tr_size, double eb0) {
    
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
        eb = fmax((*(ef_tr+0) + *(ef_tr+1))/2+eb0,*(ef_tr+0)) - *(ef_tr+0);
    else if (max_id == ef_tr_size-1)
        eb = fmax((*(ef_tr+ef_tr_size-2) + *(ef_tr+ef_tr_size-1))/2+eb0,*(ef_tr+0)) - *(ef_tr+0);
    else
        eb = fmax(fmax((*(ef_tr+max_id-1) + *(ef_tr+max_id))/2,(*(ef_tr+max_id) + *(ef_tr+max_id+1))/2)+eb0,*(ef_tr+0)) - *(ef_tr+0);
//    printf("calcualted eb=%f\n",eb);
    return(eb);
}

