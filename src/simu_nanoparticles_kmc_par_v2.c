#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <time.h>
#include "mpi.h"

void load_double_mat(double*, int, int, char*);
double cal_eb(double*,double*,int,int,int,int,double*,int,double);
double cal_eb3(double*,double*,int,int,int,int,double*,int,double*,int);
double identify_tst(double*, int, double);
int nchoose2(int*,int,int*);

#define MASTER 0
#define FROM_MASTER 1
#define FROM_WORKER 2

int main(int argc, char **argv)
{
    int i, j, k, l, ni, nj, ctr, ctr_tmp;
    int p1, p2, n1, n2, crystal1, crystal2;
    int num_cluster, row_size;
    int cid1, cid2;
    int count;
    int KMCstep;
    double t;
    double ef;
    double tmp, tmp1, tmp2;
    FILE *fp, *fp2;
    char buff_line[200], dummy[200], param_name[100];
    int param_name_len;
    
    int Nparticles, Natoms_per_particle_ini, Natoms_instant_trans;
    int maxKMCstep;
    int dispfreq;
    double T_ini, T_end;
    double *kT;
    double *T_profile;
    double eb0[2];
    int store_range_to_finalstep;

    // mpi parameters and initialization
    int numprocs, rank, mtype;
    int row_dist_size;
    int row_ini, row_end, row_offset;
    int pr_tmp_size;
    double *pr_tmp, *pr_tmp_master;
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);  // Get # processors
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);      // Get my rank (id)

	if (numprocs > 1)
		if (rank == MASTER)
			printf("Parallel mode: %d processes are running\n",numprocs);

	if (rank == MASTER)
		srand(time(NULL));


    char paramfilename[100];
    char T_filename[100];
    T_filename[0] = '\0';
//    if (argc==1) {
//        printf("Issue: no parameter name, param.dat is used\n");
    strcpy(paramfilename,"param.dat");
//    }
//    else
//        strcpy(paramfilename,argv[1]);

    fp = fopen(paramfilename, "r");
    while(fgets(buff_line,sizeof(buff_line),fp) != NULL) {
        if (buff_line[0] == '#') {
//            printf("Comment line: %s",buff_line);
            continue;
        }
        strcpy(param_name,"Nparticles");
        param_name_len = strlen(param_name);
        strncpy(dummy,buff_line,param_name_len);
        dummy[param_name_len] = '\0';
        if (strcmp(dummy,param_name)==0) {
            sscanf(buff_line,"%s %d",dummy,&Nparticles);
            if (rank == MASTER)
                printf("Nparticles = %d\n", Nparticles);
        }
        strcpy(param_name,"Natoms_per_particle_ini");
        param_name_len = strlen(param_name);
        strncpy(dummy,buff_line,param_name_len);
        dummy[param_name_len] = '\0';
        if (strcmp(dummy,param_name)==0) {
            sscanf(buff_line,"%s %d",dummy,&Natoms_per_particle_ini);
            if (rank == MASTER)
                printf("Natoms_per_particle_ini = %d\n",Natoms_per_particle_ini);
        }
        strcpy(param_name,"Natoms_instant_trans");
        param_name_len = strlen(param_name);
        strncpy(dummy,buff_line,param_name_len);
        dummy[param_name_len] = '\0';
        if (strcmp(dummy,param_name)==0) {
            sscanf(buff_line,"%s %d",dummy,&Natoms_instant_trans);
            if (rank == MASTER)
                printf("Natoms_instant_trans = %d\n",Natoms_instant_trans);
        }
        strcpy(param_name,"maxKMCstep");
        param_name_len = strlen(param_name);
        strncpy(dummy,buff_line,param_name_len);
        dummy[param_name_len] = '\0';
        if (strcmp(dummy,param_name)==0) {
            sscanf(buff_line,"%s %d",dummy,&maxKMCstep);
            if (rank == MASTER)
                printf("maxKMCstep = %d\n",maxKMCstep);
        }
        strcpy(param_name,"T_ini");
        param_name_len = strlen(param_name);
        strncpy(dummy,buff_line,param_name_len);
        dummy[param_name_len] = '\0';
        if (strcmp(dummy,param_name)==0) {
            sscanf(buff_line,"%s %lf",dummy,&T_ini);
            if (rank == MASTER)
                printf("T_ini = %f\n",T_ini);
        }
        strcpy(param_name,"T_end");
        param_name_len = strlen(param_name);
        strncpy(dummy,buff_line,param_name_len);
        dummy[param_name_len] = '\0';
        if (strcmp(dummy,param_name)==0) {
            sscanf(buff_line,"%s %lf",dummy,&T_end);
            if (rank == MASTER)
                printf("T_end = %f\n",T_end);
        }
        strcpy(param_name,"T_filename");
        param_name_len = strlen(param_name);
        strncpy(dummy,buff_line,param_name_len);
        dummy[param_name_len] = '\0';
        if (strcmp(dummy,param_name)==0) {
            sscanf(buff_line,"%s %s",dummy,T_filename);
            if (rank == MASTER)
                printf("T_filename = %s\n",T_filename);
        }
        strcpy(param_name,"eb0");
        param_name_len = strlen(param_name);
        strncpy(dummy,buff_line,param_name_len);
        dummy[param_name_len] = '\0';
        if (strcmp(dummy,param_name)==0) {
            sscanf(buff_line,"%s %lf %lf",dummy,&tmp1,&tmp2);
            eb0[0] = tmp1;
            eb0[1] = tmp2;
            if (rank == MASTER)
                printf("eb0 = %f %f\n",eb0[0],eb0[1]);
        }
        strcpy(param_name,"dispfreq");
        param_name_len = strlen(param_name);
        strncpy(dummy,buff_line,param_name_len);
        dummy[param_name_len] = '\0';
        if (strcmp(dummy,param_name)==0) {
            sscanf(buff_line,"%s %d",dummy,&dispfreq);
            if (rank == MASTER)
                printf("dispfreq = %d\n",dispfreq);
        }
        strcpy(param_name,"store_range_to_finalstep");
        param_name_len = strlen(param_name);
        strncpy(dummy,buff_line,param_name_len);
        dummy[param_name_len] = '\0';
        if (strcmp(dummy,param_name)==0) {
            sscanf(buff_line,"%s %d",dummy,&store_range_to_finalstep);
            if (rank == MASTER)
                printf("storage_range_to_finalstep = %d\n",store_range_to_finalstep);
        }
        
//        else
//            sscanf(buf_line,"%d %d %d %d %lf %lf %d", &Nparticles, &Natoms_per_particle_ini, &Natoms_instant_trans, &maxKMCstep, &kT, &eb0, &dispfreq);
    }
//    fscanf(fp, "%d %d %d %lf %lf %d", &Nparticles, &Natoms_per_particle_ini, &maxKMCstep, &kT, &eb0, &dispfreq);
    fclose(fp);
    if (Natoms_instant_trans < 1) {
        Natoms_instant_trans = 1;
        if (rank == MASTER)
            printf("Issue: Natoms_instant_trans is < 1, adjusted to 1\n");
    }
    
    T_profile = (double*) calloc(maxKMCstep,sizeof(double));
    kT = (double*) calloc(maxKMCstep,sizeof(double));
    if (strlen(T_filename)==0) {
        printf("kT is obtained from T_ini and T_end.\n");
        tmp = (T_end-T_ini)*1.0/(maxKMCstep-1);
        for (i=0; i<maxKMCstep; i++) {
            *(T_profile+i)  = tmp*i+T_ini;
            *(kT+i) = 0.02585199/300 * (tmp*i+T_ini);
        }
    }
    else {
        printf("kT is obtained from existing profile.\n");
        load_double_mat(T_profile,maxKMCstep,1,T_filename);
        for (i=0; i<maxKMCstep; i++) {
            tmp = *(T_profile+i);
            *(kT+i) = 0.02585199/300 * tmp;
        }
    }

//    for (i=0; i<maxKMCstep; i++)
//        printf("T_profile[%d]=%f, kT[%d]=%f\n",i,T_profile[i],i,kT[i]);
    
//    printf("Nparticles=%d, Natoms_per_particle_ini=%d, Natoms_instant_trans=%d, maxKMCstep=%d, kT=%f, eb0=%f, dispfreq=%d\n",Nparticles, Natoms_per_particle_ini, Natoms_instant_trans, maxKMCstep, kT, eb0, dispfreq);

    
    int crystal_natoms_max;
    crystal_natoms_max = (int) fmin(10179, Nparticles);
    int n_crystal_type = 2;
    
    int clstr_natoms[Nparticles];
    int clstr_crystal_type[Nparticles];
    int *clstr_index;
    double *eb, *ef_tr;
    double *ef_atoms;
    int avg_natoms;

    double *pr, *pr_cumnorm;
    int *K;
    int Ksize_row;
    double lambda, rand1, rand2;
    int event_id, event_id_row;
    
    eb = (double*) calloc(n_crystal_type,sizeof(double));
    ef_tr = (double*) calloc((2*Nparticles*Natoms_per_particle_ini+1),sizeof(double));

    char datfilename[200]="ef_atoms.dat";
    ef_atoms = (double*) malloc(crystal_natoms_max*n_crystal_type*sizeof(double));
    load_double_mat(ef_atoms,crystal_natoms_max,n_crystal_type,datfilename);
//    for (i=0;i<crystal_natoms_max;i++)
//        printf("ef_atoms[%d,:] = %f, %f\n",i,*(ef_atoms+n_crystal_type*i),*(ef_atoms+n_crystal_type*i+1));
    
	char dirname[100];
    strcpy(dirname,"dir_result");
    struct stat st = {0};
    if (rank == MASTER && stat(dirname, &st) == -1) {
        mkdir(dirname,0777);
        printf("Created a new directory for result storage.\n");
    }
    
    char filename1[200];    //to store size distribution
    char filename2[200];    //to store time and event
    sprintf(filename2,"%s/event_time.dat",dirname);
    fp2 = fopen(filename2,"w");
    
    // initialization
    t = 0;
    for (i=0;i<Nparticles;i++) {
        clstr_natoms[i] = Natoms_per_particle_ini;
        clstr_crystal_type[i] = i%2;    //0: icosahedron, 1: wulff
//        clstr_crystal_type[i] = 1;    //0: icosahedron, 1: wulff
    }
    clstr_index = (int*) malloc(Nparticles*sizeof(int));
    K = (int*) malloc(Nparticles*(Nparticles-1)/2*2*sizeof(int));
    pr = (double*) malloc(Nparticles*(Nparticles-1)/2*n_crystal_type*sizeof(double));
    pr_cumnorm = (double*) malloc(Nparticles*(Nparticles-1)/2*n_crystal_type*sizeof(double));
//    pr = (double*) malloc(100*sizeof(double));
//    for (i=0;i<Nparticles*(Nparticles-1)/2*n_crystal_type;i++) {
//        *(pr+i) = 0.1;
//        printf("pr[%d] = %f\n",i,*(pr+i));
//    }
//    if (rank == MASTER)
//        printf("pr was initialized: size of %d by %d\n",Nparticles*(Nparticles-1)/2,n_crystal_type);
//    pr_norm = (double*) calloc(Nparticles*(Nparticles-1)/2, sizeof(double));
    
    // pr division for MPI
    pr_tmp_size = Nparticles*(Nparticles-1)/2*n_crystal_type;
    pr_tmp = (double*) malloc(pr_tmp_size*sizeof(double));
    pr_tmp_master = (double*) malloc(pr_tmp_size*sizeof(double));
//    pr_tmp = (double*) calloc(pr_tmp_size,sizeof(double));
//    pr_tmp_master = (double*) calloc(pr_tmp_size,sizeof(double));
    
    if (rank == MASTER) {
    
        for (KMCstep=0; KMCstep<maxKMCstep; KMCstep++) {
            for (i=0;i<Nparticles;i++)
                *(clstr_index+i) = -1;
            l = 0;
            ef = 0.0;
            for (i=0;i<Nparticles;i++)
                if (clstr_natoms[i] > 0) {
                    *(clstr_index+l) = i;
                    l++;
                    ef = ef + *(ef_atoms + n_crystal_type*(clstr_natoms[i]-1) + clstr_crystal_type[i]);
    //                printf("%d-bin is not empty: %d atoms, %d type, ef = %f\n",i,clstr_natoms[i],clstr_crystal_type[i],*(ef_atoms + n_crystal_type*(clstr_natoms[i]-1) + clstr_crystal_type[i]));
                }
            num_cluster = l;
            printf("%d step: ef = %.4e, %d cluster, t = %.4e\n",KMCstep,ef,num_cluster,t);
            if (num_cluster==1)
                break;
    //        for (i=0;i<Nparticles;i++)
    //            printf("clstr_index[%d]=%d\n",i,clstr_index[i]);
            Ksize_row = nchoose2(clstr_index,num_cluster,K);
            row_offset = (int) ceil(Ksize_row*1.0/numprocs);
            
    //        printf("Ksize_row = %d, numprocs = %d, row_offset = %d\n",Ksize_row,numprocs, row_offset);
    //        for (i=0;i<Ksize_row;i++)
    //            printf("Rank%d: K[%d,:] = [%d][%d]\n",rank,i, *(K+i*2+0), *(K+i*2+1));
    //        pr = (double*) realloc(pr, Ksize_row*n_crystal_type);
    //        pr_norm = realloc(pr_norm, Ksize_row*n_crystal_type);
            ctr = 1;
            mtype = FROM_MASTER;
            for (i=1;i<numprocs;i++) {
                MPI_Send(&Ksize_row, 1, MPI_INT, i, mtype, MPI_COMM_WORLD);
                MPI_Send(&row_offset, 1, MPI_INT, i, mtype, MPI_COMM_WORLD);
                MPI_Send(&clstr_natoms[0], Nparticles, MPI_INT, i, mtype, MPI_COMM_WORLD);
                MPI_Send(&clstr_crystal_type[0], Nparticles, MPI_INT, i, mtype, MPI_COMM_WORLD);
                //MPI_Send(&K[0], Ksize_row*2, MPI_INT, i, mtype, MPI_COMM_WORLD);
                MPI_Send(&K[i*row_offset*2], row_offset*2, MPI_INT, i, mtype, MPI_COMM_WORLD);
            }
            //printf("MPI_Send from MASTER conduceted\n");

            // pr calculation by MASTER
            row_ini = 0;
            if (Ksize_row < numprocs)
                row_end = Ksize_row;
            else
                row_end = row_offset;
            
            for (i=row_ini;i<row_end;i++) {
                p1 = *(K+i*2+0);
                p2 = *(K+i*2+1);
                n1 = clstr_natoms[p1];
                n2 = clstr_natoms[p2];
                if (n1+n2 > crystal_natoms_max) {
                    ctr = 0;
                    break;
                }  
                crystal1 = clstr_crystal_type[p1];
                crystal2 = clstr_crystal_type[p2];
    //            printf("K(%d,:)=%d,%d,n1=%d,n2=%d,cry1=%d,cry2=%d\n",i,p1,p2,n1,n2,crystal1,crystal2);
    //            cal_eb(eb,ef_tr,n1,n2,crystal1,crystal2,ef_atoms,n_crystal_type,eb0);
                cal_eb3(eb,ef_tr,n1,n2,crystal1,crystal2,ef_atoms,n_crystal_type,eb0,Natoms_instant_trans);
//                printf("eb=(%f, %f)\n",eb[0],eb[1]);
                for (j=0;j<n_crystal_type;j++) {
                    tmp = pow(cbrt(n1) + cbrt(n2),2) * sqrt(T_profile[KMCstep])*sqrt(n1+n2)/sqrt(n1*n2) * exp(-eb[j]/kT[KMCstep]);
                    *(pr+n_crystal_type*i+j) = tmp;
                }
            }
			// printf("pr was calculated on MASTER\n");
            // pr calculation result from WORKERs
            for (i=1;i<numprocs;i++) {
                mtype = FROM_WORKER;
                //MPI_Recv(&pr_tmp_master[0], pr_tmp_size, MPI_DOUBLE, i, mtype, MPI_COMM_WORLD, &status);
                MPI_Recv(&pr_tmp_master[i*row_offset*n_crystal_type], row_offset*n_crystal_type, MPI_DOUBLE, i, mtype, MPI_COMM_WORLD, &status);
                MPI_Recv(&ctr_tmp, 1, MPI_INT, i, mtype, MPI_COMM_WORLD, &status);
                if (Ksize_row >= numprocs) {
                    row_ini = i * row_offset;
                    if (i == numprocs-1)
                        row_end = Ksize_row;
                    else
                        row_end = row_ini + row_offset;
                    
                    for (k=row_ini;k<row_end;k++)
                        for (j=0;j<n_crystal_type;j++) {
                            *(pr+n_crystal_type*k+j) = *(pr_tmp_master+n_crystal_type*k+j);
                        }
                    ctr = (int) fmin(ctr,ctr_tmp);
                }
            }
//            printf("pr was received from WORKERS\n");
/*			for (i=0;i<Ksize_row;i++)
				printf("pr[%d,:] = [%f\t%f]\n",i,*(pr+n_crystal_type*i+0),*(pr+n_crystal_type*i+1));*/
            if (ctr == 0)
                break;

            *(pr_cumnorm+0) = *(pr+0);
			lambda = 0.0;
            for (i=1;i<Ksize_row*n_crystal_type;i++) {
                *(pr_cumnorm+i) = *(pr_cumnorm+i-1)+*(pr+i);
				lambda = lambda + *(pr+i);
			}
//            printf("check if cumsum was correct: %f vs. %f\n",*(pr_cumnorm+Ksize_row*n_crystal_type-1),lambda);
            rand1 = (double)rand() / RAND_MAX;
            event_id = -1;
            for (i=0;i<Ksize_row*n_crystal_type;i++)
                if (*(pr_cumnorm+i)/lambda >= rand1) {
                    event_id = i;
                    break;
                }
            
            // event proceeding
            event_id_row = (int) floor(event_id*1.0/n_crystal_type);
            cid1 = *(K+2*event_id_row+0);
            cid2 = *(K+2*event_id_row+1);
			if (cid1 < 0 || cid1 >= Nparticles) 
				printf("Error: cid1 = %d",cid1);
			if (cid2 < 0 || cid2 >= Nparticles)
            	printf("Error: cid1 = %d",cid2);
//			printf("event_id_row=%d,cid1=%d(%d,%d),cid2=%d(%d,%d)\n",event_id_row,cid1,clstr_natoms[cid1],clstr_crystal_type[cid1],cid2,clstr_natoms[cid2],clstr_crystal_type[cid2]);
            fprintf(fp2,"%d %d %d %d ",clstr_natoms[cid2],clstr_crystal_type[cid2], clstr_natoms[cid1],clstr_crystal_type[cid1]);
            clstr_natoms[cid1] = clstr_natoms[cid1] + clstr_natoms[cid2];
            clstr_natoms[cid2] = 0;
            clstr_crystal_type[cid1] = event_id - n_crystal_type*event_id_row;
            clstr_crystal_type[cid2] = -1;
    //        printf("rand1=%f, event_id=%d, event_id_row=%d, cid1=%d, cid2=%d, cid1->%d atoms, %d crystal-type\n",rand1,event_id,event_id_row,cid1,cid2,clstr_natoms[cid1],clstr_crystal_type[cid1]);
            rand2 = (double)rand() / RAND_MAX;
            t = t - log(rand2)/lambda;
            
            fprintf(fp2,"%d %d %f %f %d %d\n",clstr_natoms[cid1],clstr_crystal_type[cid1],t,rand1,cid1,cid2);
//            printf("kT=%f, %d %d %f %f %d %d\n",kT[KMCstep],clstr_natoms[cid1],clstr_crystal_type[cid1],t,rand1,cid1,cid2);

            if (KMCstep%dispfreq==0) {
                sprintf(filename1,"%s/on_the_fly_data%04d.dat",dirname,KMCstep/dispfreq);
                fp = fopen(filename1,"w");
                fprintf(fp,"%d %f\n",KMCstep,t);
                for (i=0;i<Nparticles;i++)
                    fprintf(fp,"%d %d\n",clstr_natoms[i],clstr_crystal_type[i]);
                fclose(fp);
            }
            
            if (maxKMCstep - KMCstep < store_range_to_finalstep) {
                sprintf(filename1,"%s/almost_final_data%04d.dat",dirname,KMCstep);
                fp = fopen(filename1,"w");
                fprintf(fp,"%d %f\n",KMCstep,t);
                for (i=0;i<Nparticles;i++)
                    fprintf(fp,"%d %d\n",clstr_natoms[i],clstr_crystal_type[i]);
                fclose(fp);
            }
            
            if (clstr_natoms[cid1] >= crystal_natoms_max)
                break;

        }
    }
    else {
        for (KMCstep=0; KMCstep<maxKMCstep; KMCstep++) {
            mtype = FROM_MASTER;
            MPI_Recv(&Ksize_row, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
            MPI_Recv(&row_offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
            MPI_Recv(&clstr_natoms[0], Nparticles, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
            MPI_Recv(&clstr_crystal_type[0], Nparticles, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
            //MPI_Recv(&K[0], Ksize_row*2, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
			MPI_Recv(&K[rank*row_offset*2], row_offset*2, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
			//printf("Rank%d: MPI_Recv was conducted\n",rank);
			//printf("Rank%d: Ksize_row=%d, rowoffset=%d\n",rank,Ksize_row,row_offset);
//			for (i=0;i<Ksize_row;i++)
//				printf("Rank%d: K[%d,:] = [%d][%d]\n",rank,i,*(K+2*i),*(K+2*i+1));
            if (Ksize_row >= numprocs) {
                row_ini = rank * row_offset;
                if (rank == numprocs-1)
                    row_end = Ksize_row;
                else
                    row_end = row_ini+row_offset;
                
                ctr_tmp = 1;
                
                for (i=row_ini;i<row_end;i++) {
                    p1 = *(K+i*2+0);
                    p2 = *(K+i*2+1);
                    n1 = clstr_natoms[p1];
                    n2 = clstr_natoms[p2];
                    if (n1+n2 > crystal_natoms_max) {
                        ctr_tmp = 0;
                        break;
                    }
                    crystal1 = clstr_crystal_type[p1];
                    crystal2 = clstr_crystal_type[p2];
                    //            printf("K(%d,:)=%d,%d,n1=%d,n2=%d,cry1=%d,cry2=%d\n",i,p1,p2,n1,n2,crystal1,crystal2);
                    //            cal_eb(eb,ef_tr,n1,n2,crystal1,crystal2,ef_atoms,n_crystal_type,eb0);
                    cal_eb3(eb,ef_tr,n1,n2,crystal1,crystal2,ef_atoms,n_crystal_type,eb0,Natoms_instant_trans);
//                    printf("Rank%d: eb=(%f, %f)\n",rank,eb[0],eb[1]);
                    for (j=0;j<n_crystal_type;j++) {
                        tmp = pow(cbrt(n1) + cbrt(n2),2) * sqrt(T_profile[KMCstep])*sqrt(n1+n2)/sqrt(n1*n2) * exp(-eb[j]/kT[KMCstep]);
                        *(pr_tmp+n_crystal_type*i+j) = tmp;
                        lambda = lambda + tmp;
//                        printf("Rank:%d,p.r.(%d,%d) = pr(%d) = %f, lambda = %f\n",rank,i,j,n_crystal_type*i+j,*(pr_tmp+n_crystal_type*i+j),lambda);
                    }
                }
            }
			//printf("pr was calculated WORKER %d\n",rank);
            mtype = FROM_WORKER;
            //MPI_Send(&pr_tmp[0], pr_tmp_size, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
            MPI_Send(&pr_tmp[rank*row_offset*n_crystal_type], row_offset*n_crystal_type, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
            MPI_Send(&ctr_tmp, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
//			printf("MPI_Send from WORKER %d\n",rank);
        }
    }
    //MPI_Finalize();
    
    sprintf(filename1,"%s/final_result.dat",dirname);
    fp = fopen(filename1,"w");
    fprintf(fp,"%d %f\n",KMCstep,t);
    for (i=0;i<Nparticles;i++)
        fprintf(fp,"%d %d\n",clstr_natoms[i],clstr_crystal_type[i]);
    fclose(fp);
    fclose(fp2);

	if (numprocs > 1)
		if (rank == MASTER) {
			printf("The job ended.  MPI_Abort-related errors might show up because of a few\n");
			printf("slave-workers still in loop. However, all outputs were produced correctly.\n");
			printf("--------------------------------------------------------------------------\n\n");
		}
	MPI_Abort(MPI_COMM_WORLD,1);
}

void load_double_mat(double *A, int Arow, int Acol, char *datfilename)
{
        
    FILE *fp;
    int i, j;
    double tmp;
        
    printf("filename: %s\n", datfilename);
    fp = fopen(datfilename, "r");
    for (i=0;i<Arow;i++)
        for (j=0;j<Acol;j++) {
            fscanf(fp, "%lf", &tmp);
            *(A+Acol*i+j) = tmp;
        }
    
    fclose(fp);
    
}

int nchoose2(int *clstr_index,int num_cluster, int *K) {
    
    int i, j, row_size;
    
    row_size = 0;
    for (i=0;i<num_cluster;i++)
        for (j=i+1;j<num_cluster;j++) {
            *(K+2*row_size+0) = *(clstr_index+i);
            *(K+2*row_size+1) = *(clstr_index+j);
//            printf("K(%d,:) = %d,%d\n",row_size,*(clstr_index+i),*(clstr_index+j));
            row_size++;
        }
//    printf("K matrix was produced\n");
    return(row_size);
}


        
        
        
        
        
        
        
        
        
        
        
        
