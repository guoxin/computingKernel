/**********************************************
*  This subroutine tries to compute the normalized
*  K3 kernel in parallel.
*  GK3WN for get K3 with weights and normalization.
*  
*  CK3: for constructing K3 from K3pep and K3alle.
*  
*  should be compiled with: 
*  	[OFFICIAL] Uniform so library Makefile 1.1.zip
**********************************************/
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
#include <string.h>
#include "gx_core.h"

const char **peptides;
int length;
double (*K3_messenger)[];
const double (*K1)[20];
double *weight;

int currentTask, totalTask;
const static int indexAminoAcid[26] =
    { 0, -1, 4, 3, 6, 13, 7, 8, 9, -1, 11, 10, 12, 2, -1, 14, 5, 1, 15, 16,
    -1, 19, 17, -1, 18, -1
};

double *sqroots;

void *getK3N_Single_Weight(void *Data);
void *getK3N_normalize(void *Data);
void *CK3_slave(void *Data);

const double (*K3pep_messenger)[], (*K3alle_messenger)[];
int len_pep, len_alle;
const int *pep_index, *alle_index;

void CK3(double (*p_K3)[], const int p_length[1],
	 const double (*p_K3pep)[], const int p_len_pep[1],
	 const int p_pep_index[],
	 const double (*p_K3alle)[], const int p_len_alle[1],
	 const int p_alle_index[])
{
    pthread_t thread[CORENUM];
    pthread_attr_t attr;
    int k;
    int locker[CORENUM];

    K3_messenger = p_K3;
    length = p_length[0];
    K3pep_messenger = p_K3pep;
    len_pep = p_len_pep[0];
    pep_index = p_pep_index;
    K3alle_messenger = p_K3alle;
    len_alle = p_len_alle[0];
    alle_index = p_alle_index;

    currentTask = 0;		/*indexing from zero. */
    totalTask = length;		/*indexing to "totalTask - 1". */
    gx_sfgd_info(pthread_attr_init(&attr), "attr init");
    gx_sfgd_info(pthread_attr_setdetachstate
		 (&attr, PTHREAD_CREATE_JOINABLE), "attr set joinable");

    for (k = 0; k < CORENUM; k++) {
	locker[k] = k;
	gx_sfgd_info(pthread_create(thread + k, &attr, &CK3_slave,
				    (void *) (locker + k)),
		     "creating, CK3");
    }
    for (k = 0; k < CORENUM; k++) {
	gx_sfgd_info(pthread_join(*(thread + k), NULL), "joining, CK3");
    }
}

void GK3WN(const char **p_peptides, const int p_length[1],
	   double (*p_K3)[], const double p_K1[20][20], double p_weight[])
{
    pthread_t thread[CORENUM];
    pthread_attr_t attr;
    int k;
    int locker[CORENUM];
    double roots[p_length[0]];
    double (*K3_local)[p_length[0]];

    peptides = p_peptides;
    length = p_length[0];
    K3_messenger = p_K3;
    K1 = p_K1;
    weight = p_weight;
    sqroots = roots;
    K3_local = p_K3;

    currentTask = 0;		/*indexing from zero. */
    totalTask = length;		/*indexing to "totalTask - 1". */
    gx_sfgd_info(pthread_attr_init(&attr), "attr init");
    gx_sfgd_info(pthread_attr_setdetachstate
		 (&attr, PTHREAD_CREATE_JOINABLE), "attr set joinable");

    for (k = 0; k < CORENUM; k++) {
	locker[k] = k;
	gx_sfgd_info(pthread_create
		     (thread + k, &attr, &getK3N_Single_Weight,
		      (void *) (locker + k)), "creating, getK3");
    }
    for (k = 0; k < CORENUM; k++) {
	gx_sfgd_info(pthread_join(*(thread + k), NULL), "joining, getK3");
    }

/*** the following parts do normalization. ***/
    for (k = 0; k < length; k++) {
	sqroots[k] = sqrt(K3_local[k][k]);
	K3_local[k][k] = 1;
    }

    for (k = 0; k < CORENUM; k++) {
	locker[k] = k;
	gx_sfgd_info(pthread_create(thread + k, &attr, &getK3N_normalize,
				    (void *) (locker + k)),
		     "creating, normalize");
    }
    for (k = 0; k < CORENUM; k++) {
	gx_sfgd_info(pthread_join(*(thread + k), NULL),
		     "joining, normalize");
    }
/*****delete this part simply cancels the normalization.****/
}

void *getK3N_normalize(void *taskList)
{
    int this_task;
    int s;
    double (*K3)[length] = K3_messenger;
    for (this_task = ((int *) taskList)[0];
	 this_task < length; this_task += CORENUM) {
	for (s = this_task + 1; s < length; s++) {
	    K3[this_task][s] /= (sqroots[s] * sqroots[this_task]);
	    K3[s][this_task] = K3[this_task][s];
	}
    }
    return NULL;
}

void *CK3_slave(void *taskList)
{
    int i, j;
    double (*K3)[length] = K3_messenger;
    const double (*K3pep)[len_pep] = K3pep_messenger;
    const double (*K3alle)[len_alle] = K3alle_messenger;
    for (i = ((int *) taskList)[0]; i < length; i += CORENUM) {
	for (j = i; j < length; j++) {
	    K3[i][j] = K3[j][i] =
		K3pep[pep_index[i]][pep_index[j]] *
		K3alle[alle_index[i]][alle_index[j]];
	}
    }
}

void *getK3N_Single_Weight(void *taskList)
{
    int i, j, k, s;
    double (*K3)[length] = K3_messenger;
    double sum, sum_temp;
    int this_task;
    int len_f, len_g, len_k;
    int row, col;

    for (this_task = ((int *) taskList)[0];
	 this_task < length; this_task += CORENUM) {
	for (s = this_task; s < length; s++) {
	    len_f = strlen(peptides[this_task]);
	    len_g = strlen(peptides[s]);
	    len_k = (len_f > len_g) ? len_g : len_f;

	    double slk[len_f][len_g], slk1[len_f][len_g];
	    sum = 0;
	    sum_temp = 0;
	    for (i = 0; i < len_f; i++) {
		for (j = 0; j < len_g; j++) {
		    slk[i][j] = slk1[i][j] =
			K1[indexAminoAcid
			   [(int) (peptides[this_task][i] - 'A')]]
			[indexAminoAcid[(int) (peptides[s][j] - 'A')]];
		    sum_temp += slk[i][j];
		}
	    }
	    sum += sum_temp * weight[0];
	    row = len_f;
	    col = len_g;
	    for (k = 1; k < len_k; k++) {
		sum_temp = 0;
		row--;
		col--;
		for (i = 0; i < row; i++) {
		    for (j = 0; j < col; j++) {
			slk[i][j] *= slk1[i + k][j + k];
			sum_temp += slk[i][j];
		    }
		}
		sum += sum_temp * weight[k];
	    }
	    K3[this_task][s] = K3[s][this_task] = sum;
	}
      /*****the main prog ends:***************/
    }
    return NULL;
}
