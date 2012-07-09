/****************************************************************************
*This library tries to meet all the demand on computing the string kernel, as 
* well as the normalized ones. This library uses pthread and thus should be 
* compiled on a POSIX environment. The library is designed to be called by R.
*
* Copyright: Guo Xin on July 07, 2012. Released under GNU GPLv3 license.
*
* To call the library in R, use format like (the order should be followed):
* len1 <- length(string1.set); len2 <- length(string2.set)
* maxNchar <- max(c(nchar(string1.set), nchar(string2.set)))
* .C("K3",
*   peptides1 = as.character(string1.set), 
*   peptides2 = as.character(string2.set), 
*   Length1 = as.integer(len1),
*   Length2 = as.integer(len2),
*   maxNchar = maxNchar,
*   K3 = double(length = len1 * len2), 
*   K1 = K1 ^ beta, #the beta must be handled in R!!!
*   weights = as.double(rep(1, maxNchar)),
*   computing.mode = as.integer(0),
*   DUP = TRUE, #because of the characters
*   NAOK = FALSE,
*   PACKAGE = "K3")$K3
* Here comes some explanation:
* computing.mode: add the two categories.
*    category ONE: normalize (0) or not (1)
*    category TWO: symmetric (0) or unsymmetric (2) or paired (4)
* 1. symmetric means peptides2 and Length2 would be ignored and just
*    considered same as peptides1 and Length1. However in this case the former
*    ones should still be listed as parameters.
* 2. unsymmetric means the K3 might not be a square matrix and in any case
*    peptides1 and peptides2 should be considered as they are different.
* 3. paired means one wants just to compute some special K3 (or K3_hat).
*    In this case the Length1 should be equal to Length2 (otherwise, the c code
*    would give an unpredictable result) and the code just compute
*    K3(peptides1[i],peptides2[i]) for each i.
*
* Particularly, we require in this program that as long as illegal amino acid 
* characters appear, warnings must be given at the very beginning. (just 
* warnings, the program will continue running).
*******************************************************************************/
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
#include <string.h>
#include "gx_core.h"
/* WARNING: BE CAUTIOUS WHEN USING THE FOLLOWING TWO MACROS AS MAX(i++, j++) !!*/
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MIN(x, y) (((x) < (y)) ? (x) : (y))

const static int indexAminoAcid[26] =
    { 0, -1, 4, 3, 6, 13, 7, 8, 9, -1, 11, 10, 12, 2, -1, 14, 5, 1, 15, 16,
	-1, 19, 17, -1, 18, -1
};

typedef struct {
	int (*f)[];
	int (*g)[];
	const int *f_nchar;
	const int *g_nchar;
	const double (*K1)[20];
	const double *weight;
	int listLength;		/* f_len */
	int secondStepSize;	/* g_len */
	void *K3_store;
	int maxNchar;
} TASK;

double K3_single(const int *f, const int *g, const int f_nchar,
		 const int g_nchar, const double K1[20][20],
		 const double *weight)
{
	double slk[f_nchar][g_nchar], slk1[f_nchar][g_nchar];
	int i, j, k, s, row, col, min_len;
	double sum, sum_temp;
	min_len = MIN(f_nchar, g_nchar);

	sum = sum_temp = 0;
	for (i = 0; i < f_nchar; i++) {
		for (j = 0; j < g_nchar; j++) {
			sum_temp += (slk[i][j] = slk1[i][j] = K1[f[i]][g[j]]);
		}
	}
	sum += sum_temp * weight[0];
	row = f_nchar;
	col = g_nchar;
	for (k = 1; k < min_len; k++) {
		sum_temp = 0;
		row--;
		col--;
		for (i = 0; i < row; i++) {
			for (j = 0; j < col; j++) {
				slk[i][j] *= slk1[i + k][j + k];
				sum_temp += slk[i][j];
			}
		}
		sum += (sum_temp * weight[k]);
	}
	return (sum);
}

void *K3_pairSlave(void *task)
{
	int i, j;
	TASK thisTask = ((TASK *) task)[0];
	int (*f)[thisTask.maxNchar];
	int (*g)[thisTask.maxNchar];
	f = thisTask.f;
	g = thisTask.g;
	int listLength = thisTask.listLength;
	int secondStepSize = thisTask.secondStepSize;
	double *K3_store = (double *)thisTask.K3_store;
	j = 0;
	for (i = 0; i < listLength; i++) {
		K3_store[i] = K3_single(f[i], g[j],
					thisTask.f_nchar[i],
					thisTask.g_nchar[j],
					thisTask.K1, thisTask.weight);
		j += secondStepSize;
	}
}

void char2int(const char **ff, const int len, const int *f_nchar,
	      int (*store)[], const int maxNchar)
{
	int (*st)[maxNchar];
	int i, j;
	int theInteger;
	st = store;
	for (i = 0; i < len; i++) {
		for (j = 0; j < f_nchar[i]; j++) {
			theInteger = (int)(ff[i][j] - 'A');
			if ((theInteger > 25) || (theInteger < 0)) {
				error
				    ("Bad amino acid. Please use Capital Characters");
			}
			st[i][j] = indexAminoAcid[theInteger];
			if ((st[i][j] > 19) || (st[i][j] < 0)) {
				error
				    ("Unkown amino acid. Please follow BLOSUM62-2 index.");
			}
		}
	}
}

void K3_pair(const char **ff, const char **gg, const int len,
	     const double K1[20][20], const double *weight, double *K3_store,
	     const int normalize, const int maxNchar)
{
	int i, stepSize;
	int f_nchar[len], g_nchar[len];
	int startPt;
	pthread_t thread[CORENUM];
	pthread_attr_t attr;
	stepSize = len / CORENUM + 1;
	TASK task[CORENUM];

	for (i = 0; i < len; i++) {
		f_nchar[i] = strlen(ff[i]);
		g_nchar[i] = strlen(gg[i]);
		if((g_nchar[i]<=0)||(f_nchar[i]<=0)){
			error("empty string found.");
		}
	}
	int f[len][maxNchar], g[len][maxNchar];
	char2int(ff, len, f_nchar, f, maxNchar);
	char2int(gg, len, g_nchar, g, maxNchar);

	gx_sfgd_info(pthread_attr_init(&attr), "attr init");
	gx_sfgd_info(pthread_attr_setdetachstate
		     (&attr, PTHREAD_CREATE_JOINABLE), "attr set joinable");

	for (i = 0; i < CORENUM; i++) {
		startPt = i * stepSize;
		task[i].maxNchar = maxNchar;
		task[i].listLength = len - startPt;
		task[i].listLength = MIN(task[i].listLength, stepSize);
		task[i].f = f + startPt;
		task[i].g = g + startPt;
		task[i].f_nchar = f_nchar + startPt;
		task[i].g_nchar = g_nchar + startPt;
		task[i].K1 = K1;
		task[i].weight = weight;
		task[i].secondStepSize = 1;
		task[i].K3_store = (void *)(K3_store + startPt);
		gx_sfgd_info(pthread_create(thread + i, &attr, &K3_pairSlave,
					    (void *)(task + i)),
			     "creating, task");
	}
	for (i = 0; i < CORENUM; i++) {
		gx_sfgd_info(pthread_join(*(thread + i), NULL),
			     "joining, normalize");
	}
	if (normalize == 0) {
		double diag_f[len], diag_g[len];
		K3_pair(ff, ff, len, K1, weight, diag_f, 1, maxNchar);
		K3_pair(gg, gg, len, K1, weight, diag_g, 1, maxNchar);
		for (i = 0; i < len; i++) {
			K3_store[i] /= sqrt(diag_f[i] * diag_g[i]);
		}
	}
}

void *K3_unsymmetricSlave(void *task)
{
	int i, j;
	TASK thisTask = ((TASK *) task)[0], temp;
	int (*f)[thisTask.maxNchar] = thisTask.f;
	int (*g)[thisTask.maxNchar] = thisTask.g;
	int f_len = thisTask.listLength;
	int g_len = thisTask.secondStepSize;
	double (*K3_store)[f_len] = thisTask.K3_store;

	for (i = 0; i < g_len; i += CORENUM) {
		temp = thisTask;
		temp.g = g + i;
		temp.g_nchar = thisTask.g_nchar + i;
		temp.secondStepSize = 0;
		temp.K3_store = K3_store + i;
		(*K3_pairSlave) (&temp);
	}
}

void K3_unsymmetric(const char **ff, const char **gg,
		    const int f_len, const int g_len,
		    const double K1[20][20], const double *weight,
		    double (*K3)[], const int normalize, const int maxNchar)
{
	/*we abuse TASK type a little bit. */
	int i, j;
	int f_nchar[f_len], g_nchar[g_len];
	pthread_t thread[CORENUM];
	pthread_attr_t attr;
	TASK task[CORENUM];
	double (*K3_store)[f_len] = K3;

	for (i = 0; i < f_len; i++) {
		f_nchar[i] = strlen(ff[i]);
		if(f_nchar[i]<=0){
			error("empty string found.");
		}
	}
	for (i = 0; i < g_len; i++) {
		g_nchar[i] = strlen(gg[i]);
		if(g_nchar[i]<=0){
			error("empty string found.");
		}
	}
	int f[f_len][maxNchar], g[g_len][maxNchar];
	char2int(ff, f_len, f_nchar, f, maxNchar);
	char2int(gg, g_len, g_nchar, g, maxNchar);

	gx_sfgd_info(pthread_attr_init(&attr), "attr init");
	gx_sfgd_info(pthread_attr_setdetachstate
		     (&attr, PTHREAD_CREATE_JOINABLE), "attr set joinable");

	for (i = 0; i < CORENUM; i++) {
		task[i].maxNchar = maxNchar;
		task[i].listLength = f_len;
		task[i].secondStepSize = g_len - i;
		task[i].f = f;
		task[i].g = g + i;
		task[i].f_nchar = f_nchar;
		task[i].g_nchar = g_nchar + i;
		task[i].K1 = K1;
		task[i].weight = weight;
		task[i].K3_store = (void *)(K3_store + i);
		gx_sfgd_info(pthread_create
			     (thread + i, &attr, &K3_unsymmetricSlave,
			      (void *)(task + i)), "creating, task");
	}
	for (i = 0; i < CORENUM; i++) {
		gx_sfgd_info(pthread_join(*(thread + i), NULL),
			     "joining, normalize");
	}
	if (normalize == 0) {
		double diag_f[f_len], diag_g[g_len];
		K3_pair(ff, ff, f_len, K1, weight, diag_f, 1, maxNchar);
		K3_pair(gg, gg, g_len, K1, weight, diag_g, 1, maxNchar);
		for (j = 0; j < f_len; j++) {
			diag_f[j] = 1 / sqrt(diag_f[j]);
		}
		for (i = 0; i < g_len; i++) {
			diag_g[i] = 1 / sqrt(diag_g[i]);
		}
		for (i = 0; i < g_len; i++) {
			for (j = 0; j < f_len; j++) {
				K3_store[i][j] *= (diag_f[j] * diag_g[i]);
			}
		}
	}
}

void K3(const char **p_peptidesF, const char **p_peptidesG,
	const int p_lenF[1], const int p_lenG[1], const int maxNchar[1],
	double (*p_K3)[], const double p_K1[20][20], const double p_weight[],
	const int mode[1])
{
	int normalize = mode[0] % 2;
	int format = mode[0] / 2;
	if (format == 3) {
		double *K3_store = (double *)(p_K3);
		if (p_lenF[0] != p_lenG[0]) {
			error
			    ("error: peptide lists of different lengths in pair mode.");
		}
		K3_pair(p_peptidesF, p_peptidesG, p_lenF[0], p_K1, p_weight,
			K3_store, normalize, maxNchar[0]);
	} else if (format == 2) {
		double (*K3_store)[] = p_K3;
		K3_unsymmetric(p_peptidesF, p_peptidesG, p_lenF[0], p_lenG[0],
			       p_K1, p_weight, K3_store, normalize,
			       maxNchar[0]);
	} else {
		error("Unkown mode");
	}
}
