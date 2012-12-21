/* Document: my personal notes on Nov 22, 2012 */
/* Copyright(c) Xin Guo on Nov 22, 2012,
	released under the GPLv3 license */

#include <math.h>
#include <string.h>
#include <stdlib.h>
#define PIFACTOR 0.017453292519943295769236907684886127134428718885

const static int indexAminoAcid[26] =
    { 0, -1, 4, 3, 6, 13, 7, 8, 9, -1, 11, 10, 12, 2, -1, 14, 5, 1, 15, 16,
	-1, 19, 17, -1, 18, -1
};

double K12K3(double (*p_K1)[], int crow, int ccol, int maxLen)
{
	/* the code uses a full K1 matrix and will NOT pollute it */
	double (*K1)[ccol] = p_K1;
	int len, i, j;
	double sum, sum_temp;
	double (*slk)[ccol] =
	    (double (*)[ccol])malloc(crow * ccol * sizeof(double));
	if (crow < maxLen)
		maxLen = crow;
	if (ccol < maxLen)
		maxLen = ccol;
	sum = sum_temp = 0;
	for (i = 0; i < crow; i++) {
		for (j = 0; j < ccol; j++) {
			sum_temp += slk[i][j] = K1[i][j];
		}
	}
	sum += sum_temp;
	for (len = 1; len < maxLen; len++) {
		sum_temp = 0;
		crow--;
		ccol--;
		for (i = 0; i < crow; i++) {
			for (j = 0; j < ccol; j++) {
				slk[i][j] *= K1[i + len][j + len];
				sum_temp += slk[i][j];
			}
		}
		sum += sum_temp;
	}
	free(slk);
	return (sum);
}

void StrGenK3(const char **p_Fa, const char **p_Ga,
	      const double *p_FPHI, const double *p_GPHI,
	      const double *p_FPSY, const double *p_GPSY,
	      const int p_Flength[1], const int p_Glength[1],
	      const int p_maxLen[1], const double (*p_CSigma)[],
	      const double p_AAK1[20][20], double p_K3[1])
{
	int Flength = p_Flength[0], Glength = p_Glength[0];
	int maxLen = p_maxLen[0];
	int i, j, aaInd;
	int *f = (int *)malloc(Flength * sizeof(int));
	int *g = (int *)malloc(Glength * sizeof(int));
	double *DF = (double *)malloc(Flength * sizeof(double));
	double *DG = (double *)malloc(Glength * sizeof(double));
	double (*CSigma)[4] = (double (*)[4])p_CSigma, cPHI, cPSY, sPHI, sPSY;
	double (*K1)[Glength] =
	    (double (*)[Glength])malloc(Flength * Glength * sizeof(double));

	for (i = 0; i < Flength; i++) {
		f[i] = j = indexAminoAcid[(int)((*p_Fa)[i] - 'A')];
		cPHI = cos(p_FPHI[i] * PIFACTOR);
		cPSY = cos(p_FPSY[i] * PIFACTOR);
		sPHI = sin(p_FPHI[i] * PIFACTOR);
		sPSY = sin(p_FPSY[i] * PIFACTOR);
		DF[i] = 1 / (CSigma[j][0] * cPHI * cPSY +
			     CSigma[j][1] * cPHI * sPSY +
			     CSigma[j][2] * sPHI * cPSY +
			     CSigma[j][3] * sPHI * sPSY);
	}
	for (i = 0; i < Glength; i++) {
		g[i] = j = indexAminoAcid[(int)((*p_Ga)[i] - 'A')];
		cPHI = cos(p_GPHI[i] * PIFACTOR);
		cPSY = cos(p_GPSY[i] * PIFACTOR);
		sPHI = sin(p_GPHI[i] * PIFACTOR);
		sPSY = sin(p_GPSY[i] * PIFACTOR);
		DG[i] = 1 / (CSigma[j][0] * cPHI * cPSY +
			     CSigma[j][1] * cPHI * sPSY +
			     CSigma[j][2] * sPHI * cPSY +
			     CSigma[j][3] * sPHI * sPSY);
	}

	i = strlen(*p_Fa);
	maxLen = (maxLen > i) ? i : maxLen;
	i = strlen(*p_Ga);
	maxLen = (maxLen > i) ? i : maxLen;
	for (i = 0; i < Flength; i++) {
		for (j = 0; j < Glength; j++) {
			K1[i][j] = p_AAK1[f[i]][g[j]] *
			    cos((p_FPHI[i] - p_GPHI[j]) * PIFACTOR) *
			    cos((p_FPSY[i] -
				 p_GPSY[j]) * PIFACTOR) * DF[i] * DG[j];
		}
	}
	p_K3[0] = K12K3(K1, Flength, Glength, maxLen);
	free(DF);
	free(DG);
	free(f);
	free(g);
	free(K1);
}

void StrGenK3WithD(const char **p_Fa, const char **p_Ga,
		   const double *p_FPHI, const double *p_GPHI,
		   const double *p_FPSY, const double *p_GPSY,
		   const int p_Flength[1], const int p_Glength[1],
		   const double *p_FD, const double *p_GD,
		   const int p_maxLen[1],
		   const double p_AAK1[20][20], double p_K3[1])
{
	int Flength = p_Flength[0], Glength = p_Glength[0];
	int maxLen = p_maxLen[0];
	int i, j, aaInd;
	int *f = (int *)malloc(Flength * sizeof(int));
	int *g = (int *)malloc(Glength * sizeof(int));
	double (*K1)[Glength] =
	    (double (*)[Glength])malloc(Flength * Glength * sizeof(double));

	for (i = 0; i < Flength; i++) {
		f[i] = indexAminoAcid[(int)((*p_Fa)[i] - 'A')];
	}
	for (i = 0; i < Glength; i++) {
		g[i] = indexAminoAcid[(int)((*p_Ga)[i] - 'A')];
	}

	i = strlen(*p_Fa);
	maxLen = (maxLen > i) ? i : maxLen;
	i = strlen(*p_Ga);
	maxLen = (maxLen > i) ? i : maxLen;
	for (i = 0; i < Flength; i++) {
		for (j = 0; j < Glength; j++) {
			K1[i][j] = p_AAK1[f[i]][g[j]] *
			    exp(cos((p_FPHI[i] - p_GPHI[j]) * PIFACTOR) +
				cos((p_FPSY[i] -
				     p_GPSY[j]) * PIFACTOR)) / (p_FD[i] *
								p_GD[j]);
		}
	}
	p_K3[0] = K12K3(K1, Flength, Glength, maxLen);
	free(f);
	free(g);
	free(K1);
}

void AaGenK3(const char **p_Fa, const char **p_Ga,
	     const int p_Flength[1], const int p_Glength[1],
	     const int p_maxLen[1], const double p_AAK1[20][20], double p_K3[1])
{
	int Flength = p_Flength[0], Glength = p_Glength[0];
	int maxLen = p_maxLen[0];
	int i, j, aaInd;
	double (*K1)[Glength] =
	    (double (*)[Glength])malloc(Flength * Glength * sizeof(double));

	i = strlen(*p_Fa);
	maxLen = (maxLen > i) ? i : maxLen;
	i = strlen(*p_Ga);
	maxLen = (maxLen > i) ? i : maxLen;
	for (i = 0; i < Flength; i++) {
		for (j = 0; j < Glength; j++) {
			K1[i][j] = p_AAK1[indexAminoAcid[(int)((*p_Fa)[i] -
							       'A')]]
			    [indexAminoAcid[(int)((*p_Ga)[j] - 'A')]];
		}
	}
	p_K3[0] = K12K3(K1, Flength, Glength, maxLen);
	free(K1);
}
