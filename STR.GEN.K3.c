/* Document: my personal notes on Nov 22, 2012 */
/* Copyright(c) Xin Guo on Nov 22, 2012,
	released under the GPLv3 license */

#include <math.h>
#include <string.h>
#include <stdlib.h>

const static int indexAminoAcid[26] =
    { 0, -1, 4, 3, 6, 13, 7, 8, 9, -1, 11, 10, 12, 2, -1, 14, 5, 1, 15, 16,
	-1, 19, 17, -1, 18, -1
};

double K12K3(double (*p_K1)[], int crow, int ccol, int maxLen)
{
	/* the code uses a full K1 matrix and will NOT pollute it */
	int len, i, j;
	double (*K1)[ccol] = p_K1;
	double sum, sum_temp;
	double (*slk)[ccol] =
	    (double (*)[ccol])malloc(crow * ccol * sizeof(double));
	if (crow < maxLen)
		maxLen = crow;
	if (ccol < maxLen)
		maxLen = ccol;
	for (i = 0; i < crow; i++) {
		for (j = 0; j < ccol; j++) {
			slk[i][j] = 1.0;
		}
	}
	sum = 0;
	for (len = 0; len < maxLen; len++) {
		sum_temp = 0;
		for (i = 0; i < crow; i++) {
			for (j = 0; j < ccol; j++) {
				slk[i][j] *= K1[i + len][j + len];
				sum_temp += slk[i][j];
			}
		}
		sum += sum_temp;
		crow--;
		ccol--;
	}
	free(slk);
	return (sum);
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
