#pragma once
#include "defs.h"

int write_prim(FILE* fp) {
	size_t total_elements = N1 * N2 * N3 * NPRIM;
	double psave[N1][N2][N3][NPRIM];
	for (int i = 0; i < N1; i++)
		for (int j = 0; j < N2; j++)
			for (int k = 0; k < N3; k++)
				for (int m = 0; m < NPRIM; m++)
					psave[i][j][k][m] = prim[i][j][k][m];
	if (fwrite(&psave, sizeof(double), total_elements, fp) != total_elements) {
		return -1;
	}
	return 0;
}