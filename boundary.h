#pragma once

#include "defs.h"
#include "metric.h"
#include "mks.h"

void bound_prim(double prim[N1M][N2M][N3M][NPRIM])
{
	double gcov_orig[4][4];
	double gcon_orig[4][4];
	double gdet_orig;
	double gcov_dest[4][4];
	double gcon_dest[4][4];
	double gdet_dest;
	// 1. 鬼化X1方向
	for (int i = -1; i >= - NG; i--)
	{
		for (int j = 0; j < N2; j++)
		{
			for (int k = 0; k < N3; k++)
			{
				get_metric(i, j, k, gcov_orig, gcon_orig, gdet_orig);
				get_metric(i + 1, j, k, gcov_dest, gcon_dest, gdet_dest);
				prim[i][j][k][RHO] = prim[i + 1][j][k][RHO] * sqrt(fabs(gdet_dest)) / sqrt(fabs(gdet_orig));
				prim[i][j][k][UU] = prim[i + 1][j][k][UU] * sqrt(fabs(gdet_dest)) / sqrt(fabs(gdet_orig));
				prim[i][j][k][B1] = prim[i + 1][j][k][B1] * sqrt(fabs(gdet_dest)) / sqrt(fabs(gdet_orig));

				prim[i][j][k][U2] = prim[i + 1][j][k][U2] * (1 - dx1);
				prim[i][j][k][U3] = prim[i + 1][j][k][U3] * (1 - dx1);
				prim[i][j][k][B2] = prim[i + 1][j][k][B2] * (1 - dx1);
				prim[i][j][k][B3] = prim[i + 1][j][k][B3] * (1 - dx1);

				prim[i][j][k][U1] = prim[i + 1][j][k][U1] * (1 + dx1);
			}
		}
	}

	for (int i = N1; i < NG + N1; i++)
	{
		for (int j = 0; j < N2; j++)
		{
			for (int k = 0; k < N3; k++)
			{
				get_metric(i - 1, j, k, gcov_orig, gcon_orig, gdet_orig);
				get_metric(i, j, k, gcov_dest, gcon_dest, gdet_dest);
				prim[i][j][k][RHO] = prim[i - 1][j][k][RHO] * sqrt(fabs(gdet_dest)) / sqrt(fabs(gdet_orig));
				prim[i][j][k][UU] = prim[i - 1][j][k][UU] * sqrt(fabs(gdet_dest)) / sqrt(fabs(gdet_orig));
				prim[i][j][k][B1] = prim[i - 1][j][k][B1] * sqrt(fabs(gdet_dest)) / sqrt(fabs(gdet_orig));

				prim[i][j][k][U2] = prim[i - 1][j][k][U2] * (1 - dx1);
				prim[i][j][k][U3] = prim[i - 1][j][k][U3] * (1 - dx1);
				prim[i][j][k][B2] = prim[i - 1][j][k][B2] * (1 - dx1);
				prim[i][j][k][B3] = prim[i - 1][j][k][B3] * (1 - dx1);

				prim[i][j][k][U1] = prim[i - 1][j][k][U1] * (1 + dx1);
			}
		}
	}

	// 2. 鬼化X2方向（反射边界条件）
	for (int i = 0; i < N1; i++)
	{
		for (int j = -NG; j < 0; j++)
		{
			for (int k = 0; k < N3; k++)
			{
				// 1) 把上面的格子移动到下面
				prim[i][j][k][RHO] = prim[i][-1 - j][N3 - k - 1][RHO];
				prim[i][j][k][UU] = prim[i][-1 - j][N3 - k - 1][UU];
				prim[i][j][k][U1] = prim[i][-1 - j][N3 - k - 1][U1];
				prim[i][j][k][U2] = -prim[i][-1 - j][N3 - k - 1][U2];
				prim[i][j][k][U3] = prim[i][-1 - j][N3 - k - 1][U3];
				prim[i][j][k][B1] = prim[i][-1 - j][N3 - k - 1][B1];
				prim[i][j][k][B2] = -prim[i][-1 - j][N3 - k - 1][B2];
				prim[i][j][k][B3] = prim[i][-1 - j][N3 - k - 1][B3];
				// 2) 把下面的格子移动到上面
				prim[i][j + N2 + NG][k][RHO] = prim[i][N2 + NG - 1 + j][N3 - k - 1][RHO];
				prim[i][j + N2 + NG][k][UU] = prim[i][N2 + NG - 1 + j][N3 - k - 1][UU];
				prim[i][j + N2 + NG][k][U1] = prim[i][N2 + NG - 1 + j][N3 - k - 1][U1];
				prim[i][j + N2 + NG][k][U2] = -prim[i][N2 + NG - 1 + j][N3 - k - 1][U2];
				prim[i][j + N2 + NG][k][U3] = prim[i][N2 + NG - 1 + j][N3 - k - 1][U3];
				prim[i][j + N2 + NG][k][B1] = prim[i][N2 + NG - 1 + j][N3 - k - 1][B1];
				prim[i][j + N2 + NG][k][B2] = -prim[i][N2 + NG - 1 + j][N3 - k - 1][B2];
				prim[i][j + N2 + NG][k][B3] = prim[i][N2 + NG - 1 + j][N3 - k - 1][B3];
			}
		}
	}

	// 3. 鬼化X3方向（周期边界条件）
	for (int i = 0; i < N1; i++)
	{
		for (int j = 0; j < N2; j++)
		{
			for (int k = -NG; k < 0; k++)
			{
				// 1) 把上面的格子移动到下面
				prim[i][j][k][RHO] = prim[i][j][N3 + k][RHO];
				prim[i][j][k][UU] = prim[i][j][N3 + k][UU];
				prim[i][j][k][U1] = prim[i][j][N3 + k][U1];
				prim[i][j][k][U2] = prim[i][j][N3 + k][U2];
				prim[i][j][k][U3] = prim[i][j][N3 + k][U3];
				prim[i][j][k][B1] = prim[i][j][N3 + k][B1];
				prim[i][j][k][B2] = prim[i][j][N3 + k][B2];
				prim[i][j][k][B3] = prim[i][j][N3 + k][B3];
				// 2) 把下面的格子移动到上面
				prim[i][j][N3 + NG + k][RHO] = prim[i][j][k + NG][RHO];
				prim[i][j][N3 + NG + k][UU] = prim[i][j][k + NG][UU];
				prim[i][j][N3 + NG + k][U1] = prim[i][j][k + NG][U1];
				prim[i][j][N3 + NG + k][U2] = prim[i][j][k + NG][U2];
				prim[i][j][N3 + NG + k][U3] = prim[i][j][k + NG][U3];
				prim[i][j][N3 + NG + k][B1] = prim[i][j][k + NG][B1];
				prim[i][j][N3 + NG + k][B2] = prim[i][j][k + NG][B2];
				prim[i][j][N3 + NG + k][B3] = prim[i][j][k + NG][B3];
			}
		}
	}
}