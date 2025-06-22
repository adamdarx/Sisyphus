#pragma once
#include "defs.h"
#include "mks.h"
#include "la.h"

void get_metric(int i, int j, int k, double gcov[4][4], double gcon[4][4], double &gdet)
{
	// 计算协变度规
	for(int row = 0; row < 4; row++)
		for (int col = 0; col < 4; col++)
			gcov[row][col] = metricFunc[row][col](0, X1min + i * dx1, X2min + j * dx2, X3min + k * dx3);
	// 计算逆变度规
	inverseMatrix(gcov, gcon);
	gdet = determinant(gcov);
}

void get_conn(int i, int j, int k, double conn[4][4][4])
{
	const double epsilon = 1e-3;
	double tmp[4][4][4];
	double gcovh[4][4], gconh[4][4], gdeth;
	double gcovl[4][4], gconl[4][4], gdetl;
	double gcov[4][4], gcon[4][4], gdet;
	double dx[4] = { 1, dx1, dx2, dx3 };

	for (int kk = 0; kk < 4; kk++) {
		int idel = kk == 1 ? 1 : 0;
		int jdel = kk == 2 ? 1 : 0;
		int kdel = kk == 3 ? 1 : 0;
		get_metric(i + idel * epsilon, j + jdel * epsilon, k + kdel * epsilon, gcovh, gconh, gdeth);
		get_metric(i - idel * epsilon, j - jdel * epsilon, k - kdel * epsilon, gcovl, gconl, gdetl);
		for (int ii = 0; ii < 4; ii++)
			for (int jj = 0; jj < 4; jj++)
				conn[ii][jj][kk] = (gcovh[ii][jj] - gcovl[ii][jj]) / (2 * epsilon * dx[kk]);
	}

	/* now rearrange to find \Gamma_{ijk} */
	for (int ii = 0; ii < 4; ii++)
		for (int jj = 0; jj < 4; jj++)
			for (int kk = 0; kk < 4; kk++)
				tmp[ii][jj][kk] = 0.5 * (conn[jj][ii][kk] + conn[kk][ii][jj] - conn[kk][jj][ii]);

	get_metric(i, j, k, gcov, gcon, gdet);
	/* finally, raise index */
	for (int ii = 0; ii < 4; ii++)
		for (int jj = 0; jj < 4; jj++)
			for (int kk = 0; kk < 4; kk++) {
				conn[ii][jj][kk] = 0.;
				for (int ll = 0; ll < 4; ll++) conn[ii][jj][kk] += gcon[ii][ll] * tmp[ll][jj][kk];
			}
}
