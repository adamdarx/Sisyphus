#pragma once

double dot(double a[4], double b[4])
{
    return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3]);
}

double dot(double rowVec[3], double colVec[3], double metric[3][3]) {
    double vec[3] = { 0 , 0, 0 };
    double res = 0;
    for (int col = 0; col < 3; col++)
        for (int row = 0; row < 3; row++)
            vec[col] += rowVec[row] * metric[row][col];
    for (int index = 0; index < 3; index++)
        res += vec[index] * colVec[index];
    return res;
}

double square(double vec[3], double metric[3][3]) {
    return dot(vec, vec, metric);
}

double contract(double A[3][3], double B[3][3]) {
    double sum = 0;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            sum += A[i][j] * B[i][j];
    return sum;
}

double determinant(double matrix[4][4]) {
    return matrix[0][0] * (matrix[1][1] * (matrix[2][2] * matrix[3][3] - matrix[2][3] * matrix[3][2]) -
        matrix[1][2] * (matrix[2][1] * matrix[3][3] - matrix[2][3] * matrix[3][1]) +
        matrix[1][3] * (matrix[2][1] * matrix[3][2] - matrix[2][2] * matrix[3][1]))
        - matrix[0][1] * (matrix[1][0] * (matrix[2][2] * matrix[3][3] - matrix[2][3] * matrix[3][2]) -
            matrix[1][2] * (matrix[2][0] * matrix[3][3] - matrix[2][3] * matrix[3][0]) +
            matrix[1][3] * (matrix[2][0] * matrix[3][2] - matrix[2][2] * matrix[3][0]))
        + matrix[0][2] * (matrix[1][0] * (matrix[2][1] * matrix[3][3] - matrix[2][3] * matrix[3][1]) -
            matrix[1][1] * (matrix[2][0] * matrix[3][3] - matrix[2][3] * matrix[3][0]) +
            matrix[1][3] * (matrix[2][0] * matrix[3][1] - matrix[2][1] * matrix[3][0]))
        - matrix[0][3] * (matrix[1][0] * (matrix[2][1] * matrix[3][2] - matrix[2][2] * matrix[3][1]) -
            matrix[1][1] * (matrix[2][0] * matrix[3][2] - matrix[2][2] * matrix[3][0]) +
            matrix[1][2] * (matrix[2][0] * matrix[3][1] - matrix[2][1] * matrix[3][0]));
}

int inverseMatrix(double matrix[4][4], double inverse[4][4], int n = 4) {
    double augmented[4][4 * 2];
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            augmented[i][j] = matrix[i][j];
            augmented[i][j + n] = (i == j) ? 1.0 : 0.0;
        }
    }

    for (int i = 0; i < n; i++) {
        if (augmented[i][i] == 0) {
            int swapRow = -1;
            for (int j = i + 1; j < n; j++) {
                if (augmented[j][i] != 0) {
                    swapRow = j;
                    break;
                }
            }
            if (swapRow == -1) {
                return 0;
            }
            for (int j = 0; j < 2 * n; j++) {
                double temp = augmented[i][j];
                augmented[i][j] = augmented[swapRow][j];
                augmented[swapRow][j] = temp;
            }
        }

        double pivot = augmented[i][i];
        for (int j = 0; j < 2 * n; j++)
            augmented[i][j] /= pivot;

        for (int j = 0; j < n; j++)
            if (i != j) {
                double factor = augmented[j][i];
                for (int k = 0; k < 2 * n; k++)
                    augmented[j][k] -= factor * augmented[i][k];
            }
    }

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            inverse[i][j] = augmented[i][j + n];
    return 1;
}
