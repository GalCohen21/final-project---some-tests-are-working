#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "symnmf.h"

#define EPSILON 1e-4
#define MAX_ITER 300
#define BETA 0.5

/*
create a new matrix of size n x m
input: n (rows), m (columns)
output: pointer to created matrix
*/
double** allocate_matrix(int n, int m) {
    int i;
    double** matrix = (double**)malloc(n * sizeof(double*));
    if (matrix == NULL) {
        return NULL;
    }
    for (i = 0; i < n; i++) {
        matrix[i] = (double*)calloc(m, sizeof(double));
        if (matrix[i] == NULL) {
            while (--i >= 0) {
                free(matrix[i]);
            }
            free(matrix);
            return NULL;
        }
    }
    return matrix;
}

/*
Free the memory allocated for a matrix
input: matrix pointer, number of rows n
output: none
*/
void free_matrix(double** matrix, int n) {
    int i;
    if (matrix == NULL) {
        return;
    }
    for (i = 0; i < n; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

/*
Copy a matrix
input: A matrix, dim n, dim m
output: A copy of the matrix
*/
double** copy_matrix(double** src, int n, int m) {
    int i, j;
    double** dst = allocate_matrix(n, m);
    if (dst == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            dst[i][j] = src[i][j];
        }
    }
    return dst;
}

/*
calculate squared euclidean distance 
input: two points, dimention of the points
output: squared euclidean distance between point1, point2
*/
double euclidean_dist_squared(double* point1, double* point2, int dim) {
    double dist = 0.0;
    int i;
    for (i = 0; i < dim; i++) {
        dist += (point1[i] - point2[i]) * (point1[i] - point2[i]);
    }
    return dist;
}

/*
calculate similarity matrix
input: set of points, size of matrix, dim of points
output: similarity matrix
*/
double** calculate_similarity_matrix(double** data, int n, int dim) {
    double** A = allocate_matrix(n, n);
    int i, j;
    double dist;
    if (A == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i != j) {
                dist = euclidean_dist_squared(data[i], data[j], dim);
                A[i][j] = exp(-dist / 2);
            }
        }
    }
    return A;
}

/*
calculate diagonal matrix
input: similarity matrix, size of matrix
output: diagonal matrix
*/
double** calculate_diagonal_degree_matrix(double** similarity, int n) {
    double** D = allocate_matrix(n, n);
    int i, j;
    double sum;
    if (D == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i = 0; i < n; i++) {
        sum = 0.0;
        for (j = 0; j < n; j++) {
            sum += similarity[i][j];
        }
        D[i][i] = sum;
    }
    return D;
}

/*
calculate normalized similarity matrix
input: similarity matrix, size of matrix
output: normalized similarity matrix
*/
double** calculate_normalized_similarity(double** A, int n) {
    double** D = calculate_diagonal_degree_matrix(A, n);
    double** W = allocate_matrix(n, n);
    int i, j;
    if (W == NULL) {
        free_matrix(D, n);
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (D[i][i] > EPSILON && D[j][j] > EPSILON) {
                W[i][j] = A[i][j] / (sqrt(D[i][i]) * sqrt(D[j][j]));
            }
        }
    }
    free_matrix(D, n);
    return W;
}

/*
calculate Frobenius norm of difference between two matrices
input: 2 matrices , size n (rows), size m (cols)
output: Frobenius norm of (A-B)
*/
double frobenius_norm_dist(double** A, double** B, int n, int m) {
    int i, j;
    double sum = 0.0;
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            sum += (A[i][j] - B[i][j]) * (A[i][j] - B[i][j]);
        }
    }
    return sqrt(sum);
}

/*
multiply matrices (n x m) and (m x l)
input: 2 matrices, sizes (n, m, l)
output: product matrix (n x l)
*/
double** matrix_mul(double** A, double** B, int n, int m, int l) {
    int i, j, k;
    double** C = allocate_matrix(n, l);
    if (C == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < l; j++) {
            for (k = 0; k < m; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

/*
transpose a matrix
input: a matrix, sizes (n, m)
output: the src matrix transposed (m x n)
*/
double** matrix_transpose(double** A, int n, int m) {
    int i, j;
    double** B = allocate_matrix(m, n);
    if (B == NULL) {
        printf("An Error Has Occurred\n");
        exit(1);
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            B[j][i] = A[i][j];
        }
    }
    return B;
}

/* symmetric NMF iterative update
input: initial matrix H, normalized similarity W, size n x k
output: optimized matrix H */
double** symnmf(double** H, double** W, int n, int k) {
    int t, i, j;
    double diff;
    double numerator, denominator;
    double **WH, **Ht, **HHt, **HHtH;
    double** H_next = copy_matrix(H, n, k);
    for (t = 0; t < MAX_ITER; t++) {
        WH = matrix_mul(W, H, n, n, k);
        Ht = matrix_transpose(H, n, k);
        HHt = matrix_mul(H, Ht, n, k, n);
        HHtH = matrix_mul(HHt, H, n, n, k);
        for (i = 0; i < n; i++) {
            for (j = 0; j < k; j++) {
                numerator = WH[i][j];
                denominator = HHtH[i][j];
                if (denominator == 0) {
                    printf("An Error Has Occurred\n");
                    exit(1);
                }
                else{
                    H_next[i][j] = H[i][j] * (1 - BETA + (BETA * (numerator / denominator)));
                }
            }
        }
        diff = frobenius_norm_dist(H_next, H, n, k);
        free_matrix(WH, n);
        free_matrix(Ht, k);
        free_matrix(HHt, n);
        free_matrix(HHtH, n);
        for (i = 0; i < n; i++) {
            for (j = 0; j < k; j++) {
                H[i][j] = H_next[i][j];
            }
        }
        if (diff < EPSILON) {
            break;
        }
    }
    free_matrix(H_next, n);
    return H;
}