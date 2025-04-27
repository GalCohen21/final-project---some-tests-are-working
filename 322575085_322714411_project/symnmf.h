#ifndef SYMNMF_H
#define SYMNMF_H

/*
create a new matrix of size n x m
input: n (rows), m (columns)
output: pointer to created matrix
*/
double** allocate_matrix(int n, int m);

/*
Free the memory allocated for a matrix
input: matrix pointer, number of rows n
output: none
*/
void free_matrix(double** matrix, int n);

/*
Copy a matrix
input: A matrix, dim n, dim m
output: A copy of the matrix
*/
double** copy_matrix(double** src, int n, int m);

/*
calculate squared euclidean distance 
input: two points, dimention of the points
output: squared euclidean distance between point1, point2
*/
double euclidean_dist_squared(double* point1, double* point2, int dim);

/*
calculate similarity matrix
input: set of points, size of matrix, dim of points
output: similarity matrix
*/
double** calculate_similarity_matrix(double** data, int n, int dim);

/*
calculate diagonal matrix
input: similarity matrix, size of matrix
output: diagonal matrix
*/
double** calculate_diagonal_degree_matrix(double** similarity, int n);

/*
calculate normalized similarity matrix
input: similarity matrix, size of matrix
output: normalized similarity matrix
*/
double** calculate_normalized_similarity(double** A, int n);

/*
calculate Frobenius norm of difference between two matrices
input: 2 matrices , size n (rows), size m (cols)
output: Frobenius norm of (A-B)
*/
double frobenius_norm_dist(double** A, double** B, int n, int m);

/*
multiply matrices (n x m) and (m x l)
input: 2 matrices, sizes (n, m, l)
output: product matrix (n x l)
*/
double** matrix_mul(double** A, double** B, int n, int m, int l);

/*
transpose a matrix
input: a matrix, sizes (n, m)
output: the src matrix transposed (m x n)
*/
double** matrix_transpose(double** A, int n, int m);

/*
symmetric NMF iterative update
input: initial matrix H, normalized similarity W, size n x k
output: optimized matrix H
*/
double** symnmf(double** H, double** W, int n, int k);

#endif