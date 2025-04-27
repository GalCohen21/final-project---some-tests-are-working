#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "symnmf.h"

int main(int argc, char* argv[]) {
    char *goal, *file_name;
    FILE *fp;
    int i, j, n = 0, d = 0, lines = 0;
    char buffer[1024];
    char *token;
    double **data, **result;
    
    /* Check argument count */
    if (argc != 3) {
        printf("An Error Has Occurred\n");
        return 1;
    }
    
    goal = argv[1];
    file_name = argv[2];
    
    /* Count number of points and dimensions */
    fp = fopen(file_name, "r");
    if (fp == NULL) {
        printf("An Error Has Occurred\n");
        return 1;
    }
    
    /* Count lines and dimensions */
    while (fgets(buffer, sizeof(buffer), fp)) {
        if (lines == 0) {
            token = strtok(buffer, ",");
            while (token != NULL) {
                d++;
                token = strtok(NULL, ",");
            }
        }
        lines++;
    }
    n = lines;
    
    /* Read the data */
    rewind(fp);
    data = allocate_matrix(n, d);
    if (data == NULL) {
        printf("An Error Has Occurred\n");
        fclose(fp);
        return 1;
    }
    
    for (i = 0; i < n; i++) {
        fgets(buffer, sizeof(buffer), fp);
        token = strtok(buffer, ",");
        j = 0;
        while (token != NULL && j < d) {
            data[i][j] = atof(token);
            token = strtok(NULL, ",");
            j++;
        }
    }
    fclose(fp);
    
    /* Process based on goal */
    if (strcmp(goal, "sym") == 0) {
        result = calculate_similarity_matrix(data, n, d);
    } else if (strcmp(goal, "ddg") == 0) {
        double **sim = calculate_similarity_matrix(data, n, d);
        result = calculate_diagonal_degree_matrix(sim, n);
        free_matrix(sim, n);
    } else if (strcmp(goal, "norm") == 0) {
        double **sim = calculate_similarity_matrix(data, n, d);
        result = calculate_normalized_similarity(sim, n);
        free_matrix(sim, n);
    } else {
        printf("An Error Has Occurred\n");
        free_matrix(data, n);
        return 1;
    }
    
    /* Print results */
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%.4f", result[i][j]);
            if (j < n - 1) {
                printf(",");
            }
        }
        printf("\n");
    }
    
    /* Free memory */
    free_matrix(data, n);
    free_matrix(result, n);
    
    return 0;
}