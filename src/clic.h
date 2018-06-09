#ifndef LEGO_CLIC
#define LEGO_CLIC

void matrix_mult(double** m1, int m1_rows, int m1_cols, double** m2, int m2_rows, int m2_cols, double** results);

double KL_to_lnL(double KL, double* p_matrix, int p_matrix_size, double sum);

#endif
