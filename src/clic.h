#ifndef LEGO_CLIC
#define LEGO_CLIC

void matrix_mult(double** m1, int m1_rows, int m1_cols, double** m2, int m2_rows, int m2_cols, double** results);

double** get_fit_param_array(char* title, int num_files, int num_params);
double** get_fit_param_array_num_unkown(char* title, int num_params);

double** make_covar_matrix(double** array, int files, int params);

double KL_to_lnL(double KL, double* p_matrix, int p_matrix_size, double sum);

#endif
