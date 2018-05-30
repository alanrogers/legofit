#ifndef LEGO_CLIC
#define LEGO_CLIC

void matrix_mult(double** m1, int m1_rows, int m1_cols, double** m2, int m2_rows, int m2_cols, double** results);

char* itoa(int data, char* buffer);
char* tri_cat(char* a, char* b, char* c, char* str);

double** get_fit_param_array(char* title, double** array, int num_files, int num_params);
double** get_fit_param_array_num_unkown(char* title, double** array, int num_params);

#endif
