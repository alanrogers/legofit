#ifndef LEGO_CLIC
#define LEGO_CLIC

void matrix_mult(float** m1, int m1_rows, int m1_cols, float** m2, int m2_rows, int m2_cols, float** results);

char* itoa(int data, char* buffer);
char* tri_cat(char* a, char* b, char* c, char* str);

float** get_fit_param_array(char* title, int num_files);
float** get_fit_param_array_num_unkown(char* title)

float** get_fit_param_array_num_unkown(char* title);
float** get_fit_param_array(char* title, int num_files);

#endif
