#include <stdio.h>
#include <string.h>
#include <errno.h>
int main(void) {
    char buff[100], str[100]={'\0'}, dummy[100]={'\0'};
    int i, j, status;
    strcpy(buff, " 12 34 \n");
    status = sscanf(buff, "%d %d %s %s\n", &i, &j, str, dummy);
    errno=0;
    printf("status=%d i=%d j=%d str=%s dummy=%s\n",
           status, i, j, str, dummy);
    printf("errno=%d\n",errno);
    return 0;
}
