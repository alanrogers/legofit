#include <stdio.h>
#include <string.h>
#include <errno.h>

int main(int argc, char **argv) {

    char buff[10];
    int status, err = ENXIO;

    status = strerror_r(err, buff, sizeof(buff));

    printf("status=%d\n", status);
    printf("buff=%s\n", buff);
    
	return 0;
}
