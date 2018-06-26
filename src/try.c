#include <libgen.h>
#include <string.h>
#include <stdio.h>
int main(void) {
  char buff[100];
  strcpy(buff, "/a/b/c/foo.txt");
  printf("basename: %s\n", basename(buff));
  return 0;
}
