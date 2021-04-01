#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>




int main() {

    // file
    FILE *fp;
    fp = fopen("/Users/ligk2e/Desktop/test.txt","w");
    fprintf(fp,"Hello");
    fseek(fp,0,SEEK_END);
    printf("%ld\n",ftell(fp));

    fp = fopen("/Users/ligk2e/Desktop/test.bin","w");
    char *mem = calloc (10,1);
    fwrite(mem,1,10,fp);
    free(mem);

    fp = fopen("/Users/ligk2e/Desktop/test.bin","r");
    char *newMem = malloc (10);
    fread(newMem,1,10,fp);
    

    // pointer
    // all the pointer are of size 8, size_t object is similar to int
    printf("%zu\n",sizeof(fp));
    char *str = "hey";   // string literral is read-only
    printf("%p\n",str);
    printf("%s\n",str);
    printf("%c\n",str[0]);
    printf("%p\n",&str);

    char *proxyStr = (char *) str;  // safeway to do
    printf("%s\n",proxyStr);  

    int a = 5;
    int *pA = &a;
    printf("%d\n",*pA);  // int pointer has to be dereferenced

    int *arr = malloc (sizeof(int)*3);   // int array, or any type array you name it
    arr[0] = 10;
    printf("%d\n",arr[0]);
    printf("%p\n",arr);   // arr name refer to the address of the first element, same for string literal, or any other type array.
    







    
    return 0;

    

}