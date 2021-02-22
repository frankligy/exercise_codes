#include <stdio.h>
#include <stdlib.h>

typedef char* sm;

void test_func(char** arg)
{
    char result = **arg;
    printf("%c\n",result);
}


int main()
{
char a = 'A';
// sm pa = &a;
// printf("%p\n",pa);
// printf("%p",&pa);

char* pa = &a;
printf("%p\n",pa);
printf("%p\n",&pa);

test_func(&pa);

char buffer[10];
buffer[0] = (1 % 10) + '0';
printf("%s",buffer);

}



