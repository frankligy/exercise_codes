#include <stdio.h>
#include <stdlib.h>


struct Student {
    int age;
    double gpa;
    char major[20];
    char name[20];
}

int main()
{
    char letter = 'A';  // char needs to be single quotation
    char name[] = "Frank";   // in C, no string type, only have char type, so need to declare a string array.
    printf("my name is %s\n", name);

    const int num = 5; // this number is immutable.

    /*
    int age;
    printf("input your age\n");
    scanf("%d", &age);
    printf("my age is %d",age);
    */

    /*
    char all_names[10];
    printf("input your name ok?\n");
    fgets(all_names,10,stdin);
    printf("my name is %s",all_names);
    */

    char all_names[10];
    printf("input your name ok?\n");
    scanf("%s",all_names);
    printf("my name is %s",all_names);

    struct Student student1;
    student1.age=15;
    student1.gpa=3.4;
    strcpy(student1.major,"Biology");
    strcpy(student1.name,"Frank");

    // pointer
    int age = 5;
    printf("%p",&age);  // ampersand means the physical address of a variable
    int * pAge = &age;  // think bout pAge is a pointer variable, pointer type is determinied by *, ignore int, that just means the type of variable you are storing
    printf("%p",pAge);
    printf("%d",*pAge);  // dereference

    //files
    FILE * fpointer = fopen("test.txt","w");  // the mode can be r,w,a  // fpointer points to/contains address of a new created file
    fprintf(fpointer,"hi,I am here");
    fclose(fpointer);

    char line[255];
    FILE * fpointer = fopen("test.txt","r");
    fgets(line,255,fpointer);  // get first line, and move fpointer to next line
    fgets(line,255,fpointer); // get second line, and move fpointer to next line
    printf(line);

    return 0;
}