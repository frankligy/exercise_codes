#include <iostream>
#include <cmath>

using namespace std;

/*
check if gcc or clang is installed, usually they are
gcc -c
clang --version

Then in VS code
install c++ extension and code runner. then you can run .cpp by just clicking the run button 

*/



classs Book {
    private:
        string note;  // only code in this class code block can access, so it's necessary to construct a setter and getter public function to communicate with private variable
    public:
        string title;
        string author;
        int pages;
        Book(string atitle){   // this constructor function will be called when we instantiate the object
            title = atitle; 
        }

        bool hasRead(){
            if(a > 1){
                return true;
            }
            return false;
        }
};



void sayHi(string name, int age);    // declare a function, define it later

int main()
{
    // data types
    char grade = 'A';
    string phrase = "Frank Li";
    int age = 50;
    float gpa = 4.0;
    double gpaD = 4.0;
    bool isMale = true;

    // working with string
    string phrase = "Frank Li";
    phrase[2] == 'B';     // mutable, different from python
    phrase.length();
    phrase.find("Li",3);   // Start searching from index 3
    phrase.substr(8,3)   // 3 means length, 8 means from index 8

    // working with number
    cout << 10 / 3;      // will get 3, different from python
    cout << pow(2,5);   // call cmath package
    cout << fmax(3,10);   // return bigger number

    // get users input
    int age;
    cout << "Enter your age: ";
    cin >> age;
    cout << "You are " << age << " years old";

    string name;
    cout << "Enter your name: ";
    getline(cin, name);
    cout << "Hello " << name;

    // array
    int luckyNums[20] = {4,8,15,16,23,42};    # specify length as 20
    luckyNum[1] = 89;    // mutable
    luckyNum[10] = 78;
    cout << luckyNums[0];

    // if statement
    bool isMale = true;
    bool isTall = true;
    if(isMale && isTall){
        cout << "You are a male";
    } else if(!isMale) {
        cout << "You are not a male";
    } else {
        cout << "hello";
    }

    // switch
    switch(x){
        case 0:
            y = "hi";
            break;
        case 1:
            y = "wow";
            break;
    }

    // while
    while(index <= 5){
        cout << index << endl;
        index++;
    }

    do{
        cout << index << endl;
        index++;
    }while(index <= 5);

    // for loop
    for(int i = 1; i <= 5; i++){
        cout << i << endl;
    }

    int nums = {1,2,3,4,5}
    for(int i = 0; i < 5; i++){
        cout << nums[i] << endl;
    }

    // 2D array
    int numberGrid[3][2] = {
        {1,2},
        {3,4},
        {5,6}
    };

    cout << numberGrid[0][1]

    // pointer
    double gpa = 2.7;
    cout << &gpa;   // will print out its memory address, a hexadecimal figure, aka, a pointer
    double *pgpa = &gpa // declare and define a pointer variable to contain this pointer/memory address
    cout << pgpa   // will print out memory address
    cout << *pgpa   // dereference, will print out the value stored in this memory address
    cout << *&gpa 

    // class and object
    Book book1("Harry Potter");
    book1.author = "JK Rowling";
    book1.page = 500;

    cout << book1.title;

    
    cout << "Hello world!" << endl;   // endl means printing a newline
    return 0;
}


void sayHi(string name, int age){
    cout << "Hello " << name << age
}