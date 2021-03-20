'use strict';
// in chrome, use cmd+opt+J to open console window
console.log('hello,world');

/*
primitive datatype:
1. number
2. string
3. boolean
4. undefined (variable that only declared or hasn't been declared)


complex datatype:
1. object (typeof will return "object")
2. array  (typeof will return "object")
3. null   (typeof will return "object")
4. function  (typeof will return "function")
*/

// typeof operator will return a string object with its type
// strict equal === which will not do type cast implicitly.
// js has block scope (if and for) in addtion to global scope and function scope, so you can use let in the block scope

/*
const, you can not change const primitive value, but you can change const complex type (only a const reference), but you can not reassign const complex type
const has to be declared and initialized, can not assigned later.
*/

//arrow function (differences of this)

//template literal

//spread operator

//rest operator

//destructuring

// import, export

// class, constructor, getter and setter

// HTML DOM, Event listener
var x = document.getElementById("main");
x = x.getElementsByTagName("p");
console.log(x[1].innerHTML);