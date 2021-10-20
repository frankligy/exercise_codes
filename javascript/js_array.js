
const fruits = ['banana','orange','apple'];


// access and length
console.log(fruits[0]);
console.log(fruits[fruits.length-1]);

// iterate (two way, one using numbered index, another just use raw item value)
for (let i=0; i<fruits.length; i++) {
    console.log('the',i,'one is:',fruits[i]);
}

fruits.forEach(atomFunc);

function atomFunc(value) {
    console.log('it is:',value);
}

// type determination
console.log(typeof fruits);
console.log(fruits instanceof Array);

// set
fruits[1] = 'lemon';
console.log(fruits);

// push
x = fruits.push('mango');
console.log(fruits,x);

// pop
y = fruits.pop();
console.log(fruits,y);

// shift
z = fruits.shift();
console.log(fruits,z);

// unshirt
a = fruits.unshift('kiwi');
console.log(fruits,a);

// splice

// concat, can also take two other arrays as argument, or just one string to merge
const newFruits = ['coconut','strawberry'];
const concatArray = fruits.concat(newFruits);
console.log(concatArray);

// slice
const nowFruits = concatArray.slice(1,3); // not including the end one
console.log(nowFruits);

// to string
console.log(fruits.toString());
console.log(fruits.join('*'));

// sort
concatArray.sort()
console.log(concatArray);
concatArray.reverse()
console.log(concatArray);

// if you want to compare number
const points = [40, 100, 1, 5, 25, 10];
points.sort(function(a, b){return a - b});
console.log(points);