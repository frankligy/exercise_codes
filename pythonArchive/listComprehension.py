###### list comprehension
# example 1
M = [[1,2,3],[4,5,6],[7,8,9]] 
N = [[2,2,2],[3,3,3],[4,4,4]] 
​
[M[row][col]*N[row][col] for row in range(3) for col in range(3)] 
​
# [2, 4, 6, 12, 15, 18, 28, 32, 36] 

[[M[row][col]*N[row][col] for col in range(3)] for row in range(3)] 
# [[2, 4, 6], [12, 15, 18], [28, 32, 36]] 
​
[[M[row][col]*N[row][col] for row in range(3)] for col in range(3)] 
​
# [[2, 12, 28], [4, 15, 32], [6, 18, 36]]
# you can also range(len()) or in a direct list i.e. [1,2,3]

# example 2
bob = {'pay': 3000, 'job': 'dev', 'age': 42, 'name': 'Bob Smith'} 
​
sue = {'pay': 4000, 'job': 'hdw', 'age': 45, 'name': 'Sue Jones'} 
​
people = [bob, sue] 
​
[rec['age']+100 if rec['age'] >= 45 else rec['age'] for rec in people] 

# example 3
[item for item in list_ if item >= 1]   # compare with example 2, depending on whether you have else condition.

# example 4
condition = True
a = 1 if condition else 2    # just if statement

# example 5
dicNode = {node:cluster for cluster,nodes in enumerate(clusters) for node in nodes}




###### dict comprehension
# example 1
strings = ['import','is','with','if','file','exception']
D = {key: val for val,key in enumerate(strings)} # enumerater types

# {'exception': 5, 'file': 4, 'if': 3, 'import': 0, 'is': 1, 'with': 2}



#### set comprehension
strings = ['a','is','with','if','file','exception'] 

{len(s) for s in strings}  # only retain the unique value
# {1, 2, 4, 9}


#######  generator
# method 1
g = (x * x for x in range(10))   # g is a generator 
next(g)


# method 2
def fib(max):
    n, a, b = 0, 0, 1
    while n < max:
        yield b  # same as print, but lead the function to a generator
        a, b = b, a + b
        n = n + 1


f = fib(6)      # <generator object fib at 0x104feaaa0>
next(f)
next(f)
