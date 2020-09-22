'''
in order to set up the environment:
    conda create -n tensorflow python=3.6
    pip install tensorflow   (this is for cpu version i think)

in pycharm:
    firstly download the community version, install it.
    create new project, choose interpreter as conda env, the interpreter path should be like:
    /Users/ligk2e/opt/anaconda3/env/tensorflow/bin/python

    in scratch, create files, save as into project, then edit the python files. Done



# some keyboard shortcut
execute selection: Opt + shift + E
'''

import tensorflow as tf

# initialization of Tensor
x = tf.constant(4.0)
x = tf.constant([[1,2,3],[4,5,6]])
x = tf.ones([3,3])
x = tf.zeros([2,3])
x = tf.eye(3)
x = tf.random.normal([3,3],mean=0,stddev=1)
x = tf.random.uniform([3,3],minval=0,maxval=1)
x = tf.range(9)
x = tf.cast(x,dtype=tf.float64)   # type conversion


# mathematical operations
x = tf.constant([1,2,3])
y = tf.constant([9,8,7])
z = x + y
z = x - y
z = x / y
z = x * y

z = tf.reduce_sum(x*y,axis=0)
z = tf.tensordot(x,y,axes=1)

z = x ** 5

x = tf.random.normal([2,3])
y = tf.random.normal([3,4])
z = tf.matmul(x,y)
z = x @ y


# indexing

x = tf.constant([0,1,1,2,3,1,2,3])
print(x[::2])
print(x[::-1])

indices = tf.constant([0,3])
x_ind = tf.gather(x,indices)

# Reshaping
x = tf.range(9)
x = tf.reshape(x,(3,3))

x = tf.transpose(x,perm=[1,0])
