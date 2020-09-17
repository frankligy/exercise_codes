'''
numba:

1. function complier
2. type-specializing
3. just-in-time: compile when the function get called, so jit will know the input type
4. numerically-focused: not compatible if you have string or dict like that

'''

'''
Using numba or running on GPU not always good due to the fact that it will add a lot overhead when copy data to GPU an
if your input data is small, you are even not able to keep GPU busy.

'''


'''
4 layer of speeding up python code:

1. use pandas, numpy, scipy etc whose implemention is written by C
2. multi-processing or multi-threading on CPU
3. using numba to complie to machine code (limited to numerics)
4. on GPU (limited to numerics)

'''

# just complie to machine code using jit

from numba import jit
import math            # use scalar function in math module when using numba

@jit
def hypot(x,y):
    x = abs(x);
    y = abs(y);
    t = min(x,y);
    x = max(x,y);
    t = t/x;
    return x * math.sqrt(1+t*t)


hypot(3.0,4.0)
hypot.py_func(3.0,4.0)    # original python function



# run on GPU
from numba import vectorize

@vectorize(['int64(int64,int64)'],target='cuda')
def add_ufunc(x,y):
    return x+y


from numba import cuda

@cuda.jit(device=True)     # this function is just a normal python function, cuda.jit is different from numba.jit on CPU
def polar_to_cartesian(rho,theta):
    x = rho * math.cos(theta)(
    y = rho * math.sin(theta)
    return x,y

@vectorize(['float32(float32,float32,float32,float32)'],target='cuda')
def polar_distance(rho1,theta1,rho2,theta2):
    x1,y1 = polar_to_cartesian(rho1,theta1)
    x2,y2 = polar_to_cartesian(rho2,theta2)

    return ((x1-x2)**2 + (y1-y2)**2) ** 0.5


# memory management
# in the above example, we create array on CPU and need to copy to cuda and then return back to cpu

@vectorize(['int64(int64,int64)'],target='cuda')
def add_ufunc(x,y):
    return x+y

n = 100000
x = np.arange(n).astype(np.float32)
y = 2 * x

# avoid allocating from CPU to GPU by converting them to devicearray beforehand
x_device = cuda.to_device(x)
y_device = cuda.to_device(y)

print(x_device)    # numba.cuda.cudadrv.devicearray.DeviceNDArray object
print(x_device.shape)   # (100000,)
print(x_device.dtype)   # float32

# avoid allocatin from GPU to CPU by creating a output buffer
out_device = cuda.device_array(shape=(n,),dtype=np.float32)
add_ufunc(x_device,y_device,out=out_device)

out_device.copy_to_host()    # convert back to normal ndarray








