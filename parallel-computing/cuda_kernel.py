'''
CUDA kernel confer much more flexibility

'''

from numba import cuda

@cuda.jit          # no need to specify type, will automatically infer
def add_kernel(x,y,out):
    tx = cuda.threadIdx.x    # we have a 1d block with length 128
    ty = cuda.blockIdx.x        # we have a 1d grid with length 30
    block_size = cuda.blockDim.x   # 128
    grid_size = cuda.gridDim.x     # 30

    start = tx + ty * block_size    # for each thread, the starting point differ
    stride = block_size * grid_size

    for i in range(start,x.shape[0],stride):    # x.shape[0] = 100000
        out[i] = x[i] + y[i]  


import numpy as np

n = 100000
x = np.arange(n).astype(np.float32)
y = 2*x
out = np.empty_like(x)

threads_per_block = 128
blocks_per_grid = 30

add_kernel[blocks_per_grid,threads_per_block](x,y,out)
print(out[:10])


# simplify a little
@cuda.jit
def add_kernel(x,y,out):
    start = cuda.grid(1)      # automatically infer the starting point for each thread
    stride = cuda.gridsize(1)   # automatically infer the stride
    for i in range(start,x.shape[0],stride):
        out[i] = x[i] + y[i]



# synchronization
# if you use x_device, y_device, it will be asybchronous
# in order to explicitly implement synchronization
cuda.synchronize()
add_kernel[blocks_per_grid,threads_per_block](x_device,y_device,out_device)
cuda.synchronize()



# race condition
# atomic operations, sort of like lock, just combine 3 step as 1 indivisible step

@cuda.jit
def thread_counter_safe(global_conter):
    cuda.atomic.add(global_counter,0,1)