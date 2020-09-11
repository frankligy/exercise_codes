#!/usr/bin/env python
# coding: utf-8

# ## Platform
# 
# 1. OSU super computer, Owen cluster, access through ondemand desktop.
# 2. Launching jupyter notebook and requests 28 gpu codes with latest cuda version
# 3. here I modified the code from rendering mandelbrot to julia set
# 4. excecution time comparion is at the bottom
# 5. the detailed implemetation is shown below

# ## Let's render Julia set using basic Python syntax

# In[11]:


import numpy as np
from pylab import imshow, show    # show plot from 2D array
from timeit import default_timer as timer    # execution time for a small snippet of code, otherwise use time.time()


# In[40]:


def julia(x, y, max_iters):

  z = complex(x, y)   # here 
  c = complex(-0.35, 0.65)
  for i in range(max_iters):
    z = z*z + c
    if abs(z) >= 4:
      return i

  return max_iters


# In[42]:


def create_fractal(min_x, max_x, min_y, max_y, image, iters):
  height = image.shape[0]
  width = image.shape[1]

  pixel_size_x = (max_x - min_x) / width
  pixel_size_y = (max_y - min_y) / height
    
  for x in range(width):
    real = min_x + x * pixel_size_x
    for y in range(height):
      imag = min_y + y * pixel_size_y
      color = julia(real, imag, iters)
      image[y, x] = color


# In[43]:


image = np.zeros((1024, 1536), dtype = np.uint8)
start = timer()
create_fractal(-2.0, 1.0, -1.0, 1.0, image, 20) 
dt = timer() - start

print("julia created in %f s" % dt)
imshow(image)
show()


# ## Let's render julia set using compiled python codes

# In[44]:


from numba import jit

@jit
def julia(x, y, max_iters):
  """
    Given the real and imaginary parts of a complex number,
    determine if it is a candidate for membership in the Mandelbrot
    set given a fixed number of iterations.
  """
  z = complex(x, y)
  c = complex(-0.35, 0.65)
  for i in range(max_iters):
    z = z*z + c
    if (z.real*z.real + z.imag*z.imag) >= 4:
      return i

  return max_iters

@jit
def create_fractal(min_x, max_x, min_y, max_y, image, iters):
  height = image.shape[0]
  width = image.shape[1]

  pixel_size_x = (max_x - min_x) / width
  pixel_size_y = (max_y - min_y) / height
    
  for x in range(width):
    real = min_x + x * pixel_size_x
    for y in range(height):
      imag = min_y + y * pixel_size_y
      color = julia(real, imag, iters)
      image[y, x] = color


# In[45]:


image = np.zeros((1024, 1536), dtype = np.uint8)
start = timer()
create_fractal(-2.0, 1.0, -1.0, 1.0, image, 20) 
dt = timer() - start

print("julia created in %f s" % dt)
imshow(image)
show()


# ## Let's render julia set using python CUDA on GPU

# In[46]:


from numba import cuda
from numba import *


# In[47]:


julia_gpu = cuda.jit(restype=uint32, argtypes=[f8, f8, uint32], device=True)(julia)


# In[48]:


@cuda.jit(argtypes=[f8, f8, f8, f8, uint8[:,:], uint32])
def julia_kernel(min_x, max_x, min_y, max_y, image, iters):
  height = image.shape[0]
  width = image.shape[1]

  pixel_size_x = (max_x - min_x) / width
  pixel_size_y = (max_y - min_y) / height

  startX, startY = cuda.grid(2)
  gridX = cuda.gridDim.x * cuda.blockDim.x;
  gridY = cuda.gridDim.y * cuda.blockDim.y;

  for x in range(startX, width, gridX):
    real = min_x + x * pixel_size_x
    for y in range(startY, height, gridY):
      imag = min_y + y * pixel_size_y 
      image[y, x] = julia_gpu(real, imag, iters)


# In[49]:


gimage = np.zeros((1024, 1536), dtype = np.uint8)
blockdim = (32, 8)
griddim = (32,16)

start = timer()
d_image = cuda.to_device(gimage)
mandel_kernel[griddim, blockdim](-5.0, 1.0, -4.0, 1.0, d_image, 20) 
d_image.to_host()
dt = timer() - start

print("julia created on GPU in %f s" % dt)

imshow(gimage)
show()


# 
# ## Conclusions
# 
# #### 1. basic syntax: 3.404387 s
# #### 2. compiled python: 0.250383 s
# #### 3. CUDA on GPU: 0.003787 s

# In[ ]:




