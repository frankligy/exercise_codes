#!/usr/bin/env python
# coding: utf-8

# ## Let's manipulate a real png image using python PIL library, we will implement both blurring and sharpening. This original imge is 200*200 pixel and it's a RGB image

# In[15]:


from PIL import Image, ImageDraw
from matplotlib.pyplot import show,imshow
import numpy as np


# In[2]:


input_image = Image.open("dog.png")   # PngImageFile, you can show, access size using this object
input_pixels = input_image.load()     # PixelAccess, you can access RGBalpha value using this object


# In[20]:


get_ipython().run_line_magic('matplotlib', 'inline')
imshow(np.asarray(input_image))


# ## First do blurring by using blur_kernel with size 3*3

# In[3]:


# Box Blur kernel
blur_kernel = [[1 / 9, 1 / 9, 1 / 9],
              [1 / 9, 1 / 9, 1 / 9],
              [1 / 9, 1 / 9, 1 / 9]]
# Box sharpening kernel
sharpening_kernel = [[  0  , -.5 ,    0 ],
          [-.5 ,   3  , -.5 ],
          [  0  , -.5 ,    0 ]]


# In[4]:


kernel = blur_kernel


# In[5]:


offset = len(kernel) // 2


# In[6]:


output_image = Image.new("RGB", input_image.size)   # Image
draw = ImageDraw.Draw(output_image)   # ImageDraw


# In[7]:


for x in range(offset, input_image.width - offset):
    for y in range(offset, input_image.height - offset):
        acc = [0, 0, 0]
        for a in range(len(kernel)):
            for b in range(len(kernel)):
                xn = x + a - offset
                yn = y + b - offset
                pixel = input_pixels[xn, yn]
                acc[0] += pixel[0] * kernel[a][b]
                acc[1] += pixel[1] * kernel[a][b]
                acc[2] += pixel[2] * kernel[a][b]

        draw.point((x, y), (int(acc[0]), int(acc[1]), int(acc[2])))


# In[21]:


get_ipython().run_line_magic('matplotlib', 'inline')
imshow(np.asarray(output_image))


# ## Then let's do sharpening by using sharpening_kernel of size 3*3

# In[22]:


kernel = sharpening_kernel


# In[23]:


offset = len(kernel) // 2
output_image = Image.new("RGB", input_image.size)   # Image
draw = ImageDraw.Draw(output_image)   # ImageDraw
for x in range(offset, input_image.width - offset):
    for y in range(offset, input_image.height - offset):
        acc = [0, 0, 0]
        for a in range(len(kernel)):
            for b in range(len(kernel)):
                xn = x + a - offset
                yn = y + b - offset
                pixel = input_pixels[xn, yn]
                acc[0] += pixel[0] * kernel[a][b]
                acc[1] += pixel[1] * kernel[a][b]
                acc[2] += pixel[2] * kernel[a][b]

        draw.point((x, y), (int(acc[0]), int(acc[1]), int(acc[2])))


# In[24]:


get_ipython().run_line_magic('matplotlib', 'inline')
imshow(np.asarray(output_image))


# In[ ]:




