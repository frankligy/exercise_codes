
import squidpy as sq
import scanpy as sc
import numpy as np
import pandas as pd

arr = np.ones((100, 100, 3))
arr[40:60, 40:60] = [0, 0.7, 1]
img = sq.im.ImageContainer(arr, dims=('y','x','channels'),layer="img1")

arr_seg = np.zeros((100, 100))
arr_seg[40:60, 40:60] = 1
img.add_img(arr_seg, layer="seg1")

img.show(layer='img1')
