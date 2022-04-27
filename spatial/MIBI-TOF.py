import scanpy as sc
import squidpy as sq
import numpy as np
import matplotlib.pyplot as plt

# load
adata = sq.datasets.mibitof()

# plot
for library_id in adata.uns["spatial"].keys():
    sc.pl.spatial(adata[adata.obs["library_id"] == library_id], color="Cluster", library_id=library_id, title=library_id)

# z-stack img
imgs = []
for library_id in adata.uns["spatial"].keys():
    img = sq.im.ImageContainer(adata.uns["spatial"][library_id]["images"]["hires"], library_id=library_id)
    img.add_img(adata.uns["spatial"][library_id]["images"]["segmentation"], library_id=library_id, layer="segmentation")
    img["segmentation"].attrs["segmentation"] = True
    imgs.append(img)
img = sq.im.ImageContainer.concat(imgs)

# show images
img.show("image")
img.show("image", segmentation_layer="segmentation",channelwise=True)

# convert image to CMYK
def rgb2cmyk(arr):
    R = arr[..., 0] / 255
    G = arr[..., 1] / 255
    B = arr[..., 2] / 255
    K = 1 - (np.max(arr, axis=-1) / 255)
    C = (1 - R - K) / (1 - K + np.finfo(float).eps)  # avoid division by 0
    M = (1 - G - K) / (1 - K + np.finfo(float).eps)
    Y = (1 - B - K) / (1 - K + np.finfo(float).eps)
    return np.stack([C, M, Y, K], axis=3)

img.apply(rgb2cmyk, layer="image", new_layer="image_cmyk", copy=False)
img.show("image_cmyk", channelwise=True)