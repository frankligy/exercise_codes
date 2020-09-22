import os
import matplotlib.pyplot as plt
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import tensorflow_datasets as tfds   # pip install tesnforflow_datasets, have to downgrade absl-py==0.8

(ds_train,ds_test), ds_info = tfds.load(
    'mnist',
    split=['train','test'] ,
    shuffle_files=True,
    as_supervised=False,    # true will be (img,label), false will be a dict
    with_info=True
)

print(ds_info)
fig = tfds.show_examples(ds_train,ds_info,rows=4,cols=4)   # need as_supervised=False

def normalize_img(image,label):
    # normalize images
    return tf.cast(image,tf.float32)/255.0, label


df_train = ds.train.map(normalize_img)







