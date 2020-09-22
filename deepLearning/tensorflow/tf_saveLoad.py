import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers, regularizers
from tensorflow.keras.datasets import mnist

(x_train,y_train),(x_test,y_test) = mnist.load_data()
x_train = x_train.reshape(-1,28,28,1).astype('float32') / 255.0
x_test = x_test.reshape(-1,28,28,1).astype('float32') / 255.0

# imagine we have model and model.fit

# now we save the model weights
model.save_weights('save_model/')   # you will have a save_model folder


# first initialize model
# then load weights
model.load_weights('save_model')
# compile and fit




# Serialize the entire model
model.save('complete_saved_model')
# deserialize the entire model
model = keras.model.load_model('complete_saved_model')   # now you have the whole model, no need to compile, just need to fit
