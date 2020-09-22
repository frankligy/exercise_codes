import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.datasets import mnist
import tensorflow_hub as hub

(x_train,y_train),(x_test,y_test) = mnist.load_data()
x_train  = x_train.reshape(-1,28*28).astype('float32') / 255.0
x_test = x_test.reshape(-1,28*28).astype('float32') / 255.0

model = keras.models.load_model('pretrained/')

model.trainable = False
for layer in model.layers:
    assert layer.trainable = False
    layer.trainable = False

base_inputs = model.layers[0].input
base_outputs = model.layers[-2].output
final_outputs = layers.Dense(10)(base_outputs)

new_model = keras.Model(inputs=base_inputs,outputs=final_outputs)
print(new_model.summary())




model.compile(
    loss = keras.losses.SparseCategoricalCrossentropy(from_logits=False),  # if not run softmax, set it as true
    optimizer = keras.optimizers.Adam(lr=0.001),
    metrics=['accuracy'],
)

model.fit(x_train,y_train,batch_size=32,epochs=5,verbose=2)