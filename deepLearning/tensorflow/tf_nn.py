import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.datasets import mnist

# # if you use gpu
# physical_device = tf.config.list_physical_devices('GPU')
# tf.config.experimental.set_memory_growth(physical_device[0],True)

(x_train,y_train),(x_test,y_test) = mnist.load_data()
# x_train: 60000,28,28
# y_train: 60000,
# x_test: 10000,28,28
# y_test: 10000,
# they are all ndarray dtype: np.float64

x_train  = x_train.reshape(-1,28*28).astype('float32') / 255.0
x_test = x_test.reshape(-1,28*28).astype('float32') / 255.0

x_train = tf.convert_to_tensor(x_train)   # optional

# Sequential API (very convenient, not very flexible)
model = keras.Sequential(
    [
        keras.Input(shape=(28*28)),
        layers.Dense(512,activation='relu'),
        layers.Dense(256,activation='relu'),
        layers.Dense(10),
    ]
)

print(model.summary())


## another way to build model one layer at once, so you can inspect the summary in real time as you change the model
model = keras.Sequential()
model.add(keras.Input(shape=(784)))
model.add(layers.Dense(512,activation='relu'))
print(model.summary())
model.add(layers.Dense(256,activation='relu'))
model.add(layers.Dense(10))


# Functional API (a bit more flexible)
inputs = keras.Input(shape=(784))
x = layers.Dense(512,activation='relu')(inputs)
x = layers.Dense(256,activation='relu')(x)
outputs = layers.Dense(10,activation='softmax')(x)
model = keras.Model(inputs=inputs,outputs=outputs)

# extract intermediate result from model
model = keras.Model(inputs=model.inputs,outputs=[layer.output for layer in model.layers])
features = model.predict(x_train)
for feature in features:
    print(feature.shape)
# above extract all layers, you can also specity outputs as a certain layer by use index or layer name
model = keras.Model(inputs=model.inputs,outputs=[model.layers[-2].output])
model = keras.Model(inputs=model.inputs,outputs=[model.get_layer('my_layer').output])



model.compile(
    loss = keras.losses.SparseCategoricalCrossentropy(from_logits=False),  # if not run softmax, set it as true
    optimizer = keras.optimizers.Adam(lr=0.001),
    metrics=['accuracy'],
)

model.fit(x_train,y_train,batch_size=32,epochs=5,verbose=2)
model.evaluate(x_test,y_test,batch_size=32,verbose=2)