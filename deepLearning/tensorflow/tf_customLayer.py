import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers, regularizers
from tensorflow.keras.datasets import mnist

(x_train,y_train),(x_test,y_test) = mnist.load_data()
x_train = x_train.reshape(-1,28,28,1).astype('float32') / 255.0
x_test = x_test.reshape(-1,28,28,1).astype('float32') / 255.0


class Dense(layers.Layer):
    def __init__(self,units,input_dim):
        super(Dense,self).__init__()
        self.w = self.add_weight(
            name='w',
            shape=(input_dim,units),
            initializer = 'random_normal',
            trainable=True,
        )
        self. = self.add_weight(
            name='b',
            shape=(units,),
            initializer='zeros',
            trainable=True
        )

    def call(self,inputs):
        return tf.matmul(inputs,self.w) + self.b



class MyModel(keras.Model):
    def __init__(self,num_classes=10):
        super(MyModel,self).__init__()
        self.dense1 = Dense(64)    # now we use user-defined layer
        self.dense2 = Dense(num_classes)

    def call(self,input_tensor):
        x = tf.nn.relu(self.dense1(input_tensor))
        return self.dense2(x)


model = MyModel()
model.compile(
    loss=keras.losses.SparseCategoricalCrossentropy(from_logits=True),
    optimizer=keras.optimizers.Adam(lr=0.001),
    metrics=['accuracy']
)
print(model.summary())
model.fit(x_train,y_train,batch_size=32,epochs=10,verbose=2)
model.evaluate(x_test,y_test,batch_size=32,verbose=2)