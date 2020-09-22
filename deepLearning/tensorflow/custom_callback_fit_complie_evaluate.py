import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers


class Custom(keras.Model):    # behave like a wrapper of final model
    def __init__(self,model):
        super(Custom,self).__init__()
        self.model = model

    def compile(self,optimizer,loss):
        super(Custom,self).compile()
        self.optimizer = optimizer
        self.loss = loss

    def train_step(self,data):    # will get called when model.fit
        x,y = data

        with tf.GradientTape() as tape:
            y_pred = self.model(x,training=True)
            loss = self.loss(y,y_pred)

        training_vars = self.trainable_variables
        gradients = tape.gradient(loss,training_vars)    # calculate gradient for all trainable variable

        self.optimizer.apply_gradients(zip(gradients,training_vars))  # apply calculated gradients to all trainable variable
        acc_metric.update_state(y,y_pred)   # if wanna custom complie as well
        # self.complied_metrics.update_state(y,y_pred)  # if you wanna use complied_metrics
        return {'loss':loss,'accuracy':acc_metric.result()}

    def test_step(self,data):    # will get called when model.evaluate
        x,y = data    # unpack the testing data
        y_pred = self.model(x,training=False)   # self.model defined in __init__
        loss = self.loss(y,y_pred)              # self.loss defined in compile, otherwise will inherit from parent
        acc_metric.update_state(y,y_pred)       # keras.metrics object
        return {'loss':loss, 'accuracy':acc_metric.result()}








if __name__ == '__main__':

    acc_metric = keras.metrics.SparseCategoricalAccuracy(name='accuracy')

    training = Custom(model)
    training.complie(
        optimizer = keras.optimizers.Adam(),
        loss = keras.losses.SparseCategoricalCrossentropy(from_logits=True),
    )
    training.fit(x_train,y_train,batch_size=32,epochs=2)
    training.evaluate(x_test,y_test,batch_size=32)


    for epoch in range(num_epochs):
        for batch_idx,(x_train,y_train) in enumerate(training):
            with tf.GradientTape() as tape:
                y_pred = model(x_train,training=True)
                loss = loss_fn(y_train,x_train)
            gradients = tape.gradient(loss,model.trainable_weights)
            optimizer.apply_gradients(zip(gradients,model.trainable_weights))
            acc_metric.update_state(y_train,x_train)

        train_acc = acc_metric.result()
        print('Accuracy over epoch {}'.format(train_acc))
        acc_metric.reset_states()

    for batch_idx,(x_test,y_test) in enumerate(testing):
        y_pred = model(x_test,training=False)
        acc_metric.update_state(y_test,y_pred)

    testing_acc = acc_metric.result()
    print('Accuracy over epoch {}'.format(testing_acc))
    acc_metric.reset_states()