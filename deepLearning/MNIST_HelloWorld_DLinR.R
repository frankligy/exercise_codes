if(!("keras" %in% installed.packages()[,"Package"])) devtools::install_github("rstudio/keras")
library(keras)
install_keras(method="conda")  # have to use conda, otherwise, it will create a virtualenv in ~/virtualenv/r-reticulate
                               # if using "conda", it will create a conda env named r-reticalate, and tensorflow and keras will all be installed there
library(reticulate)
reticulate::use_condaenv("r-reticulate",required = T)   # force to use r-reticulate env, required has to be TRUE
reticulate::use_python("/Users/ligk2e/opt/anaconda3/envs/r-reticulate/bin/python",required = T)
reticulate::py_config()


# Import the MNIST dataset
mnist <- dataset_mnist()

# Create variables for our test and training data
x_train <- mnist$train$x
y_train <- mnist$train$y
x_test <- mnist$test$x
y_test <- mnist$test$y

# Preprocess the data
# reshape
x_train <- array_reshape(x_train, c(nrow(x_train), 784))  # array_reshape is row-major reshaping, different from default R setting
x_test <- array_reshape(x_test, c(nrow(x_test), 784))
# rescale
x_train <- x_train / 255
x_test <- x_test / 255

# The y data is an integer vector with values ranging from 0 to 9. 
# To prepare this data for training we one-hot encode the vectors 
# into binary class matrices
y_train <- to_categorical(y_train, 10)
y_test <- to_categorical(y_test, 10)

# Build the model - 
# Building the neural network requires configuring the layers of the model, 
# then compiling the model.

# Setup the layers

model <- keras_model_sequential() 
model %>% 
  layer_dense(units = 256, activation = 'relu', input_shape = c(784)) %>% 
  layer_dropout(rate = 0.4) %>% 
  layer_dense(units = 128, activation = 'relu') %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 10, activation = 'softmax')

# Print the details of the model
model %>% summary()

# Compile the model - define loss and optimizer

model %>% compile(
  loss = 'categorical_crossentropy',
  optimizer = optimizer_rmsprop(),
  metrics = c('accuracy')
)

# Train the model

history <- model %>% fit(
  x_train, y_train, 
  epochs = 30, batch_size = 128, 
  validation_split = 0.2
)

# The history object returned by fit() includes loss and accuracy metrics 
# which we can plot
plot(history)

# Evaluate the model’s performance on the test data
model %>% evaluate(x_test, y_test)

# Generate predictions on new data:
model %>% predict_classes(x_test)