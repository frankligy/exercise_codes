# Two resources to help to understand CNN
# 1. CNN.png in my onedrive learning folder
# 2. youtube, Andrew Ng deep learning course, learn what padding, stride is.

library(keras)
library(reticulate)
reticulate::use_condaenv("r-reticulate",required = T)   # force to use r-reticulate env, required has to be TRUE
reticulate::use_python("/Users/ligk2e/opt/anaconda3/envs/r-reticulate/bin/python",required = T)


# Import the MNIST dataset
mnist <- dataset_mnist()

# Create variables for our test and training data
x_train <- mnist$train$x
y_train <- mnist$train$y
x_test <- mnist$test$x
y_test <- mnist$test$y

# Preprocess the data
# reshape
# input image dimensions of MNIST data
img_rows = 28 
img_cols = 28
x_train <- array_reshape(x_train, c(nrow(x_train), img_rows, img_cols, 1))
x_test <- array_reshape(x_test, c(nrow(x_test), img_rows, img_cols, 1))

# rescale
x_train <- x_train / 255
x_test <- x_test / 255

# The y data is an integer vector with values ranging from 0 to 9. 
# To prepare this data for training we one-hot encode the vectors 
# into binary class matrices
y_train <- to_categorical(y_train)
y_test <- to_categorical(y_test)

# Build the model - 
# Building the neural network requires configuring the layers of the model, 
# then compiling the model.

# Setup the layers

model <- keras_model_sequential() 
model %>% 
  layer_conv_2d(filters = 32, kernel_size = c(3, 3), activation = "relu",
                input_shape = c(28, 28, 1)) %>%
  #project 3 : Add the two layers should be added here to run the full code!
  # layer_conv_2d(filters = 32, kernel_size = c(3,3), activation = "relu") %>%
  # layer_max_pooling_2d(pool_size = c(2,2)) %>%
  layer_dropout(rate = 0.4) %>% 
  layer_flatten() %>%
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