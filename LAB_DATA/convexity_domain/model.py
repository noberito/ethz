# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.2
# ---

import tensorflow as tf
import numpy as np
import pandas as pd
from keras import callbacks
from keras.models import Sequential
from keras.layers import Dense, Input

n_test = 40000
k = 100
good_data = True
gseed = 6
degree = 4

def custom_activation(x):
    return 0.5 * (1 + tf.math.sign(x)) * (x + 1/100) +  0.5 * (1 - tf.math.sign(x)) * tf.math.exp(tf.math.minimum(0.0,x)) / 100

def custom2_activation(x):
    return tf.where(x < 0, tf.exp(tf.math.minimum(0.0, x) * k) / k, x + 1/k)

if good_data:
    X = np.load("X_convex_{}.npy".format(degree))
    Y = np.load("Y_convex_{}.npy".format(degree))
else :
    X = np.random.randint(0, 100, size=(n_test * 2, 22))
    Y = np.random.randint(0, 2, size=(n_test * 2))

print(X)

X_mean = np.mean(X, axis = 0)
print(X_mean)
X_std = np.std(X, axis = 0)
print(X_std)
X_normalized = (X - X_mean) / X_std

X_train = X_normalized[:n_test,:22]
print(X_train)
Y_train = Y[:n_test]
X_valid = X_normalized[n_test:,:22]
Y_valid = Y[n_test:]

model = Sequential([])
model.add(Input(shape=(22)))
model.add(Dense(200, activation=custom_activation, kernel_initializer=tf.keras.initializers.RandomUniform(minval=-7.0, maxval=7.0, seed=gseed)))
model.add(Dense(200, activation=custom_activation))
model.add(Dense(1, activation="sigmoid"))

model.compile(optimizer='adam',
              loss='binary_crossentropy',
              metrics=["accuracy"])

model.summary()

def learning_rate(epoch):
    if epoch < 30 :
        return 1e-3
    elif epoch < 35 :
        return 1e-4
    else :
        return 1e-5

learning_rate_cb = callbacks.LearningRateScheduler(learning_rate)

history = model.fit(X_train, Y_train, epochs=40, batch_size=8, callbacks=[learning_rate_cb])

print(X_valid, Y_valid)
# Evaluate the model
loss, accuracy = model.evaluate(X_valid, Y_valid)
print("Test accuracy:", accuracy)
