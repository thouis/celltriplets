from keras.layers.core import Flatten, Dense
from keras.layers.convolutional import Convolution2D
from keras.layers.normalization import BatchNormalization
from keras.layers.pooling import MaxPooling2D
from keras.layers.advanced_activations import ELU
from keras.layers import Input, merge
from keras.models import Model
from keras.optimizers import SGD
import keras.backend as K
from keras.utils.visualize_util import plot
import numpy as np
import os

from generate_triplets import generate_triplets


def residual_block(input, num_feature_maps, filter_size=3):
    conv_1 = BatchNormalization(axis=1, mode=2)(input)
    conv_1 = ELU()(conv_1)
    conv_1 = Convolution2D(num_feature_maps, filter_size, filter_size,
                           border_mode='same', bias=False)(conv_1)

    conv_2 = BatchNormalization(axis=1, mode=2)(conv_1)
    conv_2 = ELU()(conv_2)
    conv_2 = Convolution2D(num_feature_maps, filter_size, filter_size,
                           border_mode='same', bias=True)(conv_2)

    return merge([input, conv_2], mode='sum')


def residual_tower(input, num_feature_maps, filter_size=3,
                   output_dim=64,
                   num_levels=4, per_level=3):
    for l in range(num_levels):
        for p in range(per_level):
            input = residual_block(input, num_feature_maps, filter_size=filter_size)
        # expand depth, and maxpool
        num_feature_maps *= 2
        input = Convolution2D(num_feature_maps, 1, 1, border_mode='same', bias=True)(input)
        input = MaxPooling2D()(input)
        input = ELU()(input)

    tmp = Flatten()(input)

    return Dense(output_dim=output_dim, activation=None)(tmp)


def triplet_L1_margin_loss(y_true, y_pred, margin=1.0):
    # y_pred is (batch_size, num_features)
    # it is stacked i1, i2, i3 with i1 == i2 != i3
    vec1 = y_pred[0::3, :]
    vec2 = y_pred[1::3, :]
    vec3 = y_pred[2::3, :]
    dists_1_2 = K.sum(K.abs(vec1 - vec2), axis=1)
    dists_1_3 = K.sum(K.abs(vec1 - vec3), axis=1)

    hinge_losses = K.maximum(0, margin + dists_1_2 - dists_1_3) + 0 * K.sum(y_true)
    m = K.mean(hinge_losses)
    return m


def batchgen(root, batchsize=32):
    gen = generate_triplets(root)
    while True:
        batch = []
        for b in range(batchsize):
            i1, i2, i3 = next(gen)
            batch += [i1, i2, i3]
        b = np.stack(batch)
        b = b.astype(np.float32) / (2 ** 12)
        yield b, np.zeros((batchsize * 3,))


if __name__ == "__main__":
    root = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')
    num_feature_maps = 64
    INPUT_SHAPE = (64, 64, 3)

    x = Input(shape=INPUT_SHAPE)
    x.tag.test_value = np.random.uniform(0, 1, (96, 64, 64, 1)).astype(np.float32)
    pre = Convolution2D(num_feature_maps, 5, 5, bias=True, border_mode='same')(x)
    post = residual_tower(pre, num_feature_maps, output_dim=64)
    model = Model(input=x, output=post)

    plot(model, to_file='model.png', show_shapes=True)

    model.compile(optimizer=SGD(lr=0.0005, clipnorm=1., momentum=0.9),
                  loss=triplet_L1_margin_loss)

    gen = batchgen(root)
    hist = model.fit_generator(gen, 3 * 4096, 1000, verbose=2)
