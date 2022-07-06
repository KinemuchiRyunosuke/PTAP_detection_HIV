import os
import numpy as np
import tensorflow as tf

from Bio import SeqIO
from models.transformer import BinaryClassificationTransformer


num_words = 25
batch_size = 1024
hopping_num = 2
hidden_dim = 904
lr = 2.03e-5
head_num = 8
dropout_rate = 0.04
threshold = 0.5


def main():
    fasta_path = 'data/raw'
    fastafiles = os.listdir(fasta_path)

    model = create_model()
    checkpoint_path = 'models/saved_model.pb'
    model.load_weights(checkpoint_path)

    for fastafile in fastafiles:
        for (fragment, disc) in fragment_generator(fastafile, 26):
            y_pred = model.predict(fragment)
            if y_pred > threshold:
                y_pos += fragment


        y_pred = model.predict(fragment_generator(fastafile, 26))
        y_pred = np.squeeze(y_pred)
        y_pred = (y_pred > threshold).astype(int)
        y_pos = y_pred[y_pred == 1]



def fragment_generator(fastafile, length):
    with open(fastafile, 'r') as f:
        records = SeqIO.parse(f, 'fasta')

        for record in records:
            for i in range(len(record.seq) - length + 1):
                fragment = record.seq[i:(i+length)]
                yield (fragment, record.disc)

def create_model():
    """ モデルを定義する """
    model = BinaryClassificationTransformer(
                vocab_size=num_words,
                hopping_num=hopping_num,
                head_num=head_num,
                hidden_dim=hidden_dim,
                dropout_rate=dropout_rate)
    model.compile(optimizer=tf.keras.optimizers.Adam(
                                learning_rate=lr),
                 loss='binary_crossentropy')

    return model