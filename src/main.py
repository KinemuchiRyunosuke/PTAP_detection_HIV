import os
import pandas as pd
import pickle
import tensorflow as tf

from Bio import SeqIO
from models.transformer import BinaryClassificationTransformer

from features.preprocessing import Vocab, add_class_token


num_words = 25
batch_size = 1024
hopping_num = 2
hidden_dim = 904
lr = 2.03e-5
head_num = 8
dropout_rate = 0.04
threshold = 0.5

protein_dict = {\
    "HIV-1": {"gag": ["gag", "Gag", "GAG"],
              "pol": ["pol", "Pol", "POL"],
              "vif": ["vif", "Vif", "VIF", "viral infectivity factor"],
              "vpr": ["vpr", "Vpr", "VPR", "vpR"],
              "tat": ["tat", "Tat", "TAT"],
              "rev": ["rev", "Rev", "REV"],
              "vpu": ["vpu", "Vpu", "VPU"],
              "envelope glycoprotein":
                ["env", "Env", "ENV", "gp120", "gp41", "gp160"],
              "nef": ["nef", "Nef", "NEF"]},
    "HIV-2": {"gag": ["gag", "Gag", "GAG"],
              "pol": ["pol", "Pol", "POL"],
              "vif": ["vif", "Vif", "VIF"],
              "vpx": ["vpx", "Vpx", "VPX"],
              "vpr": ["vpr", "Vpr", "VPR"],
              "tat": ["tat", "Tat", "TAT"],
              "rev": ["rev", "Rev", "REV"],
              "envelope polyprotein":
                ["env", "Env", "ENV"],
              "nef": ["nef", "Nef", "NEF"]}
}


def main():
    vocab_path = 'references/vocab.pickle'
    with open(vocab_path, 'rb') as f:
        tokenizer = pickle.load(f)
    vocab = Vocab(tokenizer)

    fasta_dir = 'data/raw'
    fastafiles = os.listdir(fasta_dir)

    model = create_model()
    checkpoint_path = 'models/'
    model.load_weights(checkpoint_path)

    y_pos = []
    for fastafile in fastafiles:
        fasta_path = os.path.join(fasta_dir, fastafile)

        for (fragment, kind, protein) in fragment_generator(fasta_path, 26):
            fragment = vocab.encode(fragment)
            y_pred = model.predict(fragment)
            y_pred = y_pred[0]
            if y_pred >= threshold:
                y_pos.append([fragment, kind, protein])

    df = pd.DataFrame(y_pos, columns=['fragment', 'kind', 'protein'])
    df.to_csv('reports/results/false_positive.csv')


def fragment_generator(fastafile, length):
    with open(fastafile, 'r') as f:
        records = SeqIO.parse(f, 'fasta')

        for record in records:
            for i in range(len(record.seq) - length + 1):
                fragment = record.seq[i:(i+length)]

                kind, protein = determine(record.description, protein_dict)
                yield (fragment, kind, protein)

def determine(desc, protein_dict):
    # HIV-1とHIV-2を判別
    keywords_HIV1 = ["HIV-1", "Human immunodeficiency virus 1",
                     "Human immunodeficiency virus type 1"]
    keywords_HIV2 = ["HIV-2", "Human immunodeficiency virus 2",
                     "Human immunodeficiency virus type 2"]

    def have_keyword(desc, keywords):
        have_keyword = False

        for keyword in keywords:
            if keyword in desc:
                have_keyword = True
                break

        return have_keyword

    if have_keyword(desc, keywords_HIV1):
        kind = "HIV-1"
    elif have_keyword(desc, keywords_HIV2):
        kind = "HIV-2"

    # タンパク質名を判別
    protein_name = None
    for key, proteins in protein_dict[kind].items():
        for protein in proteins:
            if protein in desc:
                protein_name = key
                break

        if protein_name is not None:
            break

    return kind, protein_name

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


if __name__ == '__main__':
    main()