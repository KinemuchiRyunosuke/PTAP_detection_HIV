import os
from Bio import SeqIO


HIV_dict = {\
    "HIV-1": {"gag": ["gag"],
              "pol": ["pol"],
              "vif": ["vif"],
              "vpr": ["vpr"],
              "tat": ["tat"],
              "rev": ["rev"],
              "vpu": ["vpu"],
              "envelope glycoprotein": ["envelope glycoprotein"],
              "nef": ["nef"]},
    "HIV-2": {"gag": ["gag"],
              "vif": ["vif"],
              "vpx": ["vpx"],
              "vpr": ["vpr"],
              "tat": ["tat"],
              "rev": ["rev"],
              "envelope polyprotein": "envelope polyprotein",
              "nef": ["nef"]}
}


def main():
    fasta_dir = 'data/raw'
    fastafiles = os.listdir(fasta_dir)

    for fastafile in fastafiles:
        fasta_path = os.path.join(fasta_dir, fastafile)
        with open(fasta_path, 'r') as f:
            records = SeqIO.parse(f, 'fasta')

            for record in records:
                kind, protein = determine(record.desc, HIV_dict)


def determine(desc, HIV_dict):
    # HIVの種別を判定
    kind = None
    keys_HIV1 = ["HIV-1", "Human immunodeficiency virus 1"]
    keys_HIV2 = ["HIV-2", "Human immunodeficiency virus 2"]

    for key in keys_HIV1:
        if key in desc:
            kind = "HIV-1"

    if kind is not None:
        for key in keys_HIV2:
            if key in desc:
                kind = "HIV-2"
    else:
        print("kind couldn't be determined."\
              "desc: {}".format(desc))
        return 0

    # タンパク質を判定
    protein = None
    for key, vals in HIV_dict[kind].items():
        for val in vals:
            if val in desc:
                protein = key

    if protein is None:
        print("protein coundn't be determined."\
              "desc: {}".format(desc))

    return kind, protein
