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
    fasta_path = 'data/raw'
    fastafiles = os.listdir(fasta_path)

    for fastafile in fastafiles:
        with open(fastafile, 'r') as f:
            records = SeqIO.parse(f, 'fasta')

            for record in records:
                kind, protein = determine(record.desc, HIV_dict)


def determine(desc, HIV_dict):
    # HIVの種別を判定
    if ["HIV-1", "Human immunodeficiency virus 1"] in desc:
        kind = "HIV-1"
    elif ["HIV-2", "Human immunodeficiency virus 2"] in desc:
        kind = "HIV-2"
    else:
        print("kind couldn't be determined."\
              "desc: {}".format(desc))

    for key, val in HIV_dict[kind].items():
        if val in desc:
            protein = key
        else:
            print("protein coundn't be determined."\
                  "desc: {}".format(desc))

    return kind, protein
