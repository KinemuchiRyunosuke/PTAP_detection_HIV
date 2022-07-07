import os
from Bio import SeqIO


HIV_dict = {\
    "HIV-1": {"gag": ["gag", "Gag"],
              "pol": ["pol", "Pol"],
              "vif": ["vif", "Vif", "viral infectivity factor"],
              "vpr": ["vpr", "Vpr"],
              "tat": ["tat", "Tat"],
              "rev": ["rev", "Rev"],
              "vpu": ["vpu", "Vpu"],
              "envelope glycoprotein":
                ["env", "Env"],
              "nef": ["nef", "Nef"]},
    "HIV-2": {"gag": ["gag", "Gag"],
              "vif": ["vif", "Vif"],
              "vpx": ["vpx", "Vpx"],
              "vpr": ["vpr", "Vpr"],
              "tat": ["tat", "Tat"],
              "rev": ["rev", "Rev"],
              "envelope polyprotein":
                ["envelope polyprotein", "env protein"],
              "nef": ["nef", "Nef"]}
}


def main():
    fasta_dir = 'data/raw'
    fastafiles = os.listdir(fasta_dir)

    for fastafile in fastafiles:
        fasta_path = os.path.join(fasta_dir, fastafile)
        with open(fasta_path, 'r') as f:
            records = SeqIO.parse(f, 'fasta')

            for record in records:
                kind, protein = determine(record.description, HIV_dict)


def determine(desc, HIV_dict):
    # HIVの種別を判定
    kind, protein = None, None
    keys_HIV1 = ["HIV-1", "Human immunodeficiency virus 1"]
    keys_HIV2 = ["HIV-2", "Human immunodeficiency virus 2"]

    for key in keys_HIV1:
        if key in desc:
            kind = "HIV-1"

    if kind is None:
        for key in keys_HIV2:
            if key in desc:
                kind = "HIV-2"

    if kind is None:
        print("kind couldn't be determined. "\
              "desc: {}".format(desc))
        return kind, protein

    # タンパク質を判定
    for key, vals in HIV_dict[kind].items():
        for val in vals:
            if val in desc:
                protein = key

    if protein is None:
        print("protein coundn't be determined. "\
              "desc: {}".format(desc))

    return kind, protein

if __name__ == '__main__':
    main()
