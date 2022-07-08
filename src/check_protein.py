import os
from Bio import SeqIO


HIV_dict = {\
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
                ["envelope polyprotein", "env protein",
                 "envelope glycoprotein"],
              "nef": ["nef", "Nef", "NEF"]}
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
    keys_HIV1 = ["HIV-1", "Human immunodeficiency virus 1",
                 "Human immunodeficiency virus type 1"]
    keys_HIV2 = ["HIV-2", "Human immunodeficiency virus 2",
                 "Human immunodeficiency virus type 2"]

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
