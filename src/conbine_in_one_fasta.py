import os
import json
from Bio import SeqIO

from check_protein import HIV_dict

def main():
    fasta_dir = 'data/raw'
    fastafiles = os.listdir(fasta_dir)

    for fastafile in fastafiles:
        input_path = os.path.join(fasta_dir, fastafile)

        with open(input_path, 'r') as f_in:
            records = SeqIO.parse(f_in, 'fasta')

            for record in records:
                kind = determine_kind(record.description)
                whole_genome = is_whole_genome(record.description, kind)

                if whole_genome:
                    continue

                if kind == 'HIV-1':
                    out_path = 'data/one_file/hiv1'
                elif kind == 'HIV-2':
                    out_path = 'data/one_file/hiv2'
                else:
                    continue

                with open(out_path, 'a') as f_out:
                    SeqIO.write(record.seq, f_out, 'fasta')


def determine_kind(desc):
    """ HIV-1とHIV-2を判別する """
    kind = None
    keys_HIV1 = ["HIV-1", "Human immunodeficiency virus 1",
                 "Human immunodeficiency virus type 1"]
    keys_HIV2 = ["HIV-2", "Human immunodeficiency virus 2",
                 "Human immunodeficiency virus type 2"]

    for key in keys_HIV1:
        if key in desc:
            kind = 'HIV-1'

    if kind is None:
        for key in keys_HIV2:
            if key in desc:
                kind = "HIV-2"

    return kind

def is_whole_genome(desc, kind):
    json_path = 'references/HIV_dict.json'
    with open(json_path, 'r') as f:
        HIV_dict = json.load(f)

    whole_genome = True
    for key, vals in HIV_dict[kind].items():
        for val in vals:
            if val in desc:
                whole_genome = False

    return whole_genome


if __name__ == '__main__':
    main()