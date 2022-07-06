import pandas as pd
import subprocess


file_path = 'data/accession_id_list.nbr'
columns = ['Representative', 'Neighbor', 'Host', 'Selected lineage',
           'Taxonomy name', 'Segment name']
df = pd.read_table(file_path, header=2, names=columns, index_col=False,
                   usecols=[1, 4])

df_HIV = df[df['Taxonomy name'].str.contains('HIV') |
            df['Taxonomy name'].str.contains('Human immunodeficiency virus')]

print(df_HIV['Neighbor'].values)
print(type(df_HIV['Neighbor'].values))

for id in df_HIV['Neighbor'].values:
    filename = f'data/raw/{id}.fasta'
    url = f'http://getentry.ddbj.nig.ac.jp/getentry/na/{id}'\
           '/?format=trans&filetype=text&trace=true&'\
           'show_suppressed=false&limit=10'
    subprocess.run(['wget', '-O', filename, url])
