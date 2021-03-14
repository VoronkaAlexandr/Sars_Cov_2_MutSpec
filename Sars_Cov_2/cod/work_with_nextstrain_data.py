import pandas as pd
import yaml
import ast

nucleotide_mutations = []
positions = []
df = pd.read_csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data/nextstrain_data.csv')
for line in df['Mutations'].fillna('No Mutations'):
    mut = []
    pos = []
    if line != 'No Mutations' and 'nuc' in line:
        line = ast.literal_eval(line)
        for single_mut in line['nuc']:
            mut.append((single_mut[0], single_mut[-1]))
            pos.append(single_mut[1:-1])
        nucleotide_mutations.append(mut)
        positions.append(pos)
    else:
        nucleotide_mutations.append(None)
        positions.append(None)
df['NucMut'] = nucleotide_mutations
df['Pos'] = positions
df = df.drop(df[df.NucMut.isna()].index)
df = df.reset_index(drop=True)
ref_nuc = []
alt_nuc = []
single_pos = []
continent = []
country = []
time = []
nuc_len = []
for line in df['NucMut']:
    nuc_len.append(len(line))
    for i in line:
        ref_nuc.append(i[0])
        alt_nuc.append(i[1])
for line in df['Pos']:
    for i in line:
        single_pos.append(i)
for index_len in range(len(nuc_len)):
    for _ in range(nuc_len[index_len]):
        country.append(df['Place'][index_len][0])
        continent.append(df['Place'][index_len][1])
        time.append(df['Time'][index_len])
new_df = pd.DataFrame({'Pos' : single_pos, 'RefNuc' : ref_nuc, 'AltNuc' : alt_nuc, 'Time' : time, 'Country' : country, 'Continent' : continent})
new_df.to_csv('D:/Sars_Cov_2_MutSpec-main/Sars_Cov_2_MutSpec-main/Sars_Cov_2/data/modernized_nextstrain.csv')