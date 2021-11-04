import pandas as pd
import tqdm
from collections import defaultdict

df = pd.read_csv('U_ideal_table.csv')
df = df.drop(columns=['Unnamed: 0', 'X.2', 'X.1', 'X'])

ff_aa = ['S', 'P', 'S_UC', 'A', 'T', 'L_CU', 'V', 'R_CG']
ff_loc = df[(df['RefAa'].isin(ff_aa))]
ff_codon = list(ff_loc.RefCodon.unique())
list_first_nuc_pos = list(df[df['NucInCodon'] == '1']['Pos'].unique())
first_nuc_pos = set(map(lambda x: x-266, list_first_nuc_pos))

def nuc_ff_codon_count(r):
    nuc_lst = []
    ff_lst = []
    codon_lst = []
    for line in tqdm.tqdm(r, total = 20000):
        nuc_dct = defaultdict(int)
        ff_dct = defaultdict(int)
        codon_dct = defaultdict(int)
        if line.startswith('>') == False:
            nuc_line = line.strip()
            nuc_line = nuc_line[265:29674]
            for i in range(len(nuc_line)):
                if nuc_line[i] != '-':
                    nuc_dct[nuc_line[i]]+=1
                    if (i in first_nuc_pos):
                        if(nuc_line[i+1] != '-') and (nuc_line[i+2] != '-'):
                            codon = nuc_line[i:i+3]
                            codon_dct[codon] += 1
                            if codon in ff_codon:
                                ff_dct[nuc_line[i+2]] += 1
            nuc_lst.append(nuc_dct)
            ff_lst.append(ff_dct)
            codon_lst.append(codon_dct)
    return nuc_lst, ff_lst, codon_lst

with open('splits/mulal_split_1.fasta', 'r') as r:
    nuc, ff, codons = nuc_ff_codon_count(r)
nuc_df = pd.DataFrame(nuc)
nuc_df.rename(columns={"T": "U"}, inplace = True)
ff_df = pd.DataFrame(ff)
ff_df.rename(columns={"T": "U"}, inplace = True)
codons_df = pd.DataFrame(codons)
codons_df.fillna(0, inplace = True)
for i in codons_df.columns:
    if 'U' in i:
        codons_df.rename(columns={i: i.replace('T', 'U')}, inplace = True)

nuc_df.to_csv('19.Nuc_Start.csv', index = False)
ff_df.to_csv('19.FF_Start.csv', index = False)
codons_df.to_csv('19.Codon_Start.csv', index = False)


with open('splits/mulal_split_2.fasta', 'r') as r:
    nuc, ff, codons = nuc_ff_codon_count(r)
nuc_df = pd.DataFrame(nuc)
nuc_df.rename(columns={"T": "U"}, inplace = True)
ff_df = pd.DataFrame(ff)
ff_df.rename(columns={"T": "U"}, inplace = True)
codons_df = pd.DataFrame(codons)
codons_df.fillna(0, inplace = True)
for i in codons_df.columns:
    if 'U' in i:
        codons_df.rename(columns={i: i.replace('T', 'U')}, inplace = True)

nuc_df.to_csv('19.Nuc_Inter.csv', index = False)
ff_df.to_csv('19.FF_Inter.csv', index = False)
codons_df.to_csv('19.Codon_Inter.csv', index = False)

with open('splits/mulal_split_3.fasta', 'r') as r:
    nuc, ff, codons = nuc_ff_codon_count(r)

nuc_df = pd.DataFrame(nuc)
nuc_df.rename(columns={"T": "U"}, inplace = True)
ff_df = pd.DataFrame(ff)
ff_df.rename(columns={"T": "U"}, inplace = True)
codons_df = pd.DataFrame(codons)
codons_df.fillna(0, inplace = True)
for i in codons_df.columns:
    if 'U' in i:
        codons_df.rename(columns={i: i.replace('T', 'U')}, inplace = True)

nuc_df.to_csv('19.Nuc_End.csv', index = False)
ff_df.to_csv('19.FF_End.csv', index = False)
codons_df.to_csv('19.Codon_End.csv', index = False)