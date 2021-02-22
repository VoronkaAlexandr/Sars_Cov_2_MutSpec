import pandas as pd
from pathlib import Path
from Bio.Seq import Seq

rootPath = Path('..')
dataPath = rootPath / 'data'
figuresPath = rootPath / 'figure'
resultPath = rootPath /'data_obtained'

ref = (dataPath / 'Covid_ref.fasta')
data = pd.read_excel(dataPath / 'only_mutations.xls')

nucl_ref = ''
with open (ref, 'r') as ref:
    a = ref.readlines()
    for line in range(len(a)):
        if line != 0:
            nucl_ref += a[line].rstrip()

column_names = ['Pos', 'RefNuc', 'GenName', 'GenType', 'CodonNumber', 'RefCodon', 'RefAa', 'NucInCodon', 'AltNuc', 'AltCodon', 'AltAa', 'AaSub', 'NeighL',
                'NeighR', 'MutExp', 'MutObsNC']
df = pd.DataFrame(columns = column_names)
positions = []
nucleotids = []
altnuc = []
nuc_in_codon = []

for i in range(len(nucl_ref)):
    positions.append(i+1)
    positions.append(i+1)
    positions.append(i+1)
    nucleotids.append(nucl_ref[i])
    nucleotids.append(nucl_ref[i])
    nucleotids.append(nucl_ref[i])

df['Pos'] = positions
df['RefNuc'] = nucleotids

for i in range(0, len(nucleotids), 3):
    if nucleotids[i:i+3] == ['A', 'A', 'A']:
        altnuc.append('T')
        altnuc.append('G')
        altnuc.append('C')
    if nucleotids[i:i+3] == ['T', 'T', 'T']:
        altnuc.append('A')
        altnuc.append('G')
        altnuc.append('C')
    if nucleotids[i:i+3] == ['G', 'G', 'G']:
        altnuc.append('T')
        altnuc.append('A')
        altnuc.append('C')
    if nucleotids[i:i+3] == ['C', 'C', 'C']:
        altnuc.append('T')
        altnuc.append('G')
        altnuc.append('A')

df['AltNuc'] = altnuc

df['GenName'][df.Pos < 266] = '5UTR'
df['GenType'][df.Pos < 266] = 'untranslated'

for i in range(266 - 1):
    nuc_in_codon.append('None')
    nuc_in_codon.append('None')
    nuc_in_codon.append('None')

df['GenName'][(265 < df.Pos) & (df.Pos < 21554)] = 'ORF1ab'
df['GenType'][(265 < df.Pos) & (df.Pos < 21554)] = 'translated'
num_codon = 1
k = 1
for i in range(21554 - 265 - 1):
    if 1 <= k and k <= 3:
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        k += 1
    else:
        k = 1
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        k += 1

df['GenName'][(21553 < df.Pos) & (df.Pos < 21563)] = 'ORF1ab-S_space'
df['GenType'][(21553 < df.Pos) & (df.Pos < 21563)] = 'untranslated'

for i in range(21563 - 21553 - 1):
    nuc_in_codon.append('None')
    nuc_in_codon.append('None')
    nuc_in_codon.append('None')

df['GenName'][(21562 < df.Pos) & (df.Pos < 25385)] = 'S'
df['GenType'][(21562 < df.Pos) & (df.Pos < 25385)] = 'translated'

k = 1
for i in range(25385 - 21562 - 1):
    if 1 <= k and k <= 3:
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        k += 1
    else:
        k = 1
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        k += 1

df['GenName'][(25384 < df.Pos) & (df.Pos < 25393)] = 'S-ORF3a_space'
df['GenType'][(25384 < df.Pos) & (df.Pos < 25393)] = 'untranslated'

for i in range(25393 - 25384 - 1):
    nuc_in_codon.append('None')
    nuc_in_codon.append('None')
    nuc_in_codon.append('None')

df['GenName'][(25392 < df.Pos) & (df.Pos < 26221)] = 'ORF3a'
df['GenType'][(25392 < df.Pos) & (df.Pos < 26221)] = 'translated'

k = 1
for i in range(26221 - 25392 - 1):
    if 1 <= k and k <= 3:
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        k += 1
    else:
        k = 1
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        k += 1

df['GenName'][(26220 < df.Pos) & (df.Pos < 26245)] = 'ORF3a-E_space'
df['GenType'][(26220 < df.Pos) & (df.Pos < 26245)] = 'untranslated'

for i in range(26245 - 26220 - 1):
    nuc_in_codon.append('None')
    nuc_in_codon.append('None')
    nuc_in_codon.append('None')

df['GenName'][(26244 < df.Pos) & (df.Pos < 26473)] = 'E'
df['GenType'][(26244 < df.Pos) & (df.Pos < 26473)] = 'translated'

k = 1
for i in range(26473 - 26244 - 1):
    if 1 <= k and k <= 3:
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        k += 1
    else:
        k = 1
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        k += 1

df['GenName'][(26472 < df.Pos) & (df.Pos < 26523)] = 'E-M_space'
df['GenType'][(26472 < df.Pos) & (df.Pos < 26523)] = 'untranslated'

for i in range(26523 - 26472 - 1):
    nuc_in_codon.append('None')
    nuc_in_codon.append('None')
    nuc_in_codon.append('None')

df['GenName'][(26522 < df.Pos) & (df.Pos < 27192)] = 'M'
df['GenType'][(26522 < df.Pos) & (df.Pos < 27192)] = 'translated'

k = 1
for i in range(27192 - 26522 - 1):
    if 1 <= k and k <= 3:
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        k += 1
    else:
        k = 1
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        k += 1

df['GenName'][(27191 < df.Pos) & (df.Pos < 27202)] = 'M-ORF6_space'
df['GenType'][(27191 < df.Pos) & (df.Pos < 27202)] = 'untranslated'

for i in range(27202 - 27191 - 1):
    nuc_in_codon.append('None')
    nuc_in_codon.append('None')
    nuc_in_codon.append('None')

df['GenName'][(27201 < df.Pos) & (df.Pos < 27388)] = 'ORF6'
df['GenType'][(27201 < df.Pos) & (df.Pos < 27388)] = 'translated'

k = 1
for i in range(27388 - 27201 - 1):
    if 1 <= k and k <= 3:
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        k += 1
    else:
        k = 1
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        k += 1

df['GenName'][(27387 < df.Pos) & (df.Pos < 27394)] = 'ORF6-ORF7a_space'
df['GenType'][(27387 < df.Pos) & (df.Pos < 27394)] = 'untranslated'

for i in range(27394 - 27387 - 1):
    nuc_in_codon.append('None')
    nuc_in_codon.append('None')
    nuc_in_codon.append('None')

df['GenName'][(27393 < df.Pos) & (df.Pos < 27756)] = 'ORF7a'
df['GenType'][(27393 < df.Pos) & (df.Pos < 27756)] = 'translated'

k = 1
for i in range(27756 - 27393 - 1):
    if 1 <= k and k <= 3:
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        k += 1
    else:
        k = 1
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        k += 1

df['GenName'][(27755 < df.Pos) & (df.Pos < 27760)] = 'ORF7a-ORF7b_overlay'
df['GenType'][(27755 < df.Pos) & (df.Pos < 27760)] = 'translated_shift'

new_k = 1
for i in range(27760 - 27755 - 1):
    if (1 <= k and k <= 3) and (1 <= new_k and new_k <= 3):
        nuc_in_codon.append([k, new_k])
        nuc_in_codon.append([k, new_k])
        nuc_in_codon.append([k, new_k])
        k += 1
        new_k += 1
    elif k > 3:
        k = 1
        nuc_in_codon.append([k, new_k])
        nuc_in_codon.append([k, new_k])
        nuc_in_codon.append([k, new_k])
        k += 1
        new_k += 1
    elif new_k > 3:
        new_k = 1
        nuc_in_codon.append([k, new_k])
        nuc_in_codon.append([k, new_k])
        nuc_in_codon.append([k, new_k])
        k += 1
        new_k += 1

df['GenName'][(27759 < df.Pos) & (df.Pos < 27888)] = 'ORF7b'
df['GenType'][(27759 < df.Pos) & (df.Pos < 27888)] = 'translated'

for i in range(27888 - 27759 - 1):
    if 1 <= new_k and new_k <= 3:
        nuc_in_codon.append(new_k)
        nuc_in_codon.append(new_k)
        nuc_in_codon.append(new_k)
        new_k += 1
    else:
        new_k = 1
        nuc_in_codon.append(new_k)
        nuc_in_codon.append(new_k)
        nuc_in_codon.append(new_k)
        new_k += 1

df['GenName'][(27887 < df.Pos) & (df.Pos < 27894)] = 'ORF7b-ORF8_space'
df['GenType'][(27887 < df.Pos) & (df.Pos < 27894)] = 'untranslated'

for i in range(27894 - 27887 - 1):
    nuc_in_codon.append('None')
    nuc_in_codon.append('None')
    nuc_in_codon.append('None')

df['GenName'][(27893 < df.Pos) & (df.Pos < 28260)] = 'ORF8'
df['GenType'][(27893 < df.Pos) & (df.Pos < 28260)] = 'translated'

k = 1
for i in range(28260 - 27893 - 1):
    if 1 <= k and k <= 3:
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        k += 1
    else:
        k = 1
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        k += 1

df['GenName'][(28259 < df.Pos) & (df.Pos < 28274)] = 'ORF8-N_space'
df['GenType'][(27887 < df.Pos) & (df.Pos < 27894)] = 'untranslated'

for i in range(28274 - 28259 - 1):
    nuc_in_codon.append('None')
    nuc_in_codon.append('None')
    nuc_in_codon.append('None')

df['GenName'][(28273 < df.Pos) & (df.Pos < 29534)] = 'N'
df['GenType'][(28273 < df.Pos) & (df.Pos < 29534)] = 'translated'

k = 1
for i in range(29534 - 28273 - 1):
    if 1 <= k and k <= 3:
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        k += 1
    else:
        k = 1
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        k += 1

df['GenName'][(29533 < df.Pos) & (df.Pos < 29558)] = 'N-ORF10_space'
df['GenType'][(29533 < df.Pos) & (df.Pos < 29558)] = 'untranslated'

for i in range(29558 - 29533 - 1):
    nuc_in_codon.append('None')
    nuc_in_codon.append('None')
    nuc_in_codon.append('None')

df['GenName'][(29557 < df.Pos) & (df.Pos < 29675)] = 'ORF10'
df['GenType'][(29557 < df.Pos) & (df.Pos < 29675)] = 'translated'

k = 1
for i in range(29675 - 29557 - 1):
    if 1 <= k and k <= 3:
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        k += 1
    else:
        k = 1
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        nuc_in_codon.append(k)
        k += 1

df['GenName'][29674 < df.Pos] = '3UTR'
df['GenType'][29674 < df.Pos] = 'untranslated'

for i in range(29904 - 29674 - 1):
    nuc_in_codon.append('None')
    nuc_in_codon.append('None')
    nuc_in_codon.append('None')

df['NucInCodon'] = nuc_in_codon

neigh_l = []
neigh_r = []
neigh_l.append(None)
neigh_l.append(None)
neigh_l.append(None)

for_neighb = nucleotids[0:-3]
for i in for_neighb:
    neigh_l.append(i)

for_neighb = nucleotids[3:]
for i in for_neighb:
    neigh_r.append(i)
neigh_r.append(None)
neigh_r.append(None)
neigh_r.append(None)

df['NeighR'] = neigh_r
df['NeighL'] = neigh_l

codon_number = 1
for index, row in df.iterrows():
    if df['GenType'][index] == 'translated_shift':
        if df['NucInCodon'][index][0] == 3 and df['NucInCodon'][index][1] == 1:
            if df['Pos'][index] == 27756:
                df['RefCodon'][index] = [df['RefCodon'][index-6], df['RefNuc'][index]+df['RefNuc'][index+3]+df['RefNuc'][index+6]]
            elif df['Pos'][index] == 27759:
                df['RefCodon'][index] = [df['RefCodon'][index-6][0], df['RefNuc'][index]+df['RefNuc'][index+3]+df['RefNuc'][index+6]]
        if df['NucInCodon'][index][0] == 1 and df['NucInCodon'][index][1] == 2:
            df['RefCodon'][index] = [df['RefNuc'][index]+df['RefNuc'][index+3]+df['RefNuc'][index+6], df['RefCodon'][index-3][1]]
        if df['NucInCodon'][index][0] == 2 and df['NucInCodon'][index][1] == 3:
            df['RefCodon'][index] = [df['RefCodon'][index-3][0], df['RefCodon'][index-6][1]]
        coding_dna_a = Seq(str(df['RefCodon'][index][0]))
        coding_dna_b = Seq(str(df['RefCodon'][index][1]))
        df['RefAa'][index] = [coding_dna_a.translate(), coding_dna_b.translate()]
        new_codon_a = df['RefCodon'][int(index)][0][0:(int(df['NucInCodon'][int(index)][0])-1)] + str(df['AltNuc'][int(index)]) + df['RefCodon'][int(index)][0][int(df['NucInCodon'][int(index)][0]):4]
        new_codon_b = df['RefCodon'][int(index)][1][0:(int(df['NucInCodon'][int(index)][1])-1)] + str(df['AltNuc'][int(index)]) + df['RefCodon'][int(index)][1][int(df['NucInCodon'][int(index)][1]):4]
        df['AltCodon'][index] = [new_codon_a, new_codon_b]
        altaa_a = Seq(df['AltCodon'][index][0])
        altaa_b = Seq(df['AltCodon'][index][1])
        df['AltAa'][index] = [altaa_a.translate(), altaa_b.translate()]
    if index == 83277 or index == 83278 or index == 83279:
        df['RefCodon'][index] = df['RefCodon'][index-3][1]
    if index == 83280 or index == 83281 or index == 83282:
        df['RefCodon'][index] = df['RefCodon'][index-6][1]
    if (df['GenType'][index] == 'translated') and (df['Pos'][index] != 27760) and (df['Pos'][index] != 27761):
        if df['NucInCodon'][index] == 1:
            df['RefCodon'][index] = df['RefNuc'][index]+df['RefNuc'][index+3]+df['RefNuc'][index+6]
        elif df['NucInCodon'][index] == 2:
            df['RefCodon'][index] = df['RefCodon'][index-3]
        elif df['NucInCodon'][index] == 3:
            df['RefCodon'][index] = df['RefCodon'][index-6]
    if df['GenType'][index] == 'translated':
        coding_dna = Seq(str(df['RefCodon'][index]))
        df['RefAa'][index] = coding_dna.translate()
        new_codon = df['RefCodon'][int(index)][0:(int(df['NucInCodon'][int(index)])-1)] + str(df['AltNuc'][int(index)]) + df['RefCodon'][int(index)][int(df['NucInCodon'][int(index)]):4]
        df['AltCodon'][index] = new_codon
        AltAa = Seq(str(df['AltCodon'][index]))
        df['AltAa'][index] = AltAa.translate()
        if df['NucInCodon'][index] == 1 and df['NucInCodon'][index+5] == 2 and index+6 != 83265:
            df['CodonNumber'][index] = codon_number
            df['CodonNumber'][index+1] = codon_number
            df['CodonNumber'][index+2] = codon_number
            df['CodonNumber'][index+3] = codon_number
            df['CodonNumber'][index+4] = codon_number
            df['CodonNumber'][index+5] = codon_number
            df['CodonNumber'][index+6] = codon_number
            df['CodonNumber'][index+7] = codon_number
            df['CodonNumber'][index+8] = codon_number
            codon_number+=1
        elif index+6 == 83265:
            df['CodonNumber'][index] = codon_number
            df['CodonNumber'][index+1] = codon_number
            df['CodonNumber'][index+2] = codon_number
            df['CodonNumber'][index+3] = codon_number
            df['CodonNumber'][index+4] = codon_number
            df['CodonNumber'][index+5] = codon_number
    if df['GenType'][index] == 'translated_shift':
        if index == 83265:
            df['CodonNumber'][index] = [codon_number, codon_number+2]
            df['CodonNumber'][index+1] = [codon_number, codon_number+2]
            df['CodonNumber'][index+2] = [codon_number, codon_number+2]
            codon_number+=1
        elif index == 83274:
            df['CodonNumber'][index] = [codon_number, codon_number+2]
            df['CodonNumber'][index+1] = [codon_number, codon_number+2]
            df['CodonNumber'][index+2] = [codon_number, codon_number+2]
            df['CodonNumber'][index+3] = codon_number+2
            df['CodonNumber'][index+4] = codon_number+2
            df['CodonNumber'][index+5] = codon_number+2
            df['CodonNumber'][index+6] = codon_number+2
            df['CodonNumber'][index+7] = codon_number+2
            df['CodonNumber'][index+8] = codon_number+2
            codon_number = codon_number + 3
        if df['NucInCodon'][index][0] == 1 and df['NucInCodon'][index][1] == 2 and index == 83268:
            df['CodonNumber'][index] = [codon_number, codon_number+1]
            df['CodonNumber'][index+1] = [codon_number, codon_number+1]
            df['CodonNumber'][index+2] = [codon_number, codon_number+1]
            df['CodonNumber'][index+3] = [codon_number, codon_number+1]
            df['CodonNumber'][index+4] = [codon_number, codon_number+1]
            df['CodonNumber'][index+5] = [codon_number, codon_number+1]

for index, row in data.iterrows():
    a = df.loc[(df.Pos == data['bp'][index]) & (df['RefNuc'] == data['ref'][index]) & (df['AltNuc'] == data['snp'][index])]
    df['MutObsNC'][a.index.tolist()] = data['snp_count'][index]

df['AaSub'][df.RefAa == df.AltAa] = 'S'
df['AaSub'][df.RefAa != df.AltAa] = 'NS'
df['AaSub'][df.GenType == 'untranslated'] = None

df.to_csv(resultPath / 'ideal_table.csv')