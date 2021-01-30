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

column_names = ['Pos', 'RefNuc', 'GenName', 'CodonNumber', 'RefCodon', 'RefAa', 'NucInCodon', 'AltNuc', 'AltCodon', 'AltAa', 'AaSub', 'NeighL',
                'NeighR', 'MutExp', 'MutObsNC']
df = pd.DataFrame(columns = column_names)
positions = []
nucleotids = []
altnuc = []
codons = []
for i in range(1, len(nucl_ref), 3):
    if i == 1:
        codons.append(None)
    if len(nucl_ref[i:i+3]) < 3:
        codons.append(None)
    else:
        codons.append(nucl_ref[i:i + 3])
nine_codons = []
codon_number = []
for i in range(len(codons)):
    if i == 0:
        nine_codons.append(codons[i])
        nine_codons.append(codons[i])
        nine_codons.append(codons[i])
        codon_number.append(None)
        codon_number.append(None)
        codon_number.append(None)
    else:
        if i == len(codons)-1:
            nine_codons.append(codons[i])
            nine_codons.append(codons[i])
            nine_codons.append(codons[i])
            codon_number.append(None)
            codon_number.append(None)
            codon_number.append(None)
        else:
            nine_codons.append(codons[i])
            nine_codons.append(codons[i])
            nine_codons.append(codons[i])
            nine_codons.append(codons[i])
            nine_codons.append(codons[i])
            nine_codons.append(codons[i])
            nine_codons.append(codons[i])
            nine_codons.append(codons[i])
            nine_codons.append(codons[i])
            codon_number.append(i)
            codon_number.append(i)
            codon_number.append(i)
            codon_number.append(i)
            codon_number.append(i)
            codon_number.append(i)
            codon_number.append(i)
            codon_number.append(i)
            codon_number.append(i)
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
df['GenName'][df.Pos<266] = '5UTR'
df['GenName'][(265 < df.Pos)&(df.Pos < 13468)] = 'orf1a'
df['GenName'][(13467 < df.Pos)&(df.Pos < 21563)] = 'orf1b'
df['GenName'][(21562 < df.Pos)&(df.Pos < 25393)] = 'S'
df['GenName'][(25392 < df.Pos)&(df.Pos < 26245)] = 'orf3a'
df['GenName'][(26244 < df.Pos)&(df.Pos < 26523)] = 'E'
df['GenName'][(26522 < df.Pos)&(df.Pos < 27202)] = 'M'
df['GenName'][(27201 < df.Pos)&(df.Pos < 27394)] = 'orf6'
df['GenName'][(27393 < df.Pos)&(df.Pos < 27756)] = 'orf7a'
df['GenName'][(27755 < df.Pos)&(df.Pos < 27894)] = 'orf7b'
df['GenName'][(27893 < df.Pos)&(df.Pos < 28274)] = 'orf8'
df['GenName'][(28273 < df.Pos)&(df.Pos < 29558)] = 'N'
df['GenName'][(29557 < df.Pos)&(df.Pos < 29675)] = 'orf10'
df['GenName'][(29674 < df.Pos)&(df.Pos <= 29903)] = '3UTR'

df['NucInCodon'][df.Pos % 3 == 0] = 2
df['NucInCodon'][(df.Pos+1) % 3 == 0] = 1
df['NucInCodon'][(df.Pos+2) % 3 == 0] = 3
df['RefCodon'] = nine_codons
df['CodonNumber'] = codon_number
for index, row in df.iterrows():
    if df['RefCodon'][index] != None:
        if df['NucInCodon'][index] == 1:
            if codons[int(df['CodonNumber'][index]) - 1] != None and df['Pos'][index] != 1:
                df['NeighR'][index] = codons[int(df['CodonNumber'][index])-1][1]
        elif df['NucInCodon'][index] == 2:
            if codons[int(df['CodonNumber'][index]) - 1] != None and df['Pos'][index] != 1:
                df['NeighR'][index] = codons[int(df['CodonNumber'][index]) - 1][2]
                df['NeighL'][index] = codons[int(df['CodonNumber'][index]) - 1][0]
        elif df['NucInCodon'][index] == 3:
            if codons[int(df['CodonNumber'][index]) - 1] != None and df['Pos'][index] != 1:
                df['NeighL'][index] = codons[int(df['CodonNumber'][index]) - 1][1]

    if df['RefCodon'][index] != None:
        coding_dna = Seq(df['RefCodon'][index])
        df['RefAa'][index] = coding_dna.translate()
        new_codon = df['RefCodon'][index][0:(df['NucInCodon'][index]-1)] + str(df['AltNuc'][index]) + df['RefCodon'][index][df['NucInCodon'][index]:4]
        df['AltCodon'][index] = new_codon
        AltAa = Seq(df['AltCodon'][index])
        df['AltAa'][index] = AltAa.translate()
    # if any(data['bp'] == df['Pos'][index]):
    #     NC_nucl_data = data.loc[(data['bp'] == df['Pos'][index])]
    #     NC_nucl_data = NC_nucl_data.loc[(data['snp'] == df['AltNuc'][index])]
    #     NC_nucl_data = NC_nucl_data.loc[(data['pos_in_codon'] == df['NucInCodon'][index])]
    #     df['MutObsNC'][index] = int(NC_nucl_data['snp_count'])
df['AaSub'][df.RefAa == df.AltAa] = 'S'
df['AaSub'][df.RefAa != df.AltAa] = 'NS'
df['AaSub'][df.RefAa == None] = None
df.to_csv(resultPath / 'ideal_table.csv')
print(df)