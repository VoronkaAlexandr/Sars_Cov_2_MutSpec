import pandas as pd
import json
from tqdm import tqdm

codontab = {
    'TCA': 'Ser',    # Serina
    'TCC': 'Ser',    # Serina
    'TCG': 'Ser',    # Serina
    'TCT': 'Ser',    # Serina
    'TTC': 'Phe',    # Fenilalanina
    'TTT': 'Phe',    # Fenilalanina
    'TTA': 'Leu',    # Leucina
    'TTG': 'Leu',    # Leucina
    'TAC': 'Tyr',    # Tirosina
    'TAT': 'Tyr',    # Tirosina
    'TAA': '*',    # Stop
    'TAG': '*',    # Stop
    'TGC': 'Cys',    # Cisteina
    'TGT': 'Cys',    # Cisteina
    'TGA': '*',    # Stop
    'TGG': 'Trp',    # Triptofano
    'CTA': 'Leu',    # Leucina
    'CTC': 'Leu',    # Leucina
    'CTG': 'Leu',    # Leucina
    'CTT': 'Leu',    # Leucina
    'CCA': 'Pro',    # Prolina
    'CCC': 'Pro',    # Prolina
    'CCG': 'Pro',    # Prolina
    'CCT': 'Pro',    # Prolina
    'CAC': 'His',    # Histidina
    'CAT': 'His',    # Histidina
    'CAA': 'Gln',    # Glutamina
    'CAG': 'Gln',    # Glutamina
    'CGA': 'Arg',    # Arginina
    'CGC': 'Arg',    # Arginina
    'CGG': 'Arg',    # Arginina
    'CGT': 'Arg',    # Arginina
    'ATA': 'Ile',    # Isoleucina
    'ATC': 'Ile',    # Isoleucina
    'ATT': 'Ile',    # Isoleucina
    'ATG': 'Met',    # Methionina
    'ACA': 'Thr',    # Treonina
    'ACC': 'Thr',    # Treonina
    'ACG': 'Thr',    # Treonina
    'ACT': 'Thr',    # Treonina
    'AAC': 'Asn',    # Asparagina
    'AAT': 'Asn',    # Asparagina
    'AAA': 'Lys',    # Lisina
    'AAG': 'Lys',    # Lisina
    'AGC': 'Ser',    # Serina
    'AGT': 'Ser',    # Serina
    'AGA': 'Arg',    # Arginina
    'AGG': 'Arg',    # Arginina
    'GTA': 'Val',    # Valina
    'GTC': 'Val',    # Valina
    'GTG': 'Val',    # Valina
    'GTT': 'Val',    # Valina
    'GCA': 'Ala',    # Alanina
    'GCC': 'Ala',    # Alanina
    'GCG': 'Ala',    # Alanina
    'GCT': 'Ala',    # Alanina
    'GAC': 'Asp',    # Acido Aspartico
    'GAT': 'Asp',    # Acido Aspartico
    'GAA': 'Glu',    # Acido Glutamico
    'GAG': 'Glu',    # Acido Glutamico
    'GGA': 'Gly',    # Glicina
    'GGC': 'Gly',    # Glicina
    'GGG': 'Gly',    # Glicina
    'GGT': 'Gly'     # Glicina
}

mut = pd.read_csv("../data/u_mutation_dists.filtered.csv")
mut['pos'] = mut['pos'] -1
aa_info = mut[['pos', 'NucInCodon', 'CodonNumber']].drop_duplicates()
polym_df_ns_pos = mut.loc[(mut['pos'] >= 13442) & (mut['pos'] <= 16236) & (mut['AaSub'] == "NS"),'pos'].unique().tolist()
del mut
node_allele = {}
with open("../data/mulal.filtered.fasta.prank.anc.fas") as fasta_file:
    for line in tqdm(fasta_file):
        line = line.rstrip()
        if not line.startswith('>'):
            for pos in polym_df_ns_pos:
                if aa_info.loc[aa_info['pos'] == pos, 'NucInCodon'].values[0] == '1':
                    codon = line[pos:pos + 3]
                elif aa_info.loc[aa_info['pos'] == pos, 'NucInCodon'].values[0] == '2':
                    codon = line[pos - 1:pos + 2]
                elif aa_info.loc[aa_info['pos'] == pos, 'NucInCodon'].values[0] == '3':
                    codon = line[pos - 2:pos + 1]
                aa_num = aa_info.loc[aa_info['pos'] == pos, 'CodonNumber']
                if ("_" not in codon) and ("-" not in codon):
                    if node_name in node_allele.keys():
                        node_allele[node_name].append(str(aa_num) + codontab[codon])

                    # if line[pos] != 'T':
                    #    node_allele[node_name].append(str(pos) + line[pos])
                    # else:
                    #    node_allele[node_name].append(str(pos) + 'U')
                    elif node_name not in node_allele.keys():
                        node_allele[node_name] = [str(pos) + codontab[codon]]
                    # if line[pos] != 'T':
                    #    node_allele[node_name] = [str(pos) + line[pos]]
                    # else:
                    #    node_allele[node_name] = [str(pos) + 'U']
        else:
            node_name = line[1:]
del polym_df_ns_pos, aa_info, codontab

#file = open("../data_obtain/node_allele.txt", 'w')
#file.write(json.dumps(node_allele))
#file.close()

with open("../data_obtain/node_allele.txt", 'w') as write_file:
    write_file.write(json.dumps(node_allele))
