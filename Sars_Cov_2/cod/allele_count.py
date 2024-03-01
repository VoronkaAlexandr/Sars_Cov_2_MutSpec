import pandas as pd
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
polym_df_ns_pos = mut.loc[(mut['pos'] >= 11842) & (mut['pos'] <= 21551) & (mut['AaSub'] == "NS"),'pos'].unique().tolist()
del mut
nuc_in_codon_dict = dict(zip(aa_info['pos'], aa_info['NucInCodon']))
with open("../data/mulal.filtered.fasta.prank.anc.fas", 'r') as fasta_file, \
        open("../data_obtain/node_allele.txt", 'w') as write_file, open("../data_obtain/one_nod_per_allel.txt", 'w') as node_file:
    write_file.write('node_name' + ',' + 'node_value'+'\n')
    node_file.write('allel_name' + ',' + 'node_name' + '\n')
    for line in tqdm(fasta_file):
        line = line.rstrip()
        if line.startswith('>'):
            node_name = line[1:]
        elif not line.startswith('>'):
            node_value = []
            for pos in polym_df_ns_pos:
                if nuc_in_codon_dict[pos] == '1':
                    codon = line[pos:pos + 3]
                elif nuc_in_codon_dict[pos] == '2':
                    codon = line[pos - 1:pos + 2]
                elif nuc_in_codon_dict[pos] == '3':
                    codon = line[pos - 2:pos + 1]
                aa_num = aa_info.loc[aa_info['pos'] == pos, 'CodonNumber']
                aa_num = int(float(aa_num.values[0]))
                if ("_" not in codon) and ("-" not in codon):
                    node_file.write(str(str(aa_num) + codontab[codon]) + ',' + str(node_name) + '\n')
                    node_value.append(str(aa_num) + codontab[codon])
            write_file.write(node_name+','+';'.join(node_value)+'\n')