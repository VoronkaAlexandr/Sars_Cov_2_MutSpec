import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
import numpy as np

data = pd.read_excel('/home/bufnita/PycharmProjects/Sars_Cov_2/data/only_mutations.xls')
ref = ('/home/bufnita/PycharmProjects/Sars_Cov_2/data/Covid_ref.fasta')

def twelve_comp_mutspec(data):
    mutspec = {}
    for index, row in data.iterrows():
        if row['ref'] + row['snp'] in mutspec.keys():
            mutspec[row['ref'] + row['snp']] += row['snp_count']
        else:
            mutspec[row['ref'] + row['snp']] = row['snp_count']
    surface_mut_spec = pd.DataFrame([mutspec])
    surface_mut_spec.to_csv('/home/bufnita/PycharmProjects/Sars_Cov_2/data_obtained/surface_mut_spec.csv')

    fig = plt.figure()
    plt.title('Surface MutSpec')
    plt.bar(mutspec.keys(), mutspec.values())
    fig.savefig('/home/bufnita/PycharmProjects/Sars_Cov_2/figure/surface_mut_spec.png')

def twelve_comp_mutspec_NS_S(data):
    mutspec_S = {}
    mutspec_NS = {}
    for index, row in data.iterrows():
        if row['NS_S'] == 'S':
            if row['ref'] + row['snp'] in mutspec_S.keys():
                mutspec_S[row['ref'] + row['snp']] += row['snp_count']
            else:
                mutspec_S[row['ref'] + row['snp']] = row['snp_count']
        if row['NS_S'] == 'NS':
            if row['ref'] + row['snp'] in mutspec_NS.keys():
                mutspec_NS[row['ref'] + row['snp']] += row['snp_count']
            else:
                mutspec_NS[row['ref'] + row['snp']] = row['snp_count']
    NS_mut_spec = pd.DataFrame([mutspec_NS])
    S_mut_spec = pd.DataFrame([mutspec_S])
    NS_mut_spec.to_csv('/home/bufnita/PycharmProjects/Sars_Cov_2/data_obtained/NS_mut_spec.csv')
    S_mut_spec.to_csv('/home/bufnita/PycharmProjects/Sars_Cov_2/data_obtained/S_mut_spec.csv')

    fig_S = plt.figure()
    plt.title('Synonymous MutSpec')
    plt.bar(mutspec_S.keys(), mutspec_S.values())
    fig_S.savefig('/home/bufnita/PycharmProjects/Sars_Cov_2/figure/S_mut_spec.png')

    fig_NS = plt.figure()
    plt.title('NonSynonymous MutSpec')
    plt.bar(mutspec_NS.keys(), mutspec_NS.values())
    fig_NS.savefig('/home/bufnita/PycharmProjects/Sars_Cov_2/figure/NS_mut_spec.png')

def mutspec_per_gen(data, gen):
    mutspec = {}
    new_data = data.loc[data['ORF'] == gen]
    for index, row in new_data.iterrows():
        if row['ref'] + row['snp'] in mutspec.keys():
            mutspec[row['ref'] + row['snp']] += row['snp_count']
        else:
            mutspec[row['ref'] + row['snp']] = row['snp_count']
    gen_mut_spec = pd.DataFrame([mutspec])
    gen_mut_spec.to_csv('/home/bufnita/PycharmProjects/Sars_Cov_2/data_obtained/gen_{}_mut_spec.csv'.format(gen))

    fig = plt.figure()
    plt.title('{} MutSpec'.format(gen))
    plt.bar(mutspec.keys(), mutspec.values())
    fig.savefig('/home/bufnita/PycharmProjects/Sars_Cov_2/figure/gen_{}_mut_spec.png'.format(gen))

def mutspec_192(data, ref):
    mutspec = {}
    second_nucl_data = data.loc[data['pos_in_codon'] == 2]
    for index, row in second_nucl_data.iterrows():
        if ref[int(row['bp'])-2 : int(row['bp'])+1] + ':' + ref[int(row['bp'])-2] + row['snp'] + ref[int(row['bp'])] in mutspec.keys():
            mutspec[ref[int(row['bp'])-2 : int(row['bp'])+1] + ':' + ref[int(row['bp'])-2] + row['snp'] + ref[int(row['bp'])]] += row['snp_count']
        else:
            mutspec[ref[int(row['bp'])-2 : int(row['bp'])+1] + ':' + ref[int(row['bp'])-2] + row['snp'] + ref[int(row['bp'])]] = row['snp_count']
    gen_mut_spec = pd.DataFrame([mutspec])
    gen_mut_spec.to_csv('/home/bufnita/PycharmProjects/Sars_Cov_2/data_obtained/mut_spec_192.csv')

    fig = plt.figure()
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_title('MutSpec192')
    ax.bar(mutspec.keys(), mutspec.values(), width = 0.4)
    fig.set_figwidth(12)
    fig.set_figheight(6)
    fig.savefig('/home/bufnita/PycharmProjects/Sars_Cov_2/figure/mutspec192.png')



twelve_comp_mutspec(data)
twelve_comp_mutspec_NS_S(data)
for gen in data['ORF'].unique():
     mutspec_per_gen(data, gen)
nucl_ref = ''
with open (ref, 'r') as ref:
    a = ref.readlines()
    for line in range(len(a)):
        if line != 0:
            nucl_ref += a[line].rstrip()
mutspec_192(data, nucl_ref)

