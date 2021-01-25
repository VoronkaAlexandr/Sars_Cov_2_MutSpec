import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
import numpy as np
from pathlib import Path
import seaborn as sns
from tqdm import tqdm

sns.set_style('whitegrid')

rootPath = Path('..')
dataPath = rootPath / 'data'
figuresPath = rootPath / 'figure'
resultPath = rootPath /'data_obtained'

data = pd.read_excel(dataPath / 'only_mutations.xls')
ref = (dataPath / 'Covid_ref.fasta')

def twelve_comp_mutspec(data):
    mutspec = {}
    for index, row in tqdm(data.iterrows()):
        subst = row['ref'] + row['snp']
        if subst in mutspec.keys():
            mutspec[subst] += row['snp_count']
        else:
            mutspec[subst] = row['snp_count']
    return mutspec

def twelve_comp_mutspec_NS_S(data):
    mutspec_S= twelve_comp_mutspec(data[data.NS_S ==  'S'])
    mutspec_NS = twelve_comp_mutspec(data[data.NS_S ==  'NS'])
    NS_mut_spec = pd.DataFrame([mutspec_NS])
    S_mut_spec = pd.DataFrame([mutspec_S])
    NS_mut_spec.to_csv(resultPath / 'NS_mut_spec.csv')
    S_mut_spec.to_csv(resultPath / 'S_mut_spec.csv')

    fig_S = plt.figure()
    plt.title('Synonymous MutSpec')
    plt.bar(mutspec_S.keys(), mutspec_S.values())
    fig_S.savefig(figuresPath / 'S_mut_spec.png')

    fig_NS = plt.figure()
    plt.title('NonSynonymous MutSpec')
    plt.bar(mutspec_NS.keys(), mutspec_NS.values())
    fig_NS.savefig(figuresPath / 'NS_mut_spec.png')


def _twelve_comp_mutspec_NS_S(data):
    mutspec_S = {}
    mutspec_NS = {}
    for index, row in data.iterrows():
        snp = row['ref'] + row['snp']
        if row['NS_S'] == 'S':
            if snp in mutspec_S.keys():
                mutspec_S[snp] += row['snp_count']
            else:
                mutspec_S[snp] = row['snp_count']
        if row['NS_S'] == 'NS':
            if snp in mutspec_NS.keys():
                mutspec_NS[snp] += row['snp_count']
            else:
                mutspec_NS[snp] = row['snp_count']
    NS_mut_spec = pd.DataFrame([mutspec_NS])
    S_mut_spec = pd.DataFrame([mutspec_S])
    NS_mut_spec.to_csv(resultPath / 'NS_mut_spec.csv')
    S_mut_spec.to_csv(resultPath / 'S_mut_spec.csv')

    fig_S = plt.figure()
    plt.title('Synonymous MutSpec')
    plt.bar(mutspec_S.keys(), mutspec_S.values())
    fig_S.savefig(figuresPath/ 'S_mut_spec.png')

    fig_NS = plt.figure()
    plt.title('NonSynonymous MutSpec')
    plt.bar(mutspec_NS.keys(), mutspec_NS.values())
    fig_NS.savefig(figuresPath/ 'NS_mut_spec.png')

def mutspec_per_gen(data, gen):
    mutspec = {}
    new_data = data.loc[data['ORF'] == gen]
    mutspec = twelve_comp_mutspec(new_data)
    # for index, row in new_data.iterrows():
    #     snp = row['ref'] + row['snp']
    #     if snp in mutspec.keys():
    #         mutspec[snp] += row['snp_count']
    #     else:
    #         mutspec[snp] = row['snp_count']
    gen_mut_spec = pd.DataFrame([mutspec])
    gen_mut_spec.to_csv(resultPath / 'gen_{}_mut_spec.csv'.format(gen))

    fig = plt.figure()
    plt.title('{} MutSpec'.format(gen))
    plt.bar(mutspec.keys(), mutspec.values())
    fig.savefig(figuresPath / 'gen_{}_mut_spec.png'.format(gen))

def mutspec_192(data, ref):
    mutspec = {}
    second_nucl_data = data.loc[data['pos_in_codon'] == 2]
    for index, row in second_nucl_data.iterrows():
        codon = ref[int(row['bp'])-2 : int(row['bp'])+1]
        snp_codon = ref[int(row['bp'])-2] + row['snp'] + ref[int(row['bp'])]
        k = codon + ':' + snp_codon
        if k in mutspec.keys():
            mutspec[k] += row['snp_count']
        else:
            mutspec[k] = row['snp_count']
    gen_mut_spec = pd.DataFrame([mutspec])
    gen_mut_spec.to_csv(resultPath / 'mut_spec_192.csv')

    fig, ax = plt.subplots(figsize=(10,30))
    ax.set_title('MutSpec192')
    # plt.bar(mutspec.keys(), mutspec.values(), width = 0.4)
    df= pd.DataFrame([(k,mutspec[k]) for k  in mutspec])
    ax.barh(df[0], df[1])
    plt.yticks((np.arange((len(df)))), df[0])
    ax.axes.set_ylim(-0.5,len(df)+0.5)

    plt.savefig(figuresPath / 'mutspec192.png')


#
# mutspec = twelve_comp_mutspec(data)
# surface_mut_spec = pd.DataFrame([mutspec])
# surface_mut_spec.to_csv(resultPath / 'surface_mut_spec.csv')
#
# fig = plt.figure()
# plt.title('Surface MutSpec')
# plt.bar(mutspec.keys(), mutspec.values())
# fig.savefig(figuresPath / 'surface_mut_spec.png')
#
# twelve_comp_mutspec_NS_S(data)
# for gen in data['ORF'].unique():
#      mutspec_per_gen(data, gen)
nucl_ref = ''
with open (ref, 'r') as ref:
    a = ref.readlines()
    for line in range(len(a)):
        if line != 0:
            nucl_ref += a[line].rstrip()
mutspec_192(data, nucl_ref)


