import pandas as pd
from tqdm import tqdm
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform

mut = pd.read_csv('../data/u_mutation_dists.filtered.csv', index_col=0)
mut_syn = mut[mut['AaSub'].isin(['S', 'FF'])]
mut_syn['mut'] = mut_syn['parent_nucl'] + '>' + mut_syn['child_nucl']

subtrees = pd.read_csv('../data_obtain/allel_st.txt', sep=',')


allel_mutspec_full = pd.DataFrame(columns=('allel','A>C','A>G','A>U','C>A','C>G','C>U','G>A','G>C','G>U','U>A','U>C','U>G','sum_mutations'))
cap = pd.DataFrame({'mut': ['A>C', 'A>G', 'A>U',
                                    'C>A', 'C>G', 'C>U',
                                    'G>A', 'G>C', 'G>U',
                                    'U>A', 'U>C', 'U>G']})
allel_mutspec_bs = pd.DataFrame(columns=('allel','A>C','A>G','A>U','C>A','C>G','C>U','G>A','G>C','G>U','U>A','U>C','U>G'))
for index, row in tqdm(subtrees.iterrows()):
    allel = subtrees['allel'][index]
    subtree_nodes = subtrees['nodes'][index].split(';')
    allel_mut = mut_syn[mut_syn['child_node'].isin(subtree_nodes)][['mut']]
    allel_mutspec = pd.DataFrame(allel_mut['mut'].value_counts()).reset_index().rename(columns={'mut':'pos','index':'mut'})
    sum_mut = allel_mutspec['pos'].sum()
    if sum_mut>200:
        allel_mutspec['pos'] = allel_mutspec['pos']/sum_mut
        full_allel_mutspec = cap.merge(allel_mutspec, how='left', on='mut').fillna(0)
        t_allel_spec = full_allel_mutspec.set_index('mut').T
        t_allel_spec.columns.name = None
        t_allel_spec = t_allel_spec.reset_index(drop=True)
        t_allel_spec['allel'] = allel
        t_allel_spec['sum_mutations'] = sum_mut
        allel_mutspec_full = pd.concat([allel_mutspec_full, t_allel_spec])

        bs_inter_df = pd.DataFrame(columns=('A>C','A>G','A>U','C>A','C>G','C>U','G>A','G>C','G>U','U>A','U>C','U>G'))

        for n_bs in range(200):
            mut_sampled = allel_mut.sample(replace=True, frac=1)
            s_mutspec = pd.DataFrame(mut_sampled['mut'].value_counts()).reset_index().rename(
                    columns={'mut': 'pos', 'index': 'mut'})
            s_mutspec['pos'] = s_mutspec['pos'] / s_mutspec['pos'].sum()

            full_mutspec = cap.merge(s_mutspec, how='left', on='mut').fillna(0)

            t_allel_bs = full_mutspec.set_index('mut').T
            t_allel_bs.columns.name = None
            t_allel_bs = t_allel_bs.reset_index(drop=True)
            bs_inter_df = pd.concat([bs_inter_df, t_allel_bs])

        mean_bs_df = pd.DataFrame({'A>C': (bs_inter_df['A>C'].mean()),
                                       'A>G': (bs_inter_df['A>G'].mean()),
                                       'A>U': (bs_inter_df['A>U'].mean()),
                                       'C>A': (bs_inter_df['C>A'].mean()),
                                       'C>G': (bs_inter_df['C>G'].mean()),
                                       'C>U': (bs_inter_df['C>U'].mean()),
                                       'G>A': (bs_inter_df['G>A'].mean()),
                                       'G>C': (bs_inter_df['G>C'].mean()),
                                       'G>U': (bs_inter_df['G>U'].mean()),
                                       'U>A': (bs_inter_df['U>A'].mean()),
                                       'U>C': (bs_inter_df['U>C'].mean()),
                                       'U>G': (bs_inter_df['U>G'].mean())}, index=[0])
        mean_bs_df = mean_bs_df.reset_index(drop=True)
        mean_bs_df['allel'] = allel
        if (mean_bs_df==0).astype(int).sum(axis=1).values[0] <= 4:
            allel_mutspec_bs = pd.concat([allel_mutspec_bs, mean_bs_df])

allel_mutspec_full.to_csv('../data_obtain/allel_mutspec_12comp.csv')

allel_mutspec_bs = allel_mutspec_bs.set_index('allel')
table_cos_dist = pdist(allel_mutspec_bs, 'cosine')

df_cosdist = pd.DataFrame(squareform(table_cos_dist), index=allel_mutspec_bs.index, columns= allel_mutspec_bs.index)

df_cosdist.to_csv('../data_obtain/bs_mutspecs_allel.csv')


