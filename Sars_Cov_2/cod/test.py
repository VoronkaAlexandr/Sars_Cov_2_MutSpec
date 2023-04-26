import pandas as pd
import json
from tqdm.contrib import tzip
from numpy import dot
from numpy. linalg import norm

mut = pd.read_csv('../data/u_mutation_dists.filtered.csv', index_col=0)
mut_syn = mut[mut['AaSub'].isin(['S', 'FF'])]
mut_syn['mut'] = mut_syn['parent_nucl'] + '>' + mut_syn['child_nucl']

with open('../data_obtain/allel_st.txt', 'r') as file:
    s_t = file.read()

subtrees = json.loads(s_t)

allel_mutspec_full = pd.read_csv('../data_obtain/allel_mutspec_12comp.csv')

allel_mutspec_full = pd.DataFrame(columns=('allel','n_subtree','A>C','A>G','A>U','C>A','C>G','C>U','G>A','G>C','G>U','U>A','U>C','U>G','sum_mutations'))
cap = pd.DataFrame({'mut': ['A>C', 'A>G', 'A>U',
                                    'C>A', 'C>G', 'C>U',
                                    'G>A', 'G>C', 'G>U',
                                    'U>A', 'U>C', 'U>G']})

allel_list = list(allel_mutspec_full['allel'])
n_subtree_list = list(allel_mutspec_full['n_subtree'])

allel_allel_dict = {}
for allel_one, n_subtree_one in tzip(allel_list, n_subtree_list):
    if n_subtree_one == 0:
        subtree_nodes_one = subtrees[allel_one][n_subtree_one]
    else:
        subtree_nodes_one = subtrees[allel_one][n_subtree_one][0]

    allel_mut_one = mut_syn[mut_syn['child_node'].isin(subtree_nodes_one)][['pos', 'mut']]
    for allel_two, n_subtree_two in zip(allel_list, n_subtree_list):
        if n_subtree_two == 0:
            subtree_nodes_two = subtrees[allel_two][n_subtree_two]
        else:
            subtree_nodes_two = subtrees[allel_two][n_subtree_two][0]

        allel_mut_two = mut_syn[mut_syn['child_node'].isin(subtree_nodes_two)][['pos', 'mut']]

        dif = []
        for _ in range(1000):
            one_sampled = allel_mut_one.sample(replace=True, frac=1)
            two_sampled = allel_mut_two.sample(replace=True, frac=1)

            s_one_mutspec = one_sampled[['pos', 'mut']].groupby(by=["mut"]).count()
            s_two_mutspec = two_sampled[['pos', 'mut']].groupby(by=["mut"]).count()

            s_one_mutspec['pos'] = s_one_mutspec['pos'] / s_one_mutspec['pos'].sum()
            full_one_mutspec = cap.merge(s_one_mutspec, how='left', on='mut').fillna(0)
            s_two_mutspec['pos'] = s_two_mutspec['pos'] / s_two_mutspec['pos'].sum()
            full_two_mutspec = cap.merge(s_two_mutspec, how='left', on='mut').fillna(0)
            a = list(full_one_mutspec['pos'])
            b = list(full_two_mutspec['pos'])
            dif.append(dot(a, b) / (norm(a) * norm(b)))

        mean = sum(dif) / len(dif)
        variance = sum([((x - mean) ** 2) for x in dif]) / len(dif)
        std = variance ** 0.5

        allel_allel_dict[allel_one+'_'+str(n_subtree_one)+'__'+allel_two+''+str(n_subtree_two)] = [mean, std]

file = open("../data_obtain/allel_allel_dist.txt", 'w')
file.write(json.dumps(allel_allel_dict))
file.close()