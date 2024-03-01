from tqdm import tqdm
import json
import pickle
import pandas as pd
import numpy as np

node_allele = pd.read_csv('../data_obtain/node_allele.txt', sep = ',')

freq_allels = set()
for alleles in tqdm(node_allele['node_value']):
    freq_allels.update(alleles.split(';'))

with open('../data_obtain/tree_dict.txt', 'rb') as file:
    parent_child = pickle.load(file)
with open('../data_obtain/tree_dict_op.txt', 'r') as file:
    c_p = file.read()
child_parent = json.loads(c_p)
del c_p

def find_local_root(node, allel):
    current_node = node
    if child_parent[node] != '':
        node_parent = child_parent[node]
        if allel in node_allele[node_allele['node_name'] == node_parent]['node_value'].values[0].split(';'):
            return find_local_root(node=node_parent, allel=allel)
        else:
            return current_node
    else:
        return current_node


def get_subtree(parent_node, allel, comp_subtree=[]):
    allel_nodes.remove(parent_node)
    comp_subtree.append(parent_node)
    if parent_node in parent_child.keys():
        for child in parent_child[parent_node]:
            if allel in node_allele[node_allele['node_name'] == child]['node_value'].values[0].split(';'):
                get_subtree(parent_node=child,
                            allel=allel, comp_subtree=comp_subtree)
    return comp_subtree

with open("../data_obtain/allel_st.txt", 'w') as file:
    file.write('allel' + ',' + 'nodes'+'\n')
    for allel in tqdm(freq_allels):
        allel_nodes = list(node_allele[node_allele['node_value'].str.contains(allel, regex=False)]['node_name'].values)
        if (len(allel_nodes) >= 32) and (len(allel_nodes) < len(node_allele['node_name']) * 0.9):
            while len(allel_nodes) != 0:
                allel_local_root = find_local_root(allel_nodes[0], allel)
                subtree = get_subtree(parent_node=allel_local_root, allel=allel, comp_subtree = [])
                if len(subtree) >= 32:
                    file.write(allel+','+';'.join(map(str, subtree))+'\n')
        allel_paths = {}