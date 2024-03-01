from ete3 import PhyloTree
import json
import pickle


def make_opposite_dict(dict):
    opposite_dict = {}
    for parent in dict.keys():
        if len(dict[parent])>1 and parent != "":
            for child in dict[parent]:
                opposite_dict[child] = parent
    return opposite_dict

def make_three_dict(path_to_tree, file_to_save, file_to_save_op):
    tree = PhyloTree(path_to_tree, format=1)
    tree_dict = {}

    for node in tree.traverse():
        if node.is_root():
            root = node.name
            #key = str(['',root])
            tree_dict[''] = set([root])
        if len(node.get_children()) >=1:
            for child in node.get_children():
                #key = str([node.name, child.name])
                #tree_dict[key] = tree.get_distance(root, child.name)
                if node.name not in tree_dict.keys():
                    tree_dict[node.name] = set([child.name])
                else:
                    tree_dict[node.name].add(child.name)
        #elif len(node.get_children()) == 0:
        #    #key = str([node.name, ""])
        #    #tree_dict[key] = tree.get_distance(root, node.name)
        #    tree_dict[node.name] = ''

    child_parent_dict = make_opposite_dict(tree_dict)
    child_parent_dict[root] = ''

    file = open(file_to_save, 'wb')
    pickle.dump(tree_dict, file)
    file.close()

    file = open(file_to_save_op, 'w')
    file.write(json.dumps(child_parent_dict))
    file.close()

if __name__ == "__main__":
    make_three_dict("../data/mulal.filtered.fasta.prank.anc.dnd", '../data_obtain/tree_dict.txt', '../data_obtain/tree_dict_op.txt')
