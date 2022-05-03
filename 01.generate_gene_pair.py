# Usage: python

# import module
import numpy as np
import os, random, sys

out_dir = './'
out_dir = os.path.join(os.getcwd(), out_dir)

def get_gene_list(file_name):
    gene_list = {}
    files = open(file_name, 'r')
    for line in files:
        gene, symbol = line.strip().split('\t') # two columns: gene \t symbol
        gene_list[gene] = symbol
    files.close()
    return gene_list


def get_tf_target_list(file_name):
    train_gene_list = {}
    tf_target_list = {}
    files = open(file_name, 'r')
    for line in files:
        tf, target = line.strip().split('\t')
        train_gene_list[tf] = 1
        train_gene_list[target] = 1
        if tf in tf_target_list.keys():
            tf_target_list[tf] += ','+target
        else:
            tf_target_list[tf] = target
    files.close()
    return train_gene_list, tf_target_list


def generate_gene_pair(total_gene_list, train_gene_list, tf_target_list, out_file, sep_file):
    non_tf_list = [i for i in total_gene_list if i not in train_gene_list]
    files = open(out_file, 'w')
    seps = open(sep_file, 'w')
    # construct gene pair
    index = 0
    for tf, targetList in tf_target_list.items():
        seps.write('{}\n'.format(index))
        for target in targetList.split(','):
            files.write('{}\t{}\t{}\n'.format(tf, target, 1))
            files.write('{}\t{}\t{}\n'.format(tf, non_tf_list[random.randint(0, len(non_tf_list)-1)], 0))
            index += 1
    files.close()
    seps.write('{}\n'.format(index))
    seps.close()

def main():
    tf_target_file = 'Golden_standard_iGRN.txt'
    gene_file = 'ath_gene.txt'
    out_file = 'ath_gene_pairs.txt'
    out_file = out_dir+out_file
    sep_file = out_dir+'ath_sep_list.txt'

    total_gene_list = get_gene_list(gene_file)
    train_gene_list, tf_target_list = get_tf_target_list(tf_target_file)
    generate_gene_pair(total_gene_list, train_gene_list, tf_target_list, out_file, sep_file)

if __name__ == '__main__':
    main()
