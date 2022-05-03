
import numpy as np
import os, sys, math, glob
import matplotlib.pyplot as plt
import multiprocessing as mp


def get_tf_target_list(file_name):
    tf_target_list = []
    files = open(file_name, 'r')
    for line in files:
        tf, target, label = line.strip().split('\t')
        tf_target_list.append([tf, target, label])
    files.close()
    return np.array(tf_target_list)


def get_gene_expression_data(file_name):
    expr_list = []
    index = 0
    gene_index = {}
    files = open(file_name, 'r')
    for line in files:
        lines = line.strip().split('\t')
        expr_list.append(np.array(lines[2:], dtype='float32'))
        gene_index[lines[0]] = index
        index += 1
    files.close()
    return np.array(expr_list), gene_index


def get_sep_index(file_name):
    sep_list = []
    files = open(file_name, 'r')
    for line in files:
        sep_list.append(int(line))
    return np.array(sep_list)


def gene_train_data(sc_data, sc_index, bulk_data, bulk_index, sep_index, tf_target_list, save_dir):
    for i in range(len(sep_index)-1):
        #print(i)
        start_index = sep_index[i]
        end_index = sep_index[i+1]

        data_x = []
        data_y = []
        for tf, target, label in tf_target_list[start_index: end_index]:
            sc_flag = 0
            bulk_flag = 0
            #print('{}\t{}\n'.format(tf, target))
            if tf in sc_index.keys() and target in sc_index.keys():
                tf_sc = np.log10( sc_data[sc_index[tf]].reshape(-1)[2:] + 10 ** -2)
                target_sc = np.log10(sc_data[sc_index[target]].reshape(-1)[2:] + 10 ** -2)
                hist_sc = np.histogram2d(tf_sc, target_sc, bins=32)
                hist_sc_t = hist_sc[0].T
                hist_sc_norm = (np.log10(hist_sc_t / sc_data.shape[1] + 10 ** -4) + 4) / 4
                sc_flag = 1
                #print(tf_sc)

            if tf in bulk_index.keys() and target in bulk_index.keys():
                tf_bulk = np.log10( bulk_data[bulk_index[tf]].reshape(-1)[2:] + 10 ** -2)
                target_bulk = np.log10(bulk_data[bulk_index[target]].reshape(-1)[2:] + 10 ** -2)
                hist_bulk = np.histogram2d(tf_bulk, target_bulk, bins=32)
                hist_bulk_t = hist_bulk[0].T
                hist_bulk_norm = (np.log10(hist_bulk_t / bulk_data.shape[1] + 10 ** -4) + 4) / 4
                bulk_flag = 1

            #if sc_flag == 0 and bulk_flag == 1:
            #    x_data = hist_bulk_norm
            #elif sc_flag == 1 and bulk_flag == 0:
            #    x_data = hist_sc_norm
            #elif sc_flag == 1 and bulk_flag == 1:
            #    x_data = (hist_sc_norm + hist_bulk_norm)/2.0
            if sc_flag == 1 and bulk_flag == 1:
                x_data = np.concatenate((hist_sc_norm, hist_bulk_norm), axis=0)
                data_x.append(x_data)
                data_y.append(label)
        if len(data_x) > 0:
            data_xx = np.array(data_x, dtype='float32')[:, :, :, np.newaxis]
            #data_xx = np.array(data_x, dtype='float32')
            np.save(save_dir+'/xdata_tf_'+str(i)+'.npy', data_xx)
            np.save(save_dir+'/ydata_tf_'+str(i)+'.npy', data_y)


def main():
    sc_file = 'E-GEOD-141730.TPM.txt'
    bulk_file = 'GSE80744.TPM.txt'
    tf_target_file = 'ath_gene_pairs.txt'
    sep_file = 'ath_sep_list.txt'
    save_dir = "data"

    sep_index = get_sep_index(sep_file)
    tf_target_list = get_tf_target_list(tf_target_file)
    sc_data, sc_index = get_gene_expression_data(sc_file)
    bulk_data, bulk_index = get_gene_expression_data(bulk_file)
    gene_train_data(sc_data, sc_index, bulk_data, bulk_index, sep_index, tf_target_list, save_dir)


if __name__ == '__main__':
    main()

