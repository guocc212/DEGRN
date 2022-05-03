import numpy as np
import pandas as pd
import os, time
import h5py
from keras.models import load_model


def get_gene_pairs(file_name):
    tf_target_list = []
    files = open(file_name, 'r')
    for line in files:
        index, tf, target = line.strip().split('\t')
        tf_target_list.append([index, tf, target])
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


def get_tf_target_data(bulk_data, bulk_index, sc_data, sc_bulk, start_index, end_index, tf_target_list):
    data_x = []
    data_tf = []
    data_target = []
    data_index = []
    for index, tf, target in tf_target_list[start_index: end_index]:
        sc_flag = 0
        bulk_flag = 0
        #print('[INFO]{}_{}_{}...{}\n'.format(index, tf, target, time.asctime()))
        if tf in sc_index.keys() and target in sc_index.keys():
            tf_sc = np.log10(sc_data[sc_index[tf]].reshape(-1)[2:] + 10 ** -2)
            target_sc = np.log10(sc_data[sc_index[target]].reshape(-1)[2:] + 10 ** -2)
            hist_sc = np.histogram2d(tf_sc, target_sc, bins=32)
            hist_sc_t = hist_sc[0].T
            hist_sc_norm = (np.log10(hist_sc_t / sc_data.shape[1] + 10 ** -4) + 4) / 4
            sc_flag = 1
        if tf in bulk_index.keys() and target in bulk_index.keys():
            tf_bulk = np.log10(bulk_data[bulk_index[tf]].reshape(-1)[2:] + 10 ** -2)
            target_bulk = np.log10(bulk_data[bulk_index[target]].reshape(-1)[2:] + 10 ** -2)
            hist_bulk = np.histogram2d(tf_bulk, target_bulk, bins=32)
            hist_bulk_t = hist_bulk[0].T
            hist_bulk_norm = (np.log10(hist_bulk_t / bulk_data.shape[1] + 10 ** -4) + 4) / 4
            bulk_flag = 1
        if sc_flag == 1 and bulk_flag == 1:
            x_data = np.concatenate((hist_sc_norm, hist_bulk_norm), axis=0)
            data_x.append(x_data)
            data_index.append(index)
            data_tf.append(tf)
            data_target.append(target)
    if len(data_x) >0:
        data_x = np.array(data_x, dtype='float32')[:,:,:,np.newaxis]
    return data_x, data_index, data_tf, data_target


# parameters setting
sc_file = 'E-GEOD-141730.TPM.txt'
bulk_file = 'GSE80744.TPM.txt'
model_dir = 'models'
k = 1000
batch_size = 32

# main
print('[INFO] Reading SCgene expression...{}\n'.format(time.asctime()))
sc_data, sc_index = get_gene_expression_data(sc_file)
print('[INFO] Reading Bulk gene expression...{}\n'.format(time.asctime()))
bulk_data, bulk_index = get_gene_expression_data(bulk_file)

model_file = model_dir+'/whole_data_model.h5'
models = load_model(model_file)
for p in range(202,224):
    t = 'p{}'.format(p)
    predict_file = 'predict/predict_gene_pairs_{}.txt'.format(t)
    print('[INFO] Reading gene pairs...{}\n'.format(time.asctime()))
    tf_target_list = get_gene_pairs(predict_file)
    print('[INFO] End...{}\n'.format(time.asctime()))
    f = open('{}/predict-{}-result.txt'.format(model_dir, t), 'w')

    for i in range(1, len(tf_target_list), k):
    #for i in range(1, k+1, k):
        start_index = i-1
        end_index = i+k-1
        print('[INFO] {}--{} begins...{}\n'.format(start_index, end_index, time.asctime()))
        data_x, data_index, data_tf, data_target = get_tf_target_data(bulk_data, bulk_index, sc_data, sc_index, start_index, end_index, tf_target_list)
        print('[INFO] {}--{} model predict...{}\n'.format(start_index, end_index, time.asctime()))
        if len(data_x) > 0 :
            y_predict = models.predict(data_x, batch_size=32)
            print('[INFO] {}--{} model end...{}\n'.format(start_index, end_index, time.asctime()))
            for index, tf, target, predict in zip(data_index, data_tf, data_target, y_predict):
                f.write('{}\t{}\t{}\t{}\n'.format(index, tf, target, predict.astype('|S10').tobytes().decode('utf-8')))
    f.close()
