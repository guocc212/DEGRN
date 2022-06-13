# This script is used for predict the user's customized files;
# python predict_tf_target_file_from_DEGRN.py sc_file bulk_file model_file predict_file out_file
# eg : python predict_tf_target_file_from_DEGRN.py examples/scRNA-seq.tpm.txt examples/bulk_RNA-seq.tpm.txt models/whole_data_model.h5 examples/examples.tf.target.txt examples/predict-result.txt
# sc_file: scRNA-Seq expression file
# bulk_file: RNA-Seq expression file
# model_file: model file of DEGRN, default as "models/whole_data_model.h5"
# predict_file: the customized of user, including three columns: 1. ID; 2. tf ID (AGI); 3. target ID (AGI)
# out_file: results of prediction by DEGRN, including four columns.
#          The firstly three columns are same with predict_file. And the fourth column was the score by DEGRN.

# import module
import numpy as np
import sys, time
from keras.models import load_model
##
error_info = "  python predict_tf_target_file_from_DEGRN.py sc_file bulk_file model_file predict_file out_file\n" \
             "  eg : python predict_tf_target_file_from_DEGRN.py ../E-GEOD-141730.TPM.txt ../GSE80744.TPM.txt models/whole_data_model.h5 examples/examples.tf.target.txt examples/predict-result.txt\n" \
             "    sc_file: scRNA-Seq expression file\n" \
             "    bulk_file: RNA-Seq expression file\n" \
             "    model_file: model file of DEGRN, default as 'models/whole_data_model.h5'\n" \
             "    predict_file: the customized of user, including three columns: 1. ID; 2. tf ID (AGI); 3. target ID (AGI)\n" \
             "    out_file: results of prediction by DEGRN, including four columns.\n" \
             "        The firstly three columns are same with predict_file. And the fourth column was the score by DEGRN.\n"

if len(sys.argv) < 6:
    print('[Error] Need paremeters')
    print(error_info)
    sys.exit(1)


# parameters
sc_file = sys.argv[1]
bulk_file = sys.argv[2]
if sys.argv[3] == "NA":
    model_file = 'models/whole_data_model.h5'
else:
    model_file = sys.argv[3]
predict_file = sys.argv[4]
out_file = sys.argv[5]

#########
# functions
def get_gene_pairs(file_name):
    # read and store the user's self customized file
    # the file contain three columns, including ID, the AGI ID of TF and target,
    # examples are as followed:
    #   1	AT1G01010	AT1G01020
    #   2	AT1G01010	AT1G01030
    #   3	AT1G01010	AT1G01040
    tf_target_list = []
    files = open(file_name, 'r')
    for line in files:
        index, tf, target = line.strip().split('\t')
        tf_target_list.append([index, tf, target])
    files.close()
    return np.array(tf_target_list)


def get_gene_expression_data(file_name):
    # read and store the file of expression data
    # the file needs the headers, containing the sample name...
    # the first column was the gene name
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


def get_tf_target_data(bulk_data, bulk_index, sc_data, sc_index, start_index, end_index, tf_target_list):
    # transform the input data for DEGRN

    # defined the lists of data
    data_x = []
    data_tf = []
    data_target = []
    data_index = []
    for index, tf, target in tf_target_list[start_index: end_index]:
        sc_flag = 0
        bulk_flag = 0
        # determine whether if tf and target expressed in scRNA-Seq
        if tf in sc_index.keys() and target in sc_index.keys():
            tf_sc = np.log10(sc_data[sc_index[tf]].reshape(-1)[2:] + 10 ** -2)
            target_sc = np.log10(sc_data[sc_index[target]].reshape(-1)[2:] + 10 ** -2)
            hist_sc = np.histogram2d(tf_sc, target_sc, bins=32)
            hist_sc_t = hist_sc[0].T
            hist_sc_norm = (np.log10(hist_sc_t / sc_data.shape[1] + 10 ** -4) + 4) / 4
            sc_flag = 1
        # determine whether if tf and target expressed in bulk RNA-Seq
        if tf in bulk_index.keys() and target in bulk_index.keys():
            tf_bulk = np.log10(bulk_data[bulk_index[tf]].reshape(-1)[2:] + 10 ** -2)
            target_bulk = np.log10(bulk_data[bulk_index[target]].reshape(-1)[2:] + 10 ** -2)
            hist_bulk = np.histogram2d(tf_bulk, target_bulk, bins=32)
            hist_bulk_t = hist_bulk[0].T
            hist_bulk_norm = (np.log10(hist_bulk_t / bulk_data.shape[1] + 10 ** -4) + 4) / 4
            bulk_flag = 1
        # if the expression of tf and target exist in both scRNA-Seq and bulk RNA-Seq, we combined the two matrixs into a large matrix.
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
#sc_file = '../E-GEOD-141730.TPM.txt'
#bulk_file = '../GSE80744.TPM.txt'
#out_file = 'examples/predict-result.txt'
#model_dir = 'models'
#model_file = model_dir+'/whole_data_model.h5'
#predict_file = "examples/examples.tf.target.txt"
#
k = 1           # split setting

# main
# Step 1. read the expression data
print('[INFO] Reading scRNA-Seq expression data...{}\n'.format(time.asctime()))
sc_data, sc_index = get_gene_expression_data(sc_file)
print('[INFO] Reading bulk RNA-Seq expression data...{}\n'.format(time.asctime()))
bulk_data, bulk_index = get_gene_expression_data(bulk_file)
# Step 2. load model from DEGRN
models = load_model(model_file)
print('[INFO] Reading gene pairs...{}\n'.format(time.asctime()))
tf_target_list = get_gene_pairs(predict_file)
print('[INFO] End...{}\n'.format(time.asctime()))
out = open(out_file, 'w')

for i in range(1, len(tf_target_list), 1):
    start_index = i-1
    end_index = i+k-1
    print('[INFO] {}--{} begins...{}\n'.format(start_index, end_index, time.asctime()))
    data_x, data_index, data_tf, data_target = get_tf_target_data(bulk_data, bulk_index, sc_data, sc_index, start_index, end_index, tf_target_list)
    print('[INFO] {}--{} model predict...{}\n'.format(start_index, end_index, time.asctime()))
    if len(data_x) > 0 :
        y_predict = models.predict(data_x, batch_size=32)
        print('[INFO] {}--{} model end...{}\n'.format(start_index, end_index, time.asctime()))
        for index, tf, target, predict in zip(data_index, data_tf, data_target, y_predict):
            out.write('{}\t{}\t{}\t{}\n'.format(index, tf, target, predict.astype('|S10').tobytes().decode('utf-8')))
out.close()


