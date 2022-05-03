import numpy as np
import pandas as pd
import datetime
import os, time
import keras
from sklearn.model_selection import train_test_split
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve, auc
from keras.models import Sequential
from keras.callbacks import EarlyStopping,ModelCheckpoint
from keras.layers import Dense, Dropout, Activation, Flatten, Conv2D, MaxPooling2D
from keras.optimizers import SGD, Adam
import matplotlib.pyplot as plt


def get_sep_index(file_name):
    sep_list = []
    files = open(file_name, 'r')
    for line in files:
        sep_list.append(int(line))
    return np.array(sep_list)


def get_TF_data(sep_index, data_dir):
    data_x_list = []
    data_y_list = []
    for i in sep_index:
        x_file = data_dir + '/xdata_tf_'+str(i)+'.npy'
        y_file = data_dir + '/ydata_tf_'+str(i)+'.npy'
        if os.path.exists(x_file):
            data_x = np.load(data_dir + '/xdata_tf_'+str(i)+'.npy', allow_pickle=True)
            data_y = np.load(data_dir + '/ydata_tf_' + str(i) + '.npy', allow_pickle=True)
            data_x_list.append(data_x)
            data_y_list.append(data_y)
    data_x_list = np.concatenate(data_x_list).astype("float32")
    data_y_list = np.concatenate(data_y_list).astype("int")
    return data_x_list, data_y_list


# parameters setting
batch_size = 32
epochs = 150
sep_file = 'ath_sep_list.txt'
data_dir = 'data'
out_dir = 'models'
stat_file ='three_fold_cross_result.txt'

sep_index = get_sep_index(sep_file)
stats = open(stat_file, 'w')
xx, yy = get_TF_data(sep_index, data_dir)
k = np.floor(len(xx)/10).astype("int")
print(k)
print(len(xx))
print(len(yy))
print(xx.shape)

for rep in range(1, 21):
    for fold in range(1,11):
        test_data = [i for i in range((fold-1)*k, fold*k)]
        train_data = [i for i in range(0, len(xx)-1) if i not in test_data]
        x_test = xx[test_data]
        y_test = yy[test_data]
        x_train = xx[train_data]
        y_train = yy[train_data]
        print("test sample number: {}\n".format(x_test.shape))
        print("train sample number: {}\n".format(x_train.shape))

        model = Sequential()
        model.add(Conv2D(filters=32,kernel_size=(3, 3),
                padding='same',
                activation='relu',
                input_shape=x_train.shape[1:])
        )
        model.add(Conv2D(32, (3, 3), activation="relu", padding="same"))
        model.add(MaxPooling2D(pool_size=(2, 2)))
        model.add(Dropout(0.5))
        model.add(Conv2D(64, (3, 3), activation="relu", padding="same"))
        model.add(Conv2D(64, (3, 3), activation="relu", padding="same"))
        model.add(MaxPooling2D(pool_size=(2, 2)))
        model.add(Dropout(0.5))
        model.add(Conv2D(128, (3, 3), activation="relu", padding="same"))
        model.add(Conv2D(128, (3, 3), activation="relu", padding="same"))
        model.add(MaxPooling2D(pool_size=(2, 2)))
        model.add(Dropout(0.5))
        model.add(Flatten())
        model.add(Dropout(0.5))
        model.add(Dense(512, activation='relu'))
        model.add(Dropout(0.5))
        model.add(Dense(1, activation='sigmoid'))
        print(model.summary())
        sgd = SGD(lr=0.01, decay=1e-5, momentum=0.9, nesterov=True)
        model.compile(optimizer=sgd, loss='binary_crossentropy', metrics=['accuracy'])
        early_stopping = keras.callbacks.EarlyStopping(monitor='val_accuracy', patience=300, verbose=0, mode='auto')
        checkpoint = ModelCheckpoint(filepath='{}/fold_{}_weights.hdf5'.format(out_dir, str(fold)),
                                     monitor='val_acc', verbose=1,save_best_only=True, mode='auto', period=1)
        res=model.fit(x_train,y_train,batch_size=batch_size,epochs=epochs,validation_split=0.2,
                      verbose=1,callbacks=[checkpoint, early_stopping])
        model.save('{}/fold_{}_rep{}_model.h5'.format(out_dir, fold, rep))

        # AUC
        y_predict = model.predict(x_test)
        fpr = dict()
        tpr = dict()
        roc_auc = dict()
        fpr, tpr, _ = roc_curve(y_test.reshape(len(y_test), 1), y_predict)
        roc_auc = auc(fpr, tpr)
        training_loss = res.history['loss']
        test_loss = res.history['val_loss']
        acc = res.history['accuracy']
        val_acc = res.history['val_accuracy']
        # Create count of the number of epochs
        epoch_count = range(1, len(training_loss) + 1)

        # get each curve data
        graph_data = pd.DataFrame(columns=['acc', 'val_acc', 'loss', 'val_loss'])
        graph_data['acc'] = acc
        graph_data['val_acc'] = val_acc
        graph_data['loss'] = training_loss
        graph_data['val_loss'] = test_loss
        graph_data.to_csv('{}/fold_{}_rep{}_acc_loss.csv'.format(out_dir, fold, rep), index=False, header=True)

        # draw graph
        plt.figure(figsize=(10, 6))
        plt.subplot(121)
        plt.plot(epoch_count, training_loss)
        plt.plot(epoch_count, test_loss)
        plt.legend(['Training Loss', 'Validation Loss'])
        plt.xlabel('Epoch')
        plt.ylabel('Loss')
        plt.title('Model Loss')
        plt.grid()
        plt.legend(['train', 'val'], loc='upper left')
        plt.subplot(122)
        plt.plot(epoch_count, acc, 'r--')
        plt.plot(epoch_count, val_acc, 'b-')
        plt.legend(['Training Accuracy', 'Validation Accuracy'])
        plt.xlabel('Epoch')
        plt.ylabel('Accuracy')
        plt.title('Model Accuracy')
        plt.grid()
        plt.legend(['train', 'val'], loc='upper left')
        plt.savefig('{}/fold_{}_rep{}_acc_los.pdf'.format(out_dir, fold, rep))
        # AUC
        fig = plt.figure(figsize=(5, 5))
        roc_auc = auc(fpr, tpr)
        plt.plot(fpr, tpr, color='#24a8ac', lw=2, label='ROC curve (area = %0.2f)' % roc_auc)
        plt.plot([0, 1], [0, 1], color='gray', lw=2, linestyle='--')
        plt.xlim([0, 1])
        plt.ylim([0, 1])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('Receiver operating characteristic example')
        plt.legend(loc="lower right")
        plt.grid()
        plt.savefig('{}/fold_{}_rep{}_ROC.pdf'.format(out_dir, fold, rep))
        rocV = np.trapz(tpr, fpr)
        stats.write('{}\t{}\t{}\t{}\n'.format(rep, fold, roc_auc, rocV))
        roc_data = pd.DataFrame(columns=['fpr', 'tpr'])
        roc_data['fpr'] = fpr
        roc_data['tpr'] = tpr
        roc_data.to_csv('{}/fold_{}_rep{}_roc.csv'.format(out_dir, fold, rep), index=False, header=True)
        ##
stats.close()
