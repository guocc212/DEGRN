import numpy as np
import pandas as pd
import os, time
import keras
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc
from keras.models import Sequential,model_from_json
from keras.callbacks import EarlyStopping,ModelCheckpoint
from keras.layers import Dense, Dropout, Activation, Flatten, Conv2D, MaxPooling1D, MaxPooling2D
from keras.optimizers import SGD, Adam
import matplotlib.pyplot as plt


def get_sep_index(file_name):
    sep_list = []
    files = open(file_name, 'r')
    for line in files:
        sep_list.append(int(line))
    return np.array(sep_list)


def get_train_test_data(sep_index, data_dir):
    data_x_list = []
    data_y_list = []
    for i in sep_index:
        x_file = data_dir + '/xdata_tf_'+str(i)+'.npy'
        if os.path.exists(x_file):
            data_x = np.load(data_dir + '/xdata_tf_'+str(i)+'.npy', allow_pickle=True)
            data_y = np.load(data_dir + '/ydata_tf_' + str(i) + '.npy', allow_pickle=True)
            data_x_list.append(data_x)
            data_y_list.append(data_y)
    data_x_list = np.concatenate(data_x_list).astype("float32")
    data_y_list = np.concatenate(data_y_list).astype("int")
    x_train, x_test, y_train, y_test = train_test_split(data_x_list, data_y_list, train_size=0.7, random_state=2022)
    return x_train, x_test, y_train, y_test

# parameters setting
batch_size = 32
epochs = 150
sep_file = 'ath_sep_list.txt'
data_dir = "data"
out_dir = "models_bulk_sc"
type = 'whole'
roc_file = '{}/{}_roc.file.txt'.format(out_dir, type)

# main
rc = open(roc_file, 'w')
for rep in range(1, 6):
    sep_index = get_sep_index(sep_file)
    x_train, x_test, y_train, y_test = get_train_test_data(sep_index, data_dir)
    model = Sequential()
    model.add(
        Conv2D(filters=32,
               kernel_size=(3, 3),
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
    checkpoint = ModelCheckpoint(filepath='{}/whole_data_weights.hdf5'.format(out_dir),
                                 monitor='val_acc', verbose=1,save_best_only=True, mode='auto', period=1)
    res=model.fit(x_train,y_train,batch_size=batch_size,epochs=epochs,validation_split=0.2,
                  verbose=1,callbacks=[checkpoint, early_stopping])
    model.save('{}/{}_rep{}_model.h5'.format(out_dir, type, rep))

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
    graph_data.to_csv('{}/{}_rep{}_acc_loss.csv'.format(out_dir, type, rep), index=False, header=True)

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
    plt.plot(epoch_count, acc)
    plt.plot(epoch_count, val_acc)
    plt.legend(['Training Accuracy', 'Validation Accuracy'])
    plt.xlabel('Epoch')
    plt.ylabel('Accuracy')
    plt.title('Model Accuracy')
    plt.grid()
    plt.legend(['train', 'val'], loc='upper left')
    plt.savefig('{}/{}_rep{}_acc_los.pdf'.format(out_dir, type, rep))

    # plot AUC
    fig = plt.figure(figsize=(5, 5))
    roc_auc = auc(fpr, tpr)
    plt.plot(fpr, tpr, color='darkorange', lw=2, label='ROC curve (area = %0.2f)' % roc_auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0, 1])
    plt.ylim([0, 1])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic example')
    plt.legend(loc="lower right")
    plt.grid()
    plt.savefig('{}/{}_rep{}_ROC.pdf'.format(out_dir, type, rep))
    roc_data = pd.DataFrame(columns=['fpr', 'tpr'])
    roc_data['fpr'] = fpr
    roc_data['tpr'] = tpr
    roc_data.to_csv('{}/{}_rep{}_roc.csv'.format(out_dir, type, rep), index=False, header=True)
    rc.write('{}\t{}\t{}\n'.format(type, rep, roc_auc))

rc.close()