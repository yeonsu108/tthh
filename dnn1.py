# ssh gpu-0-X ; conda activate py36

import os
import sys
import time
#os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"] = "1"
import uproot
import pandas as pd
import numpy as np
import tensorflow as tf
import pickle
from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras.models import load_model
from utils.plots import *
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split

start_time = time.time()
###################################################
#                     I/O                         #
###################################################
indir = "./samples2/"; PRE = "Test1105_bJet" # This sample should be apart from one for event classification.
outdir = "./DNN_result/" + PRE + "/bJetCassification/bCat_higgs_2/" # Higgs_2 bJet Categorization. 
os.makedirs(outdir, exist_ok=True)
class_names = ["Cat1","Cat2", "Cat3", "NoCat"]

inputvars = [
     "bJet_size",
     "bJet1_pt", "bJet1_eta", "bJet1_phi", "bJet1_mass",
     "bJet2_pt", "bJet2_eta", "bJet2_phi", "bJet2_mass",
     "bJet3_pt", "bJet3_eta", "bJet3_phi", "bJet3_mass",
     "bJet4_pt", "bJet4_eta", "bJet4_phi", "bJet4_mass",
     "b1b2_dr", "b1b3_dr", "b1b4_dr", "b2b3_dr", "b2b4_dr", "b3b4_dr", 

#     "Lep_size",
#     "Lep1_pt", "Lep1_eta", "Lep1_phi", "Lep1_t",
#     "Lep2_pt", "Lep2_eta", "Lep2_phi", "Lep2_t",

#     "bb_avg_dr", "bb_max_dr", "bb_min_dr", "b_ht", "bb_dEta_WhenMaxdR", "b_cent", "bb_max_deta", "bb_max_mass", "bb_twist",
#     "close_Higgs_pt", "close_Higgs_eta", "close_Higgs_phi", "close_Higgs_mass"
            ]

catvar = ["bCat_higgs_2"]
openvars = inputvars + catvar

###################################################
#                 PreProcessing                   #
###################################################
pd_data = uproot.open(indir+PRE+"_tthh_di.root")["Delphes"].arrays(openvars,library="pd")
pd_cat1 = pd_data.loc[pd_data["bCat_higgs_2"] == 0]
pd_cat2 = pd_data.loc[pd_data["bCat_higgs_2"] == 1]
pd_cat3 = pd_data.loc[pd_data["bCat_higgs_2"] == 2]
pd_cat4 = pd_data.loc[pd_data["bCat_higgs_2"] == 3]
#pd_cat5 = pd_data.loc[pd_data["bCat_higgs_1"] == 4]
#pd_cat6 = pd_data.loc[pd_data["bCat_higgs_1"] == 5]
#pd_cat7 = pd_data.loc[pd_data["bCat_higgs_1"] == 6]
#pd_cat8 = pd_data.loc[pd_data["bCat_top_1"] == 7]
#pd_cat9 = pd_data.loc[pd_data["bCat_top_1"] == 8]
#pd_cat10 = pd_data.loc[pd_data["bCat_top_1"] == 9]
#pd_cat11 = pd_data.loc[pd_data["bCat_top_1"] == 10]

nCat1 = len(pd_cat1)
nCat2 = len(pd_cat2)
nCat3 = len(pd_cat3)
nCat4 = len(pd_cat4)
#nCat5 = len(pd_cat5)
#nCat6 = len(pd_cat6)
#nCat7 = len(pd_cat7)
#nCat8 = len(pd_cat8)
#nCat9 = len(pd_cat9)
#nCat10 = len(pd_cat10)
#nCat11 = len(pd_cat11)
ntrain = min(nCat1, nCat2, nCat3, nCat4) #nCat5, nCat6, nCat7) #nCat8, nCat9, nCat10, nCat11)
print("ntrain = ", ntrain)

pd_cat1 = pd_cat1.sample(n=ntrain).reset_index(drop=True)
pd_cat2 = pd_cat2.sample(n=ntrain).reset_index(drop=True)
pd_cat3 = pd_cat3.sample(n=ntrain).reset_index(drop=True)
pd_cat4 = pd_cat4.sample(n=ntrain).reset_index(drop=True)
#pd_cat5 = pd_cat5.sample(n=ntrain).reset_index(drop=True)
#pd_cat6 = pd_cat6.sample(n=ntrain).reset_index(drop=True)
#pd_cat7 = pd_cat7.sample(n=ntrain).reset_index(drop=True)
#pd_cat8 = pd_cat8.sample(n=ntrain).reset_index(drop=True)
#pd_cat9 = pd_cat9.sample(n=ntrain).reset_index(drop=True)
#pd_cat10 = pd_cat10.sample(n=ntrain).reset_index(drop=True)
#pd_cat11 = pd_cat11.sample(n=ntrain).reset_index(drop=True)
print("pd_cat1", pd_cat1)

pd_data = pd.concat([pd_cat1, pd_cat2, pd_cat3, pd_cat4])#pd_cat5,pd_cat6, pd_cat7]) #pd_cat8, pd_cat9, pd_cat10, pd_cat11])
pd_data = pd_data.sample(frac=1).reset_index(drop=True)
print("pd_data", pd_data)
x_total = np.array(pd_data.filter(items = inputvars))
y_total = np.array(pd_data.filter(items = ['bCat_higgs_2'])) # Set your Target Column.

# Training and Cross-Validation Set
x_train, x_val, y_train, y_val = train_test_split(x_total, y_total, test_size=0.3)
print("x_train: ",len(x_train),"x_val: ", len(x_val),"y_train: ", len(y_train),"y_val", len(y_val))

###################################################
#                      Model                      #
###################################################
epochs = 1000; patience_epoch = 10; batch_size = 1024; print("batch size :", batch_size)
activation_function='relu'
weight_initializer = 'random_normal'
es = EarlyStopping(monitor='val_loss', mode='min', verbose=1, patience=patience_epoch)
mc = ModelCheckpoint(outdir+'/best_model.h5', monitor='val_loss', mode='min', save_best_only=True)

model = tf.keras.models.Sequential()
###############    Input Layer      ###############
model.add(tf.keras.layers.Flatten(input_shape = (x_train.shape[1],)))
###############    Hidden Layer     ###############
model.add(tf.keras.layers.BatchNormalization())
model.add(tf.keras.layers.Dense(30, activation=activation_function))
model.add(tf.keras.layers.BatchNormalization())
model.add(tf.keras.layers.Dense(30, activation=activation_function, kernel_regularizer='l2', kernel_initializer=weight_initializer))
###############    Output Layer     ###############
print("class_names : ", len(class_names))
model.add(tf.keras.layers.Dense(len(class_names), activation="softmax"))
###################################################

model.compile(optimizer=tf.keras.optimizers.Adam(clipvalue=0.5), loss="sparse_categorical_crossentropy", metrics = ["accuracy"])
model.summary()

hist = model.fit(x_train, y_train, batch_size=batch_size, epochs=epochs,
                                validation_data=(x_val,y_val), callbacks=[es, mc])

###################################################
#                  Prediction                     #
###################################################
print("#             PREDICTION                 #")
pred_train = model.predict(x_train); pred_train = np.argmax(pred_train, axis=1)
pred_val = model.predict(x_val) ; pred_val = np.argmax(pred_val, axis=1)
print("Is it similar?")
print("Prediction for validation set: ", pred_val)
print("Answer for train set:         ", y_val.T)

###################################################
#                Confusion Matrix                 #
###################################################
print("#           CONFUSION MATRIX             #")
plot_confusion_matrix(y_val, pred_val, classes=class_names,
                    title='Confusion matrix, without normalization', savename=outdir+"/confusion_matrix_val.pdf")
plot_confusion_matrix(y_val, pred_val, classes=class_names, normalize=True,
                    title='Normalized confusion matrix', savename=outdir+"/norm_confusion_matrix_val.pdf")
plot_confusion_matrix(y_train, pred_train, classes=class_names,
                    title='Confusion matrix, without normalization', savename=outdir+"/confusion_matrix_train.pdf")
plot_confusion_matrix(y_train, pred_train, classes=class_names, normalize=True,
                    title='Normalized confusion matrix', savename=outdir+"/norm_confusion_matrix_train.pdf")

###################################################
#                    Accuracy                     #
###################################################
print("#               ACCURACY                  #")
train_loss, train_acc = model.evaluate(x_train, y_train)
print(f"Train accuracy: {train_acc * 100:.2f}%")
test_loss, test_acc = model.evaluate(x_val, y_val)
print(f"Test accuracy: {test_acc * 100:.2f}%")

end_time = time.time()
execution_time = end_time - start_time
print(f"execution time: {execution_time} second")

###################################################
#              Feature Importance                 #
###################################################
print("#          FEATURE IMPORTANCE             #")
model_dir = outdir + '/best_model.h5'
model = tf.keras.models.load_model(model_dir)
model.summary()

input_data = tf.convert_to_tensor(x_val, dtype=tf.float32)
print(input_data)
name_inputvar = inputvars
n_evts = len(x_val)
n_var = len(name_inputvar)
mean_grads = n_var*[0.0]
all_grads = []
#mean_jacobian = np.zeros(len(name_inputvar))
#jacobian_matrix = np.zeros((len(name_inputvar),len(name_inputvar)))
for i, event in enumerate(x_val):
    print(i,"/",n_evts)
    with tf.GradientTape() as tape:
        inputs = tf.Variable([event])
        tape.watch(inputs)
        prediction = model(inputs)[:, 1]
    first_order_gradients = tape.gradient(prediction, inputs)
    gradiants = tf.abs(first_order_gradients)
    numpy_array = gradiants.numpy()[0]
    all_grads.append(numpy_array)
    #print(i,numpy_array , len(numpy_array))
    for n in range(len(mean_grads)):
        mean_grads[n] += abs(numpy_array[n])/n_evts

print(mean_grads) # This is just final iteration (final event), not important yet.
df = pd.DataFrame(all_grads)
print(df)
df.to_csv('data.csv', index=False)

gradiants = tf.abs(first_order_gradients)
numpy_array = gradiants.numpy()
df = pd.DataFrame(numpy_array)
print(df)
df.to_csv(outdir+'/data.csv', index=False)
feature_importance_first_order  = tf.reduce_mean(tf.abs(first_order_gradients), axis=0)
feature_importance_dict = dict(zip(name_inputvar, feature_importance_first_order.numpy())) # Mapping
#feature_importance_Secondorder  = tf.reduce_mean(tf.abs(second_order_gradients), axis=0)
feature_importance_series = pd.Series(feature_importance_dict)

print("Feature importance series?")
print(feature_importance_series)
max_importance_score = feature_importance_series.max()
print("Max: ", max_importance_score)
for i, importance_score in enumerate(feature_importance_first_order):
    print(f"Feature {i+1} , {name_inputvar[i]} Importance: {max_importance_score-importance_score:.5f}")

print(feature_importance_series.index, feature_importance_series.values)

# Order the Feature Importance
sorted_indices = np.argsort(feature_importance_series.values)[::-1]
sorted_importance = feature_importance_series.values[sorted_indices]
sorted_feature_names = feature_importance_series.index[sorted_indices]

plt.figure(figsize=(10, 10))
plt.barh(sorted_feature_names, max_importance_score - sorted_importance)
plt.xlabel('Feature Importance')
plt.ylabel('Feature Names')
plt.savefig(outdir+'/first_order_gradient_importance.png')
plt.title('First-Order Gradient Feature Importance')
plt.show()

print("---Done---")
