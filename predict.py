import sys
import dnaFeature
import numpy as np
import pickle
import pandas as pn
from sklearn import preprocessing

def predict_loc(pssm,spd3,option):
    if option == 0:
        model = "a.model"
        support = "a.sup"
        data = "a.csv"
    else:
        model = "b.model"
        support = "a.sup"
        data = "b.csv"
    features = dnaFeature.extraction(pssm,spd3)
    f= np.asarray(features)
    # scale data
    dataAll = pn.read_csv(data)
    dataMat = pn.DataFrame.as_matrix(dataAll)
    X = dataMat[:, :-1]
    y = dataMat[:, -1]
    scaler = preprocessing.StandardScaler().fit(X)
    f_scaled = scaler.transform(f.reshape(1,-1))
    # # step 2: reduce the features
    sup = pn.read_csv(support, header=None)
    sup = pn.DataFrame.as_matrix(sup)
    i = 0
    boolsup = []
    for s in sup:
        if s == 1.0:
            boolsup.append(True)
            i += 1
        else:
            boolsup.append(False)

    f_red = f_scaled.T[boolsup]
    # step 3: predict
    clf = pickle.load(open(model, "rb"))

    # step 4: return
    pred = clf.predict(f_red.reshape(1, -1))
    return pred

    return 0;

if(len(sys.argv)< 4):
    print("Not Enough Arguments Provided.")
    print("Usage:"+sys.argv[0]+ " PSSM SPD3 OPTION")
    print("PSSM: PSSM file generated from PSI-BLAST")
    print("SPD3: SPD3 file generated from SPIDER3")
    print("OPTION = 0 for identification of cellular phages")
    print("OPTION = 1 for prediction of phage location")
else:
    pssm = sys.argv[1]
    spd3 = sys.argv[2]
    option = int(sys.argv[3])
    if option!=1 and option !=0:
        print("Invalid OPTION parameter.")
        print("Usage:" + sys.argv[0] + " PSSM SPD3 OPTION")
        print("PSSM: PSSM file generated from PSI-BLAST")
        print("SPD3: SPD3 file generated from SPIDER3")
        print("OPTION = 0 for identification of cellular phages")
        print("OPTION = 1 for prediction of phage location")
    else:
        pred=predict_loc(pssm,spd3,option)
        if option==0:
            if pred == 1.0:
                print("Prediction: Phage located in Host Cell")
            else:
                print("Prediction: Phage is extra cellular")
        else:
            if pred == 1.0:
                print("Prediction: Phage located in Cell Membrane")
            else:
                print("Prediction: Phage located in CytoPlasm")
