from utils import *
from scipy.optimize import curve_fit
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from sklearn.metrics import explained_variance_score
from sklearn import preprocessing

### Calculate k values ###
def exponential(t, a, b, c):
    return a * (1-np.exp(-b * t)) + c

def sigmoid(x, L ,x0, k, b):
    y = L / (1 + np.exp(-k*(x-x0)))+b
    return (y)

def calcKValues(data, func, t, maxLimit = 20000, normalise = True, plot = False):
  if plot:
    fig, ax = plt.subplots((len(data))//3+1, 3, sharex='col', sharey='row', squeeze = False, figsize=(30,90))

  k_values = []
  errors = []
  r2 = []
  for i in range(len(data)):
    #print(i)
    y = data.iloc[i]

    if (i+1)%3 == 0 or np.max(y) < maxLimit:
      k_values.append(0)
      continue

    if normalise:
      y = y/np.max(y)

    try:
      if func == 'exponential':
        popt, pcov = curve_fit(exponential, t.astype('float64'), y.astype('float64'), p0 = [np.max(y), 1, np.min(y)], maxfev=20000)
        k_values.append(popt[1])
        errors.append(np.mean((y.astype('float64')-exponential(t.astype('float64'), *popt))**2))
        ss_res = np.dot((y - exponential(t, *popt)),(y - exponential(t, *popt)))
        ymean = np.mean(y)
        ss_tot = np.dot((y-ymean),(y-ymean))
        r2.append(1-ss_res/ss_tot)

      if func == 'sigmoid':
        popt, pcov = curve_fit(sigmoid, t.astype('float64'), y.astype('float64'), p0 = [np.max(y), np.median(t), 1, np.min(y)], maxfev=20000)
        k_values.append(popt[2])
        errors.append(np.mean((y.astype('float64')-sigmoid(t.astype('float64'), *popt))**2))
        ss_res = np.dot((y - sigmoid(t, *popt)),(y - sigmoid(t, *popt)))
        ymean = np.mean(y)
        ss_tot = np.dot((y-ymean),(y-ymean))
        r2.append(1-ss_res/ss_tot)
      
    except:
      print(i)
      k_values.append(None)

    if plot:  
      ax[i//3, i%3].plot(t, y, 'ko', label="Original Noised Data")
      if func == 'exponential':
        ax[i//3, i%3].plot(t, exponential(t, *popt), 'r-', label="Fitted Curve")
      if func == 'sigmoid':
        ax[i//3, i%3].plot(t, sigmoid(t, *popt), 'r-', label="Fitted Curve")
        
  return k_values, np.mean(errors), np.mean(r2)

def calcAverage(k_values, lowk = 0.05, highk = 120):
  averaged_kvalues = []
  for i in range(len(k_values)):
    if (i+1)%3 == 0:
      continue
    if i%3 == 0:
      #Anomaly detection (if very low or very high k values, just ignore it and take other value)
      if k_values[i] < lowk or k_values[i] > highk:
        k_values[i] = k_values[i+1]
      if k_values[i+1] < lowk or k_values[i+1] > highk:
        k_values[i+1] = k_values[i]
      averaged_kvalues.append((k_values[i]+k_values[i+1])/2)
  return averaged_kvalues

### Model creation and training ###
def createRForrest(records_gRNA,records_target,kvalues,motifs_to_look_for,seed,numTrees = 2500, treeMaxDepth = 8):
    regr = RandomForestRegressor(max_depth=treeMaxDepth, random_state=seed, n_estimators = numTrees)
    feature_names = []
    X = []
    y = np.append(kvalues)

    for i in range(len(records_gRNA)):
        k = align(records_gRNA[i].seq,records_target[i].seq)
        X.append(np.concatenate((one_hot_encode(records_gRNA[i].seq), one_hot_encode(records_gRNA[i].seq), one_hot_encode(records_target[i].seq[:k]), one_hot_encode(records_target[i].seq[k+20:]), one_hot_encode_duplets(records_gRNA[i].seq)
                            ,calcGC(records_gRNA[i].seq), calcGC(records_target[i].seq), countMotifs(records_gRNA[i].seq, motifs_to_look_for)), axis = None).ravel())
    
    X = np.array(X)
    y = np.array(y)
    mask = np.random.rand(len(y)) < 0.75
    regr.fit(X[mask], y[mask])

    map = dict(zip(range(4), "acgt"))
    for i in range(80):
        feature_names.append('grnaposition{}, {}'.format(i//20+1, map[i%4]))
    for i in range(80):
        feature_names.append('targetalignedaposition{}, {}'.format(i//20+1, map[i%4]))
    for i in range(80):
        feature_names.append('targetposition{}, {}'.format(i//20+1, map[i%4]))
    for i in range(80):
        feature_names.append('targetposition{}, {}'.format(i//20+41, map[i%4]))
        map = dict(zip(range(16), ('aa', 'ac', 'ag', 'at', 'ca', 'cc', 'cg', 'ct', 'ga', 'gc', 'gg', 'gt', 'ta', 'tc', 'tg', 'tt')))
    for i in range(304):
        feature_names.append('grnaposition{}, {}'.format(i//19+1, map[i%16]))
        feature_names.append('gc_grna')
        feature_names.append('gc_target')
    for motif in motifs_to_look_for:
        feature_names.append('Motif: ' + motif)

def trainSVMModel(X_reduced, kvalues,testsize,seed,kernel = 'linear'): #Kernel rbf for ss or ds
    X_train, X_test, y_train, y_test = train_test_split(X_reduced, kvalues, test_size=testsize,random_state=seed)

    #Create a svm Classifier
    clf = make_pipeline(StandardScaler(), SVR(kernel = kernel))

    #Train the model using the training sets
    # clf.fit(X_train, encoded_ytrain)
    clf.fit(X_train, y_train)

    #Predict the response for test dataset
    y_pred = clf.predict(X_reduced)

    for i in range(len(y_pred)):
        if y_pred[i] < 0:
            y_pred[i] = 0
        
    return clf, explained_variance_score(kvalues, y_pred)

