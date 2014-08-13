from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.cross_validation import train_test_split, cross_val_score
from sklearn.metrics import confusion_matrix, mean_absolute_error, make_scorer
from numpy import mean
import pandas as pd
import pymatgen as pm
import databasegenerator

def buildTreeClassifier(predictorColumns, structurestable = 'structures.csv',  targetcolumn = 'pointGroup', md = None):
    """
    Build a random forest-classifier model to predict some structure feature from compositional data.  Will return the model trained on all data, a confusion matrix calculated , and an average accuracy score. Also returns a label encoder object
    """
    df = pd.read_csv(structurestable)
    df = df.dropna()
    if('fracNobleGas' in df.columns):
        df = df[df['fracNobleGas'] <= 0]
    
    s = StandardScaler()
    le = LabelEncoder()
    
    X = s.fit_transform(df[predictorColumns])
    y = le.fit_transform(df[targetcolumn].values)

    rfc = RandomForestClassifier(max_depth = md)
    acc = mean(cross_val_score(rfc, X, y))

    X_train, X_test, y_train, y_test = train_test_split(X,y)
    rfc.fit(X_train,y_train)
    y_predict = rfc.predict(X_test)
    cm = confusion_matrix(y_test, y_predict)
    
    cm = pd.DataFrame(cm, columns=le.classes_, index=le.classes_)

    rfc.fit(X, y)

    return rfc, cm, acc, le

def buildTreeRegressor(predictorColumns, structurestable = 'structures.csv',  targetcolumn = 'c_a', md = None):
    """
    Build a random forest-regressor model to predict some structure feature from compositional data.  Will return the model trained on all data, a mean_absolute_error score, and a table of true vs. predicted values
    """
    df = pd.read_csv(structurestable)
    df = df.dropna()
    if('fracNobleGas' in df.columns):
        df = df[df['fracNobleGas'] <= 0]
    
    s = StandardScaler()
    
    X = s.fit_transform(df[predictorColumns])
    y = df[targetcolumn].values

    rfr = RandomForestRegressor(max_depth = md)
    acc = mean(cross_val_score(rfr, X, y, scoring=make_scorer(mean_absolute_error)))

    X_train, X_test, y_train, y_test = train_test_split(X,y)
    rfr.fit(X_train,y_train)
    y_predict = rfr.predict(X_test)
    
    t = pd.DataFrame({'True':y_test, 'Predicted':y_predict})
    
    rfr.fit(X, y)

    return rfr, t, acc

def buildCoordinationTreeRegressor(predictorColumns, element, coordinationDir = 'coordination/', md = None):
    """
    Build a coordination predictor for a given element from compositional structure data of structures containing that element. Will return a model trained on all data, a mean_absolute_error score, and a table of true vs. predicted values
    """
    df = pd.read_csv(coordinationDir + element + '.csv')
    df = df.dropna()
    if('fracNobleGas' in df.columns):
        df = df[df['fracNobleGas'] <= 0]
    
    s = StandardScaler()
    
    X = s.fit_transform(df[predictorColumns])
    y = df['avgCoordination'].values

    rfr = RandomForestRegressor(max_depth = md)
    acc = mean(cross_val_score(rfr, X, y, scoring=make_scorer(mean_absolute_error)))

    X_train, X_test, y_train, y_test = train_test_split(X,y)
    rfr.fit(X_train,y_train)
    y_predict = rfr.predict(X_test)
    
    t = pd.DataFrame({'True':y_test, 'Predicted':y_predict})
    
    rfr.fit(X, y)

    return rfr, t, acc
