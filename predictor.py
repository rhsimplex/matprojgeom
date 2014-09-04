import modelbuilder
import databasegenerator
import pandas as pd
import sys
import pymatgen as pm

def build_predictors(structure_csv='structures.csv'):
    column_headers = pd.read_csv(structure_csv).columns

    predictor_names = column_headers[1:20]

    target_pointgroup = column_headers[21]
    target_volumepersite = column_headers[23]
    target_c_a = column_headers[25]
    
    print 'Building point group classifier...',
    sys.stdout.flush()
    rfc, _, acc, le = modelbuilder.buildTreeClassifier(predictor_names, targetcolumn = target_pointgroup)
    print 'done. Accuracy: ' + str(acc)
    
    print 'Building volume/site regressor...',
    sys.stdout.flush()
    rfr_volume, _, mae = modelbuilder.buildTreeRegressor(predictor_names, targetcolumn = target_volumepersite)
    print 'done. MAE: ' + str(mae) + ' A^3'

    print 'Building c/a regressor...',
    sys.stdout.flush()
    rfr_c_a, _, mae = modelbuilder.buildTreeRegressor(predictor_names, targetcolumn = target_c_a)
    print 'done. MAE: ' + str(mae)

    return rfc, le, rfr_volume, rfr_c_a

def predictable_row(a, structure_csv='structures.csv'):
    column_headers = pd.read_csv(structure_csv).columns
    predictor_names = column_headers[1:20]
    row = []
    for header in predictor_names:
        row.append(eval('databasegenerator.' + header + '(a)'))
    print '=====Compositional Features====='
    print pd.Series(row, index = predictor_names)
    return row

def main():
    print 'Ryan\'s structure predictor 0.1!'
    sys.stdout.flush()
    rfc, le, rfr_volume, rfr_c_a = build_predictors()

    while True: 
        formula = raw_input('Enter composition to predict, or \'q\' to quit. Please include numbers for all elements (e.g. Li21Si5, Al2Ti1): ')
        if formula == 'q':
            break
        element = ''
        elements = []
        number = ''
        numbers = []
        for char in formula:
            if char.isalpha():
                if len(number) > 0:
                    numbers.append(int(number))
                    number=''
                element = element + char
            elif char.isdigit():
                if len(element) > 0:
                    elements.append(element)
                    element = ''
                number = number + char
        if len(number) > 0:
            numbers.append(int(number))
        if len(elements) != len(numbers):
            print 'Please enter numbers for all elements.'
        
        lattice = pm.Lattice.cubic(5.0)
        el_list = []
        for pair in zip(elements, numbers):
            el_list.extend([pair[0] for i in range(pair[1])])
        coords = [[0,0,0] for i in range(len(el_list))]
        try:
            a = pm.Structure(lattice, el_list, coords)
            row = predictable_row(a)
            print '=====Structure Predictions====='
            probs = rfc.predict_proba(row)
            pg_names = map(le.inverse_transform, range(len(probs[0])))
            sgs = pd.Series(probs[0], index = pg_names)
            sgs.sort(ascending=False)
            sgs = pd.DataFrame(sgs)
            sgs.index.name = 'Point Group'
            sgs.columns = ['Probability']
            print sgs[:5]
            print 'Volume/site: ' + str(round(rfr_volume.predict(row) ,2)) + ' A^3'
            print 'c/a: ' + str(round(rfr_c_a.predict(row), 2))
            print '==============================='
            
        except ValueError:
            print 'Invalid formula.'
         
