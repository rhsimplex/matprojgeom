import pymatgen as pm
import os

def scraper(dest_folder,  API_key='3QUoyXBlHChz5OLg', resume=True, verbose=False):
    last_materialsproject_id = 1
    if not os.path.exists(dest_folder):
        if verbose: print 'No folder, creating ' + dest_folder
        os.mkdir(dest_folder)
    files_list = os.listdir(dest_folder)
    if len(files_list) > 0 and resume:
        ids = map(lambda x: int(x.replace('mp-','').replace('.mson','')), files_list)
        last_materialsproject_id = sorted(ids)[-1]

    m = pm.MPRester(API_key)

    while True:
        try:
            mpname = 'mp-' + str(last_materialsproject_id)
            a = m.get_structure_by_material_id(mpname)
            a.write_to_json_file(dest_folder + '//' + mpname + '.mson')
            if verbose: print mpname + ' ' + a.composition.alphabetical_formula + ' written.'
        except IndexError:
            if verbose: print mpname + ' not found in database.'
        last_materialsproject_id += 1
