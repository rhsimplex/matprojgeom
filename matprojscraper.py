import pymatgen as pm
import os

def scraper(dest_folder,  API_key, resume=0, verbose=False):
    """
    Helper method for scraping database entries from the Materials Project database, in order of Materials Project ID.  Automatically resumes from largest ID in the destination folder.
    
    Keyword arguments:
    dest_folder -- path to store structure files
    API_key -- Materials project API key
    resume -- Materials Project ID to resume scraping from. I.e. 1234 will start with ID 'mp-1234', unless the dest_folder contains a file 'mp-XXXX' where XXXX is an integer > 1234.
    verbose -- print messages
    """
    last_materialsproject_id = 1
    
    # Create folder if none exists
    if not os.path.exists(dest_folder):
        if verbose: print 'No folder, creating ' + dest_folder
        os.mkdir(dest_folder)
    files_list = os.listdir(dest_folder)
    
    # Find last file by Materials Project ID
    if len(files_list) > 0 and resume > 0:
        ids = map(lambda x: int(x.replace('mp-','').replace('.mson','')), files_list)
        last_materialsproject_id = sorted(ids)[-1]
        if resume > last_materialsproject_id:
            last_materialsproject_id = resume
    
    #Creat MPRester object
    m = pm.MPRester(API_key)
    
    #Loop indefinitely querying database for entries, in order.  If entry doesn't exist, move to next ID. Structure objects saved in Materials Project '.mson' format
    while True:
        try:
            mpname = 'mp-' + str(last_materialsproject_id)
            a = m.get_structure_by_material_id(mpname)
            a.write_to_json_file(dest_folder + '//' + mpname + '.mson')
            if verbose: print mpname + ' ' + a.composition.alphabetical_formula + ' written.'
        except IndexError:
            if verbose: print mpname + ' not found in database.'
        last_materialsproject_id += 1
