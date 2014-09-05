matprojgeom
===========
matprojgeom is a very simple crystal structure predictor. It uses some machine learning (random forests) on data from around 20,000 compounds from the Materials Project (https://materialsproject.org) to make structural predictions for a given composition.

To try it out, clone this repository and run predictor.py (you will need the pymatgen (https://pymatgen.org) library):
```
$python predictor.py
```
or if using ipython:
```
In[1]: import predictor
In[2]: predictor.main()
```
You'll see something like this:
```
```
