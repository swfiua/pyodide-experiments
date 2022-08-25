import matplotlib
import numpy as np
matplotlib.use("module://matplotlib.backends.html5_canvas_backend")
#import matplotlib.cm as cm
#from matplotlib import pyplot as plt


print("HELLO WORLD!")
#from astroquery.mast import Observation
#from astroquery.simbad import Simbad

from gotu import jwst


from blume import magic, farm

fm = farm.Farm()

fm.add(jwst.Jwst())

await fm.start()

await fm.run()

