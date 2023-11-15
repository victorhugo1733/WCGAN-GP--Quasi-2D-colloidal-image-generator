from __future__ import division, unicode_literals, print_function  # for compatibility with Python 2 and 3

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas import DataFrame, Series  # for convenience
import pims
import trackpy as tp
from pathlib import Path
import pandas as pd


# Optionally, tweak styles.
mpl.rc('figure',  figsize=(10, 5))
mpl.rc('image', cmap='gray')

#abro carpeta de imagenes
frames = pims.open(r'D:\Images_Polimero0 NaCl0\Victor\Victor\Hugo\Hugo\800x800\Generadas\0.101\*.png')
radius=1


iterations=2000
iteration=0
plot_properties = dict(markeredgewidth=3, markersize=radius, markerfacecolor='none')

f_locate= tp.locate(frames[iteration],5 , minmass=200, percentile=95 ,separation=5)
   

plt.figure()
plt.imshow(frames[iteration])
plt.plot(f_locate['x'], f_locate['y'], 'x', markeredgecolor='red', **plot_properties)
plt.xlim(0, 600);
plt.ylim(0, 600);

for iteration in range(iterations):
    f = tp.locate(frames[iteration],5 , minmass=200,percentile=95 ,separation=5)
    tp.annotate(f, frames[iteration], plot_style={'markersize': radius}); #dibuja
    #Data   'D:\Images_Polimero0 NaCl0\Victor\Victor\Hugo\Hugo\600x600\4\frames\40\frame{i}.csv' valio keso
    filepath = Path(r'D:\Images_Polimero0 NaCl0\Victor\Victor\Hugo\Hugo\800x800\Generadas\frames\0.101\frame{i}.csv'.format(i=iteration) )
    df = pd.DataFrame(data=f)
    filepath.parent.mkdir(parents=True, exist_ok=True) 
    df.to_csv(filepath)
    df=df.drop(["mass", "size", "ecc", "signal", "raw_mass", "ep", "frame"],axis=1)
    df.to_csv(filepath)
      

    #Hago zoom
    plt.figure()
    tp.annotate(f, frames[iteration], plot_style={'markersize': radius}, ax=plt.gca())
    plt.ylim(0, 600)
    plt.xlim(0, 600);
    plt.savefig(r'D:\Images_Polimero0 NaCl0\Victor\Victor\Hugo\Hugo\800x800\Generadas\images_tracking\0.101\frame{i}.png'.format(i=iteration) ,dpi=150,bbox_inches='tight')
    plt.close()
    

