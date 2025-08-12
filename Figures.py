from astropy.visualization import ImageNormalize, ZScaleInterval, LinearStretch
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy.table import Table
from astropy.timeseries import TimeSeries
import numpy as np
from matplotlib.patches import Annulus

def before_after(image1, image2, filename):
    raw_image = fits.getdata(image1)
    science_image = fits.getdata(image2)


    norm1 = ImageNormalize(raw_image, interval=ZScaleInterval(), stretch=LinearStretch())
    norm2 = ImageNormalize(science_image, interval=ZScaleInterval(), stretch=LinearStretch())

    fig, ax = plt.subplot_mosaic(
        '''
        BA
        ''',
        figsize = (12, 4), 
        constrained_layout = True
    )

    # Plot at ['A']

    im1 = ax['A'].imshow(science_image, origin='lower', norm=norm1, cmap='YlOrBr_r')
    fig.colorbar(im1, ax = ax['A'])

    # Plot at ['B']

    im2 = ax['B'].imshow(raw_image, origin='lower', norm=norm1, cmap='YlOrBr_r')
    fig.colorbar(im1, ax = ax['B'])
    fig.savefig(filename)
    plt.show()

    return

def plot_image(image, filename):

    figure = fits.getdata(image)
    norm = ImageNormalize(figure, interval=ZScaleInterval(), stretch=LinearStretch())

    fig, ax = plt.subplots()
    plt.imshow(figure, origin='lower', norm=norm, cmap='YlOrBr_r')
    circle = Annulus((689.8, 555.7), r=12, width = 2 ,color='red', alpha=1)
    ax.add_patch(circle)

    plt.savefig(filename)
    return

def make_table():
    table = ascii.read('AE Ursa Majoris Time Series.csv', format = 'csv')

    ts = TimeSeries(table)

    ts['flux'] = np.around(ts['flux'], decimals = 3)

    ts['uncertainty'] = np.around(ts['uncertainty'], decimals = 3)


    latex = ascii.write(ts, format = 'latex')

    return latex