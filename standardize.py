import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits, ascii
from photometry import do_aperture_photometry
from photutils.detection import DAOStarFinder
from astropy.stats import sigma_clipped_stats
import os

def standardize_star(Vmag, 
                     Bmag, 
                     v, 
                     b,
                     ):
    
    v0=[]
    b0=[]


    for i in range(0,len(v)):
        v0.append(v[i])
        b0.append(b[i])

    b_v = []

    V_v = []

    B_b = []

    for i in range(0,len(v0)):
        b_v.append(b0[i]-v0[i])

        V_v.append(Vmag[i]-v0[i])
        B_b.append(Bmag[i]-b0[i])


    fig1, ax1 = plt.subplots()

    ax1.set_xlabel('b0-v0')
    ax1.set_ylabel('V-v0')
    ax1.set_title('V Standardization')

    ax1.scatter(b_v, V_v)

    m1, b1 = np.polyfit(b_v, V_v,1)

    ax1.plot(b_v, m1*np.array(b_v) + b1, linestyle = '--', label = f'(V-v0) = {m1:.4f}(b0-v0) + {b1:.4f}')

    ax1.legend(loc = 1)

    fig1.savefig('V Standardization.png')

    fig2, ax2 =plt.subplots()

    ax2.set_xlabel('b0-v0')
    ax2.set_ylabel('B-b0')
    ax2.set_title('B Standardization')

    ax2.scatter(b_v, B_b)

    m2, b2 = np.polyfit(b_v, B_b,1)

    ax2.plot(b_v, m2*np.array(b_v) + b2, linestyle = '--', label = f'(B-b0) = {m2:.4f}(b0-v0) + {b2:.4f}')

    ax2.legend(loc = 1)

    fig2.savefig('B Standardization.png')
    return m1, m2, b1, b2


def do_standardization(file, transformations):

    transforms = ascii.read(transformations)
    cwd = os.getcwd()
    path = os.path.abspath(file)
    os.chdir(path)
    images = os.listdir(path)
    

    kv = transforms['Ext'][0]
    kb = transforms['Ext'][1]
    mu_v = transforms['Mu'][0]
    mu_b = transforms['Mu'][1]
    Cv = transforms['C'][0]
    Cb = transforms['C'][1]


    t_exp = []
    X =[]
    mags = []
    filter = []
    r = [12]
    for i in range(0, len(images)):
        ap_phot = None
        reduced_image = fits.getdata(images[i])
        header = fits.getheader(images[i])
        t_exp.append(header['EXPTIME'])
        X.append(header['AIRMASS'])
        filter.append(header['FILTER'])

        mean, median, std = sigma_clipped_stats(reduced_image)
            

        dao = DAOStarFinder(fwhm = 3, threshold = 1.04*std)


        dets = dao(reduced_image)
        dets = dets[((dets['xcentroid'] > 685) & (dets['xcentroid'] < 695) & (dets['ycentroid'] > 545) & (dets['ycentroid'] < 560))]

        x_pos = dets['xcentroid'].tolist()
        y_pos = dets['ycentroid'].tolist()

        positions = list(zip(x_pos, y_pos))

        ap_phot = do_aperture_photometry(images[i], positions, r, 100, 5);
        mags.append(-2.5*np.log10(ap_phot['aperture_sum_0'][0]/t_exp[i]))

    os.chdir(cwd)
    for i in range(0,len(mags)):
        if filter[i] == 'V':
            v = mags[i]
            Xv = X[i]
        else:
            b = mags[i]
            Xb = X[i]
    
    v0 = v - (kv*Xv)

    b0 = b - (kb*Xb)

    Vmag = (v0 - mu_v*(b0-v0) + Cv)

    Bmag = (b0 - mu_b*(b0-v0) + Cb)

    print(f'The magnitude of the target star in the V-band is {Vmag}')
    print(f'The magnitude of the target star in the B-band is {Bmag}')

    return Vmag, Bmag

    





    






    
    


    
