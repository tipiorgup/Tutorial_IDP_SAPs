import os
from collections import Counter
import yaml
import pandas as pd
import numpy as np
import pickle
from sklearn.model_selection import train_test_split
from sklearn import svm
import matplotlib.pyplot as plt
import math
import qml
from qml import *
import seaborn as sns
import pandas as pd
from qml.representations import generate_slatm
import math
from scipy.interpolate import make_interp_spline



d={'A':1,'R':2,'N':3,'D':4,'C':5,'Q':6,'G':7,'E':8,'H':9,'I':10,'L':11,'K':12,'M':13,'F':14,'P':15,'S':16,'T':17,'W':18,'Y':19,'V':20}
# d = {ni: indi for indi, ni in enumerate(set(aa_list))}

numbers_n= [d[ni] for ni in list(d.keys())]
zs = []
zs.append(np.array(numbers_n))
mbtypes = qml.representations.get_slatm_mbtypes(zs)


def get_index(res1,res2):
    try:
        p=mbtypes.index([d[res1],d[res2]])
    except ValueError:
        try:
            p=mbtypes.index((d[res1],d[res2]))
        except ValueError:
            try:
                p=mbtypes.index([d[res2],d[res1]])
            except ValueError:
                p=mbtypes.index((d[res2],d[res1]))
    s=((p-20)*9)+p
    e=s+10
    return s,e

def get_index_a(res1,res2,res3):
    try:
        p=mbtypes.index([d[res1],d[res2],d[res3]])
    except ValueError:
        try:
            p=mbtypes.index([d[res1],d[res3],d[res2]])
        except ValueError:
            try:
                p=mbtypes.index([d[res2],d[res1],d[res3]])
            except ValueError:
                try:
                    p=mbtypes.index([d[res3],d[res1],d[res2]])
                except ValueError:
                    try:
                        p=mbtypes.index([d[res3],d[res2],d[res1]])
                    except ValueError:
                        p=mbtypes.index([d[res2],d[res3],d[res1]])
    s=((p-20)*9)+p
    e=s+10
    return s,e

def plot_hist(peptide_list,res1,res2):
    a=[word for word in peptide_list if res1 in word]
    sequnce_list=[word for word in a if res2 in word]
    s,e=get_index(res1,res2)
    for l in sequnce_list:
        x=np.arange(0,20,2)
        y=slatm[l][0][s:e]
        X_Y_Spline = make_interp_spline(x, y)
        X_ = np.linspace(x.min(), x.max(), 500)
        Y_ = X_Y_Spline(X_)
        plt.plot(X_, Y_,label=l)
    plt.title('Distance %s-%s'%(res1,res2))
    plt.xlabel(r'$\AA$')
    plt.ylabel(r'$\frac{1}{\sigma \sqrt{2\pi}}e^{\frac{x^2}{2\sigma}}$')
    plt.legend()
    plt.show()

def plot_hist_a(peptide_list,res1,res2,res3):
    a=[word for word in peptide_list if res1 in word]
    b=[word for word in a if res2 in word]
    sequnce_list=[word for word in b if res2 in word]
    s,e=get_index_a(res1,res2,res3)
    for l in sequnce_list:
        x=np.arange(0,360,360/10)
        y=slatm[l][0][s:e]
        X_Y_Spline = make_interp_spline(x, y)
        X_ = np.linspace(x.min(), x.max(), 500)
        Y_ = X_Y_Spline(X_)
        plt.plot(X_, Y_,label=l)
    plt.title('Angle %s-%s-%s'%(res1,res2,res3))
    plt.xlabel(r'$\theta$')
    plt.ylabel(r'$\frac{1}{\sigma \sqrt{2\pi}}e^{\frac{x^2}{2\sigma}}$')
    plt.legend()
    plt.show()
