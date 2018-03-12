import numpy as np
import segmentingFunctions as segfun
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#import time

# Prem Rachakonda (2017)
#
# This software was developed by employees of the National Institute of Standards and Technology (NIST), an agency of the Federal Government. Pursuant to title 17 United States Code Section 105, works of NIST employees are not subject to copyright protection in the United States and are considered to be in the public domain. Permission to freely use, copy, modify, and distribute this software and its documentation without fee is hereby granted, provided that this notice and disclaimer of warranty appears in all copies.
# THE SOFTWARE IS PROVIDED 'AS IS' WITHOUT ANY WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND FREEDOM FROM INFRINGEMENT, AND ANY WARRANTY THAT THE DOCUMENTATION WILL CONFORM TO THE SOFTWARE, OR ANY WARRANTY THAT THE SOFTWARE WILL BE ERROR FREE. IN NO EVENT SHALL NIST BE LIABLE FOR ANY DAMAGES, INCLUDING, BUT NOT LIMITED TO, DIRECT, INDIRECT, SPECIAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF, RESULTING FROM, OR IN ANY WAY CONNECTED WITH THIS SOFTWARE, WHETHER OR NOT BASED UPON WARRANTY, CONTRACT, TORT, OR OTHERWISE, WHETHER OR NOT INJURY WAS SUSTAINED BY PERSONS OR PROPERTY OR OTHERWISE, AND WHETHER OR NOT LOSS WAS SUSTAINED FROM, OR AROSE OUT OF THE RESULTS OF, OR USE OF, THE SOFTWARE OR SERVICES PROVIDED HEREUNDER.
# Distributions of NIST software should also include copyright and licensing statements of any third-party software that are legally bundled with the code in compliance with the conditions of those licenses.


plt.close("all")
    
trueRadius = 0.05; # This is the nominal diameter of all the spheres
NUM_FILES = 11
folder1 = '.\\DATA\\'

for kk in range(0,NUM_FILES):
    fname1 = folder1+'SPH'+ str(kk+101)+'.xyz' #concatenate the file name
    data1 = np.genfromtxt(fname1, delimiter='\t') #Read the data
    
    # Perform the segmentation per the ASTM E57-3125 algorithm
    [dataFinal, dataIgnored, finalCenter] = segfun.coneCylAlgoE57(data1,trueRadius,120) 
    
    #Extract the center and calculate the sum of squares
    [x0, y0, z0] = finalCenter; 
    [x, y, z] = dataFinal.T
    resids1 = np.sqrt((x-x0)**2 + (y-y0)**2 + (z-z0)**2) - trueRadius
    ssq1 = np.sum(resids1**2)
    len1 = len(data1)
    lenF = len(dataFinal)

    #Display results
    print 'Center = %2.9f, %2.9f, %2.9f, InitPoints = %4d, FinalPoints = %4d' % (x0,y0,z0,len1,lenF)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(dataFinal[:,0], dataFinal[:,1], dataFinal[:,2], c='r', marker='o')
    ax.scatter(dataIgnored[:,0], dataIgnored[:,1], dataIgnored[:,2], c='b', marker='o')
    plt.show()        
    #time.sleep(0.1)
    
    