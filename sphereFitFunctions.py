# Prem Rachakonda (2017)
#
# This software was developed by employees of the National Institute of Standards and Technology (NIST), an agency of the Federal Government. Pursuant to title 17 United States Code Section 105, works of NIST employees are not subject to copyright protection in the United States and are considered to be in the public domain. Permission to freely use, copy, modify, and distribute this software and its documentation without fee is hereby granted, provided that this notice and disclaimer of warranty appears in all copies.
# THE SOFTWARE IS PROVIDED 'AS IS' WITHOUT ANY WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND FREEDOM FROM INFRINGEMENT, AND ANY WARRANTY THAT THE DOCUMENTATION WILL CONFORM TO THE SOFTWARE, OR ANY WARRANTY THAT THE SOFTWARE WILL BE ERROR FREE. IN NO EVENT SHALL NIST BE LIABLE FOR ANY DAMAGES, INCLUDING, BUT NOT LIMITED TO, DIRECT, INDIRECT, SPECIAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF, RESULTING FROM, OR IN ANY WAY CONNECTED WITH THIS SOFTWARE, WHETHER OR NOT BASED UPON WARRANTY, CONTRACT, TORT, OR OTHERWISE, WHETHER OR NOT INJURY WAS SUSTAINED BY PERSONS OR PROPERTY OR OTHERWISE, AND WHETHER OR NOT LOSS WAS SUSTAINED FROM, OR AROSE OUT OF THE RESULTS OF, OR USE OF, THE SOFTWARE OR SERVICES PROVIDED HEREUNDER.
# Distributions of NIST software should also include copyright and licensing statements of any third-party software that are legally bundled with the code in compliance with the conditions of those licenses.

import numpy as np
#from scipy.optimize import leastsq
from scipy.optimize import least_squares #leastsq is about 8% faster than least_squares in this trial
    
def sphereFitLin1(data1):
    [xx,yy,zz] = data1.T
    oo = xx*0+1; #simpler than using ones()
    
    AA = [-2*xx, -2*yy , -2*zz, oo];
    BB = [-xx**2-yy**2-zz**2]
    YY = np.linalg.lstsq(np.transpose(AA),np.transpose(BB))
    
    [a,b,c,D] = YY[0][0:4]
    r = np.sqrt((a**2+b**2+c**2)-D);

    result = np.hstack([a,b,c,r])
    return result

def sphereFitLin2(data1, knownRadius):
    [xx,yy,zz] = data1.T
    r = knownRadius;
    oo = xx*0+1; #simpler than using ones()
    
    AA = [-2*xx, -2*yy , -2*zz, oo];
    BB = [-xx**2-yy**2-zz**2 + oo*(r**2)];
    YY = np.linalg.lstsq(np.transpose(AA),np.transpose(BB))
    
    [a,b,c,D] = YY[0][0:4]
    

    result = np.hstack([a,b,c,r])
    return result

def sphereFitLSQ1_uncR(data1,dummyVal=-1):   
    initGuess = sphereFitLin1(data1)

    def calcResids1(initGuess1, data1):
        [x0, y0, z0, R] = initGuess1
        [x, y, z] = data1.T
        resids1 = np.sqrt((x-x0)**2 + (y-y0)**2 + (z-z0)**2) - R
        return resids1
    
    
    #result, flag = leastsq(calcResids1, initGuess, args=(data1,),maxfev=1500,ftol=1E-15, xtol=1E-15) #Same arguments as MATLAB's lsqnonlin() for consistent results
    result = least_squares(calcResids1, initGuess, args=(data1,), method = 'lm', max_nfev=1500,ftol=1E-15, xtol=1E-15) #Same arguments as MATLAB's lsqnonlin() for consistent results
    result = result.x
    return result


def sphereFitLSQ1_conR(data1,knownRadius, initGuess = None):
    if initGuess is None:
        initGuess = sphereFitLin2(data1, knownRadius);
    initGuess = initGuess[0:3]

    def calcResids2(initGuess2, data1,knownRadius):
        [x0, y0, z0] = initGuess2
        R = knownRadius;
        [x, y, z] = data1.T
        resids1 = np.sqrt((x-x0)**2 + (y-y0)**2 + (z-z0)**2) - R
        return resids1

    
    #result, flag = leastsq(calcResids2, initGuess, args=(data1,knownRadius),maxfev=1500,ftol=1E-15, xtol=1E-15) #Same arguments as MATLAB's lsqnonlin() for consistent results
    result = least_squares(calcResids2, initGuess, args=(data1,knownRadius),method='lm',max_nfev=1500,ftol=1E-15, xtol=1E-15) #Same arguments as MATLAB's lsqnonlin() for consistent results
    result = result.x;

    result = np.append(result,knownRadius) # Sending back a 4 element array
    return result
    
def cart2sph(x,y,z): #Converts the coordinates from Cartesian to spherical coordinate system
    #import numpy as np
    hypotxy = np.hypot(x,y);
    rr = np.hypot(hypotxy,z);
    el = np.arctan2(z,hypotxy);
    az = np.arctan2(y,x);
    return az, el, rr

def sph2cart(az,el,rr): #Converts coordinates from the spherical to Cartesian coordinate system
    #import numpy as np
    z = rr * np.sin(el);
    rcosel = rr * np.cos(el);
    x = rcosel * np.cos(az);
    y = rcosel * np.sin(az);
    return x,y,z