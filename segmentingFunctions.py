# Prem Rachakonda (2017)
#
# This software was developed by employees of the National Institute of Standards and Technology (NIST), an agency of the Federal Government. Pursuant to title 17 United States Code Section 105, works of NIST employees are not subject to copyright protection in the United States and are considered to be in the public domain. Permission to freely use, copy, modify, and distribute this software and its documentation without fee is hereby granted, provided that this notice and disclaimer of warranty appears in all copies.
# THE SOFTWARE IS PROVIDED 'AS IS' WITHOUT ANY WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND FREEDOM FROM INFRINGEMENT, AND ANY WARRANTY THAT THE DOCUMENTATION WILL CONFORM TO THE SOFTWARE, OR ANY WARRANTY THAT THE SOFTWARE WILL BE ERROR FREE. IN NO EVENT SHALL NIST BE LIABLE FOR ANY DAMAGES, INCLUDING, BUT NOT LIMITED TO, DIRECT, INDIRECT, SPECIAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF, RESULTING FROM, OR IN ANY WAY CONNECTED WITH THIS SOFTWARE, WHETHER OR NOT BASED UPON WARRANTY, CONTRACT, TORT, OR OTHERWISE, WHETHER OR NOT INJURY WAS SUSTAINED BY PERSONS OR PROPERTY OR OTHERWISE, AND WHETHER OR NOT LOSS WAS SUSTAINED FROM, OR AROSE OUT OF THE RESULTS OF, OR USE OF, THE SOFTWARE OR SERVICES PROVIDED HEREUNDER.
# Distributions of NIST software should also include copyright and licensing statements of any third-party software that are legally bundled with the code in compliance with the conditions of those licenses.

import numpy as np
import sphereFitFunctions as sfun

def coneCylAlgoE57(data1,trueRadius,coneAngle=120):
    ITER = 6;
    halfConeAngle = coneAngle/2;
    centerFinal = np.zeros((ITER,3))
    for kk in  range(0,ITER):
        #print centerFinal
        if kk == 0:
            [dataFinal,dataIgnored,centerInit,_] = loc_closestPointMethod(data1,trueRadius)            
            centerFinal[kk,:] = centerInit
        else:
            newCenter = centerFinal[kk-1,:];
            [dataFinal,dataIgnored,centerInit,cfTemp] = loc_coneCylAlgo(data1,trueRadius,halfConeAngle,newCenter);
            centerFinal[kk,:] = cfTemp
            
    finalCenter = centerFinal[kk,:];
    results = [dataFinal, dataIgnored, finalCenter];
    return results;

def loc_closestPointMethod(data1,trueRadius):
    import sphereFitFunctions as sfun    
    
    POINTS  = 500;
    PERCENT = 0.05;   

    # First sort the ranging data from the closest to furthest
    rng1 = rssq3(data1,1)
    rng2 = np.sort(rng1);
    
    
    # Then take the closest N points at the start of the surface and grab
    # data that is within the trueRadius distance
    len1 = int(np.floor(min(POINTS, PERCENT*len(data1))))    
    if (len1 < 4): #Just a check, and to ensure a min. of 4 points
        len1 = int(max(4, np.ceil(0.1*len(data1))))
        
    surfaceStart = np.median(rng2[0:len1-1]); #Select the first len1 points
    
    # Extract points that are from the closest point to 0.5*radius
    # and find the center
    idx2          = rng1<(surfaceStart+0.5*trueRadius);
    dataFinal     = data1[idx2,:];
    [cx,cy,cz,rr] = sfun.sphereFitLSQ1_conR(dataFinal,trueRadius);
    centerInit    = np.asarray([cx,cy,cz]);
    
    dataIgnored   = data1[~idx2,:];
    centerFinal   = centerInit;
    results       = [dataFinal,dataIgnored,centerInit,centerFinal];
    return results;

def loc_coneCylAlgo(data1,trueRadius, halfConeAngle,centerInit):
#This is a cone-cylinder truncation algorithm, where data from the scan of
#a sphere from a laser scanner is processed to obtain a sphere center.
#
#Schematic: https://raw.githubusercontent.com/usnistgov/DMG_SphereFitting/master/Schematic_DistanceOfPointsFromAxis_CylinderRegion.png
#Calculate the radius of the cylinder corresponding to the cone angle at
#the sphere center
    cylRadius = trueRadius*np.sin(halfConeAngle*np.pi/180); #Cylinder radius (h)
    
    
    ## Now to:  finding data within a cone angle of 120 degrees
    #Let's say O is the origin, C is the center, and P is a point on the sphere
    #The angle between the vectors CO and CP should be within 60 degrees
    #here halfConeAngle = 60 degrees
    origin      = data1*0; #Easier way to get zeros instead of zeros()
    if centerInit.shape[0] == 3:
        centerInit = centerInit.T;
    vectorCO    = origin - centerInit;
    vectorCP    = data1 - centerInit;
    allAngles1  = vectorAngle3(vectorCO, vectorCP)*180/np.pi; #in degrees
    
    #Find the points in the original data set that are within the desired
    #cone angle for the segmentation
    idx4        = allAngles1 < halfConeAngle; #Points within the cone
    vectorOP    = data1 - origin;
    
    #The code below finds the distance of all the points from the axis joining
    #the initial center and the origin - to find points within a cylinder
    #O is the origin, C is the sphere center, and P is a point on the sphere
    #F is a point on the line CO such that PF is perpendicular to CO
    #Length of PF is the shortest distance between the point P and CO.
    #dOF = |OF| is the projection of OP on the unit vector along CO.
    dCO = rssq3(vectorCO[1,:],0)
    unitVectorCO = vectorCO/dCO; #All rows of vectorCO are identical
    #dOF   = np.dot(vectorOP.T, unitVectorCO.T);dOF = dOF.T;
    dOF   = vecDot(vectorOP,unitVectorCO);
    dOP   = rssq3(vectorOP,1);
    dPF   = np.sqrt(dOP*dOP - dOF*dOF); #Distance from each point to the line CO
    idx5  = dPF < cylRadius; #Points within the cylinder
    
    #Now combine the points that are both within the cone and the cylinder
    idxA = idx4&idx5;
    data2   = data1[idxA,:];
    data2I = data1[~idxA,:];
    
    #Truncate points whose residuals exceed 3sigma
    [cx0,cy0,cz0,rr0]    = sfun.sphereFitLSQ1_conR(data2,trueRadius,centerInit);
    center2 = [cx0,cy0,cz0];
    resids2  = rssq3(data2-center2,1)-rr0;
    idx8     = np.abs(resids2)<3*np.std(resids2,ddof=1);
    
    dataFinal    = data2[idx8,:];
    data1I = data2[~idx8,:];
    dataIgnored  = np.concatenate([data2I,data1I]);
    
    
    #Find the center of the final dataset
    [cxf,cyf,czf,_]= sfun.sphereFitLSQ1_conR(dataFinal,trueRadius,center2);
    centerFinal = [cxf,cyf,czf]
    
    results = [dataFinal, dataIgnored, centerInit, centerFinal]
    return results;


def vectorAngle3(vec1,vec2): #calculates the angle between two vectors. Each vector can be a Nx3 matrix
    #import numpy as np
    cp = np.cross(vec1,vec2);
    dp = sum([i*j for (i, j) in zip(vec1.T, vec2.T)])    
#    dp = []; #In MATLAB, this is a non-loop operation (i.e. matrix operation, here we have to do it in a loop)
#    for jj in xrange(len(vec1)):
#        dp = np.append(dp, np.dot(vec1[jj,:], vec2[jj,:]))
    sinValue = np.sqrt(np.sum(cp*cp,1))
    cosValue = dp
    vAngle = np.arctan2(sinValue,cosValue)
    return vAngle

def rssq3(data1,dim):    #Calculates the root-sum-square of a matrix
    sq1 = data1*data1
    vals = np.sqrt(np.sum(sq1,dim))
    return vals

def rms3(data1,dim): #Calculates the root-mean-square value of a matrix
    #import numpy as np
    sq1 = data1*data1
    vals = np.sqrt(np.mean(sq1,dim))
    return vals

def vecDot(vec1,vec2): #Calculates dot product of two Nx3 matrices
    dp = sum([i*j for (i, j) in zip(vec1.T, vec2.T)])    
#   dp = [];
#    for jj in xrange(len(vec1)):
#        dp = np.append(dp, np.dot(vec1[jj,:], vec2[jj,:]))
    return dp