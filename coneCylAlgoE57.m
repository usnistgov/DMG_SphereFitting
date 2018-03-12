function [dataFinal, dataIgnored, centerFinal2, resultsFinal] = coneCylAlgoE57(data1,trueRadius)
% function [dataFinal, dataIgnored, resultsFinal] = coneCylAlgoE57(data1,trueRadius,coneAngle)
%
%This function takes in a dataset from a TLS that has a scan of a sphere
%and obtains the data that is processed according to the ASTM E57.02
%method. This file has 3 main functions:
% 1. coneCylAlgoE57(): This is the main function that is called and this
%    functions performs the iterations using the next 2 functions.
% 2. loc_closestPointAlgo(): This is a function that obtains an approximate 
%    sphere center
% 3. loc_coneCylAlgo(): This is the actual function which truncates the 
%    data outside a conical+cylindrical region

%If this is being used under Octave, we need to first install the
% optimization package and then load it.
% Use 'pkg install -forge optim' at Octave prompt to install the package
%
% Prem Rachakonda & Bala Muralikrishnan(2017)
%
% This software was developed by employees of the National Institute of Standards and Technology (NIST), an agency of the Federal Government. Pursuant to title 17 United States Code Section 105, works of NIST employees are not subject to copyright protection in the United States and are considered to be in the public domain. Permission to freely use, copy, modify, and distribute this software and its documentation without fee is hereby granted, provided that this notice and disclaimer of warranty appears in all copies.
% THE SOFTWARE IS PROVIDED 'AS IS' WITHOUT ANY WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND FREEDOM FROM INFRINGEMENT, AND ANY WARRANTY THAT THE DOCUMENTATION WILL CONFORM TO THE SOFTWARE, OR ANY WARRANTY THAT THE SOFTWARE WILL BE ERROR FREE. IN NO EVENT SHALL NIST BE LIABLE FOR ANY DAMAGES, INCLUDING, BUT NOT LIMITED TO, DIRECT, INDIRECT, SPECIAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF, RESULTING FROM, OR IN ANY WAY CONNECTED WITH THIS SOFTWARE, WHETHER OR NOT BASED UPON WARRANTY, CONTRACT, TORT, OR OTHERWISE, WHETHER OR NOT INJURY WAS SUSTAINED BY PERSONS OR PROPERTY OR OTHERWISE, AND WHETHER OR NOT LOSS WAS SUSTAINED FROM, OR AROSE OUT OF THE RESULTS OF, OR USE OF, THE SOFTWARE OR SERVICES PROVIDED HEREUNDER.
% Distributions of NIST software should also include copyright and licensing statements of any third-party software that are legally bundled with the code in compliance with the conditions of those licenses.

sw = ver;
if (strcmp(sw(1).Name,'Octave') == 1)
    pkg load optim
end


ITER = 6; %We specify only 5 iterations. This can be increased.
coneAngle = 120;
halfConeAngle = coneAngle/2;

for kk = 1:ITER
    if kk == 1 %Perform the first iteration and get the initial center
        [dataFinal,dataIgnored,centerInit,centerFinal(kk,:)] = loc_closestPointAlgo(data1,trueRadius, halfConeAngle);
        %The cenerInit in iter#1 is the initial center for the next iter#         
    else 
        newCenter = centerFinal(kk-1,:);
        [dataFinal,dataIgnored,~,centerFinal(kk,:)] = loc_coneCylAlgo(data1,trueRadius, halfConeAngle,	newCenter);
    end
end
centerFinal2 = centerFinal(kk,:);
resultsFinal = sphereFitLSQ1_conR(dataFinal,trueRadius,centerFinal2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dataFinal,dataIgnored,centerInit,centerFinal] = loc_closestPointAlgo(data1,trueRadius, halfConeAngle)
%% This is the first pass where an approximate sphere center is calculated.
% This is called the "closest point method"
% Schematic: https://raw.githubusercontent.com/usnistgov/DMG_SphereFitting/master/Schematic_ClosestPointMethod.png
POINTS  = 500;
PERCENT = 0.05;
% First sort the data from the closest to furthest
[~,~,rng1] = cart2sph(data1(:,1),data1(:,2),data1(:,3));
[rng2,~]   = sort(rng1);

% Then take the closest N points at the start of the surface and grab
% data that is within the trueRadius distance
len1    = floor(min(POINTS, PERCENT*length(data1)));
if (len1 < 4) %For low density data, the 5% of points may be lower than 4
    len1 = max(4,ceil(0.1*length(data1)) );
end %Just a check

% Calculate the median distance of the points closest to the scanner
surfaceStart = median(rng2(1:len1));

% Extract points that are from the closest point to 0.5*radius
% and find the center
idx2        = rng1<(surfaceStart+0.5*trueRadius);
dataFinal   = data1(idx2,:);
res1        = sphereFitLSQ1_conR(dataFinal,trueRadius);
centerInit  = res1.Center';

%Separate out the data discarded; Return the center (in this case, the
%final center is also the initial center). Just kept it for consistency
%with the next function. 
idx9        = ismember(data1,dataFinal,'rows');
dataIgnored = data1(~idx9,:);
centerFinal = centerInit;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dataFinal,dataIgnored,centerInit,centerFinal] = loc_coneCylAlgo(data1,trueRadius, halfConeAngle,centerInit)
%This is a cone-cylinder truncation algorithm, where data from the scan of
%a sphere from a laser scanner is processed to obtain a sphere center.
%
%Schematic: https://raw.githubusercontent.com/usnistgov/DMG_SphereFitting/master/Schematic_DistanceOfPointsFromAxis_CylinderRegion.png
%Calculate the radius of the cylinder corresponding to the cone angle at
%the sphere center
cylRadius = trueRadius*sin(halfConeAngle*pi/180); %Cylinder radius (h)


%% Now to:  finding data within a cone angle of 120 degrees
%Let's say O is the origin, C is the center, and P is a point on the sphere
%The angle between the vectors CO and CP should be within 60 degrees
%here halfConeAngle = 60 degrees
origin      = zeros(size(data1));
if size(centerInit,1) ==3 %Just a check for the shape of the matrix
    centerInit = centerInit';
end
vectorCO    = bsxfun(@minus,origin, centerInit);
vectorCP    = bsxfun(@minus,data1 , centerInit);
allAngles1  = vectorAngle3(vectorCO, vectorCP)*180/pi; %in degrees

%Find the points in the original data set that are within the desired
%cone angle for the segmentation
idx4        = allAngles1 < halfConeAngle; %Points within the cone
vectorOP    = bsxfun(@minus,data1 , origin);

%The code below finds the distance of all the points from the axis joining
%the initial center and the origin - to find points within a cylinder
%O is the origin, C is the sphere center, and P is a point on the sphere
%F is a point on the line CO such that PF is perpendicular to CO
%Length of PF is the shortest distance between the point P and CO.
%dOF = |OF| is the projection of OP on the unit vector along CO.
unitVectorCO = vectorCO/rssq3(vectorCO(1,:),2); %All rows of vectorCO are identical
dOF   = dot(vectorOP', unitVectorCO')';
dOP   = rssq3(vectorOP,2);
dPF   = sqrt(dOP.^2 - dOF.^2); %Distance from each point to the line CO
idx5  = dPF < cylRadius; %Points within the cylinder

%Now combine the points that are both within the cone and the cylinder
idxA  = idx4 & idx5;
data2 = data1(idxA,:);
%dataIgnored = data1(~idxA,:);

%Truncate points whose residuals exceed 3sigma
res1    = sphereFitLSQ1_conR(data2,trueRadius,centerInit);
resids1 = res1.Residuals;
idx8    = abs(resids1) < 3*std(resids1);
center1 = res1.Center';

dataFinal   = data2(idx8,:);
idx8a       = ismember(data1,dataFinal,'rows');
dataIgnored = data1(~idx8a,:);

%Find the center of the final dataset
res2        = sphereFitLSQ1_conR(dataFinal,trueRadius,center1);
centerFinal = res2.Center;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vAngle = vectorAngle3(vec1, vec2)
%This function calculates the angle between two vectors. vec1 and vec2 can
%by Nx3 matrices for N vectors
sinValue = rssq3(cross(vec1,vec2),2);
cosValue = dot(vec1',vec2')';
vAngle   = atan2(sinValue,cosValue);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vals = rssq3(data1,dim)
%This function calculates the root-sum-square of the input data along a
%dimension 'dim'
sq1  = data1.*data1;
vals = sqrt(sum(sq1,dim));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vals = rms3(data1, dim)
%This function calculates the root-mean-square of the input data along a
%dimension 'dim'
sq1  = data1.*data1;
vals = sqrt(mean(sq1,dim));