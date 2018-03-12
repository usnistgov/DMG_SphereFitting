function result = sphereFitLSQ1_conR(data,knownRadius,initGuess)
% This function takes data of a sphere in a Cartesian coordinate system 
% and returns the paramters of a sphere (Center and radius)
%
%The input data is in the form of Nx3 matrix
%The output is a structure with the sphere parameters.
% Prem Rachakonda (2017)
%
% This software was developed by employees of the National Institute of Standards and Technology (NIST), an agency of the Federal Government. Pursuant to title 17 United States Code Section 105, works of NIST employees are not subject to copyright protection in the United States and are considered to be in the public domain. Permission to freely use, copy, modify, and distribute this software and its documentation without fee is hereby granted, provided that this notice and disclaimer of warranty appears in all copies.
% THE SOFTWARE IS PROVIDED 'AS IS' WITHOUT ANY WARRANTY OF ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTY THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, AND FREEDOM FROM INFRINGEMENT, AND ANY WARRANTY THAT THE DOCUMENTATION WILL CONFORM TO THE SOFTWARE, OR ANY WARRANTY THAT THE SOFTWARE WILL BE ERROR FREE. IN NO EVENT SHALL NIST BE LIABLE FOR ANY DAMAGES, INCLUDING, BUT NOT LIMITED TO, DIRECT, INDIRECT, SPECIAL OR CONSEQUENTIAL DAMAGES, ARISING OUT OF, RESULTING FROM, OR IN ANY WAY CONNECTED WITH THIS SOFTWARE, WHETHER OR NOT BASED UPON WARRANTY, CONTRACT, TORT, OR OTHERWISE, WHETHER OR NOT INJURY WAS SUSTAINED BY PERSONS OR PROPERTY OR OTHERWISE, AND WHETHER OR NOT LOSS WAS SUSTAINED FROM, OR AROSE OUT OF THE RESULTS OF, OR USE OF, THE SOFTWARE OR SERVICES PROVIDED HEREUNDER.
% Distributions of NIST software should also include copyright and licensing statements of any third-party software that are legally bundled with the code in compliance with the conditions of those licenses.

%% Check the dimensions of the data
[m,n] = size(data);

if(n<3 | m<4)
    data = data';
end
[m,n] = size(data);

if(n<3 | m<4)
    error(sprintf('Not enough parameters (size of matrix): SIZE = %d %d',m,n));
end

if nargin<2 
    error('Not enough parameters');
elseif isempty(knownRadius) | knownRadius == 0
    error('Radius = 0');
end

%% First get an initial estimate based on a linearized form of the sphere equation
if nargin<3
    [ro,ao,bo,co,residualo] = sphereFitInitGuess2(data, knownRadius); %this function returns the radius. We don't need to do a constrained fit for an initial estimate.
    initGuess = [ao bo co];
end
    

%% Now formulate the objective function and the options for lsqnonlin() 
if (exist('optimoptions')) %Some versions of MATLAB have optimset
    options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','MaxFunEvals',1500,'TolFun',1E-15,'TolX',1E-15,'display','off');
else
    options = optimset('Algorithm','levenberg-marquardt','MaxFunEvals',1500,'TolFun',1E-15,'TolX',1E-15,'display','off');
end
ObjectiveFunction = @(initGuess)calcResiduals(initGuess,knownRadius,data);
[result3] = lsqnonlin(ObjectiveFunction,initGuess,[],[],options);


   
%Since we are constraining the sphere, just to keep the parameters
%consistent, assign the 4th element to 'knownRadius'
result3(4) = knownRadius;


result.Parameters = result3;
result.Residuals  = calcResiduals(result3,knownRadius,data); %For legacy
result.Center     = result3(1:3)';
result.Radius     = knownRadius;



function residuals = calcResiduals(initGuess,knownRadius,data)
%Function to calculate the residuals based on the sphere parameters and the
%data.

ao = initGuess(1);
bo = initGuess(2);
co = initGuess(3);


xx = data(:,1);
yy = data(:,2);
zz = data(:,3);

residuals = sqrt( (xx-ao).^2 + (yy-bo).^2 + (zz-co).^2 ) - knownRadius;

function [r,a,b,c,residual] = sphereFitInitGuess2(data, knownRadius)
%This is a function to calculate the parameters of a sphere as an initial
%estimate. This calculates the center of the sphere for a given radius
xx = data(:,1);
yy = data(:,2);
zz = data(:,3);
r = knownRadius;

AA = [-2*xx, -2*yy , -2*zz , ones(size(xx))];
BB = [ -(xx.^2+yy.^2+zz.^2)] + r^2;

YY = mldivide(AA,BB);

a = YY(1);
b = YY(2);
c = YY(3);
residual = AA*YY-BB;
