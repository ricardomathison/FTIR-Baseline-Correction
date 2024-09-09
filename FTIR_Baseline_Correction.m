function linefit = FTIR_Baseline_Correction(rawdata,xwindow,xpeak,xlocmin)
%
%% Description
% 
% linefit: Corrected FTIR data with a flat baseline on both sides of the
% peak of interest
% rawdata: 2-column matrix (x values = column 1, y values = column 2)
% xwindow: matrix of two values with the initial and final values of the
% x-axis window that will be corrected
% xpeak: matrix of two values with the initial and final values of the 
% x-axis where the peak(s) of interest is(are) located
% xlocmin: matrix of two values with the initial and final values of the 
% x-axis where the local minimum of the peak(s) is(are) located. If there
% is no local minimum, write []
% 
% Algorithm steps:
% 1. Predicts a baseline by fitting third-degree polynomials to the
% baseline, one on each side of the peak(s) of interest
% 2. Calculates a sigmoidal average between the two poynomials, using an 
% S-shaped weighting function to smoothly transition from one polynomial 
% to the other across the data range
% 3. Corrects the raw data by substracting the predicted baseline
% 
%% Author Info:
% Ricardo Mathison, September 9, 2024.
%
%% sigmoid(x,c,a) function Author Info: Chad Greene, May 28, 2015. http://www.chadagreene.com
%
%% Algorithm:

xlim1 = rawdata(1,1); xlim2 = rawdata(end,1);

smooth = 1;
dataStruct = rawdata;
x = dataStruct(:,1);
[xminlimit,xminindex] = min(abs(x - xlim1));
[xmaxlimit,xmaxindex] = min(abs(x - xlim2));
x = x(xminindex:xmaxindex);
x = movmean(x,smooth);
xi = x(1):0.01:x(end);
y = dataStruct(:,2);
y = y(xminindex:xmaxindex);
yi = spline(x,y,xi);
yii = movmean(yi,smooth);

[a,xwindex1] = min(abs(xi - xwindow(1)));
[a,xwindex2] = min(abs(xi - xwindow(2)));
[a,xpeak1] = min(abs(xi - xpeak(1)));
[a,xpeak2] = min(abs(xi - xpeak(2)));
if isempty(xlocmin)
    xlocmin1 = 1;
    xlocmin2 = 0;
else
    [a,xlocmin1] = min(abs(xi - xlocmin(1)));
    [a,xlocmin2] = min(abs(xi - xlocmin(2)));
end

xi1 = xi(xwindex1:xpeak1);
yii1 = yii(xwindex1:xpeak1);
xi2 = xi(xpeak2:xwindex2);
yii2 = yii(xpeak2:xwindex2);
xipeak = xi(xpeak1:xpeak2);
yiipeak = yii(xpeak1:xpeak2);
xilocmin = xi(xlocmin1:xlocmin2);
yiilocmin = yii(xlocmin1:xlocmin2);

xii = [xi1 xi2];
yiii = [yii1 yii2];

xi1len = length(xi1);
xi2len = length(xi2);

N = 50;
xi1simple = zeros (1,N);
yii1simple = zeros (1,N);
xi2simple = zeros (1,N);
yii2simple = zeros (1,N);
ki1 = floor(xi1len/N);
ki2 = floor(xi2len/N);

for i = 1:N
    xi1simple(i) = xi1(i*ki1);
    yii1simple(i) = yii1(i*ki1);
    xi2simple(i) = xi2(i*ki2);
    yii2simple(i) = yii2(i*ki2);
end

TF = islocalmin(yiilocmin);

xi1simple = [xi1simple xilocmin(TF)];
yii1simple = [yii1simple yiilocmin(TF)];
xi2simple = [xilocmin(TF) xi2simple];
yii2simple = [yiilocmin(TF) yii2simple];

if isempty(xlocmin)
    N2 = 5;
    xpeakline = linspace(xi1simple(end),xi2simple(1),N2);
    ypeakline = linspace(yii1simple(end),yii2simple(1),N2);
    xsimpleall = [xi1simple xpeakline xi2simple];
    ysimpleall = [yii1simple ypeakline yii2simple];

    degree = 3;
    pall = polyfit(xsimpleall,ysimpleall,degree);
    ypeakline2 = polyval(pall,xpeakline);

    xi1simple = [xi1simple xpeakline];
    yii1simple = [yii1simple ypeakline2];
    xi2simple = [xpeakline xi2simple];
    yii2simple = [ypeakline2 yii2simple];
end

degree = 3;
p1 = polyfit(xi1simple,yii1simple,degree);
xi1peak = xi(xwindex1:xpeak2);
ypredict1 = polyval(p1,xi1);

p2 = polyfit(xi2simple,yii2simple,degree);
xi2peak = xi(xpeak1:xwindex2);
ypredict2 = polyval(p2,xi2);

xgradient = sigmoid(linspace(-10,10,length(xipeak)));
ypredictpeak1 = polyval(p1,xipeak);
ypredictpeak2 = polyval(p2,xipeak);

ypredictwave = [ypredict1 ((1 - xgradient).*ypredictpeak1 + xgradient.*ypredictpeak2) ypredict2];

linefit = [[xi1 xipeak xi2]; [yii1 yiipeak yii2] - ypredictwave];

end