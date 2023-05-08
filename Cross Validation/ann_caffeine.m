clc, clear all, close all;
%   XXTR - input data.
%   data_1 - target data.
XXTR = xlsread('MSC+1st DERV.xlsx');
data_1 = xlsread('trng ref data_caffeine-35x3_105 sample.xls');
x = XXTR;
t = data_1';
trainFcn = 'trainlm';  % Levenberg-Marquardt backpropagation.

hiddenLayerSize = 10;
net = fitnet(hiddenLayerSize,trainFcn);

net.divideParam.trainRatio = 80/100;
net.divideParam.valRatio = 10/100;
net.divideParam.testRatio = 10/100;

% Train the Network
[net,tr] = train(net,x,t);

% Test the Network
x_test = xlsread('MSC_1ST_DERV_test.xls');
y = net(x_test);
test_ref = xlsread('test ref data_caffeine_7x3_15 samples.xls');
test_ref = test_ref';
e = gsubtract(test_ref,y);
performance = perform(net,test_ref,y);
R = corr(test_ref',y');
% View the Network
view(net)

% Plots
% Uncomment these lines to enable various plots.
%figure, plotperform(tr)
%figure, plottrainstate(tr)
%figure, ploterrhist(e)
%figure, plotregression(t,y)
%figure, plotfit(net,x,t)

