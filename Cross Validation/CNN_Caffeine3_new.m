clc, clear all, close all;
XTrain = xlsread('MSC+1st DERV.xlsx');
% XTrain_Norm = normalize(XTrain, 'range');
XXTR = XTrain';
XXTrain = ResampledData(XTrain);
height = 7000;
width = 1;
channels = 1;
samples = 105;
CNN_TrainingData = reshape(XXTrain,[height,width,channels, samples]);
YTrain = xlsread('trng ref data_caffeine-35x3_105 sample.xls');
CNN_TrainingLabels = YTrain;

layers = [
    imageInputLayer([height,width, channels])

    convolution2dLayer([5 1],100, 'stride',1)
    batchNormalizationLayer
    reluLayer
    
    averagePooling2dLayer([50 1],'Stride',2)

%     convolution2dLayer(3,16, 'Padding', 'same')
%     batchNormalizationLayer
%     reluLayer
%     
%     averagePooling2dLayer(1,'Stride',1)
  
%     convolution2dLayer(3,32, 'Padding','same')
%     batchNormalizationLayer
%     reluLayer
%     
%     convolution2dLayer(3,32, 'Padding','same')
%     batchNormalizationLayer
%     reluLayer
    
    dropoutLayer(0.3)
    fullyConnectedLayer(1)
    regressionLayer];

miniBatchSize  = 2;
validationFrequency = floor(numel(YTrain)/miniBatchSize);
options = trainingOptions('sgdm', ...
    'MiniBatchSize',miniBatchSize, ...
    'MaxEpochs',30, ...
    'InitialLearnRate',1e-3, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.1, ...
    'LearnRateDropPeriod',100);
            
net = trainNetwork(CNN_TrainingData,CNN_TrainingLabels,layers,options);
YPrediction = predict(net,CNN_TrainingData);
[~,rmsep,~,~]= performance(YTrain,YPrediction);
R1=  corr(YTrain,YPrediction,'Type','Pearson');