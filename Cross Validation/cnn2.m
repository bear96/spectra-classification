clc, clear all, close all;
XTrain = xlsread('MSC+1st DERV.xlsx');
XXTrain = ResampledData(XTrain);
height = 1000;
width = 1;
channels = 1;
samples = 105;
CNN_TrainingData = reshape(XXTrain,[height,width,channels, samples]);
YTrain = xlsread('trng ref data_caffeine-35x3_105 sample.xls');
CNN_TrainingLabels = YTrain;
auimds = augmentedImageDatastore([1000 1],CNN_TrainingData,CNN_TrainingLabels,'ColorPreprocessing', 'none');
miniBatchSize  = 3;
validationFrequency = floor(numel(YTrain)/miniBatchSize);
layers = [
    imageInputLayer([height,width, channels])

    convolution2dLayer([5 1],100, 'stride',1)
    batchNormalizationLayer
    reluLayer
    
    averagePooling2dLayer([50 1],'Stride',2)
    
    dropoutLayer(0.3)
    fullyConnectedLayer(1)
    regressionLayer];

options = trainingOptions('sgdm', ...
    'MiniBatchSize',miniBatchSize, ...
    'MaxEpochs',50, ...
    'InitialLearnRate',1e-3, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.1, ...
    'LearnRateDropPeriod',30);
            
net = trainNetwork(auimds,layers,options);
% view(net)
YPrediction = predict(net,CNN_TrainingData);
[~,rmsep,~,~]= performance(YTrain,YPrediction);
R1=  corr(YTrain,YPrediction,'Type','Pearson');
