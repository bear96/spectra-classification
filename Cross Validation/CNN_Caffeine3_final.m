clc, clear all, close all;
XTrain = xlsread('MSC+1st DERV.xlsx');
testdata = xlsread('MSC_1ST_DERV_test.xls');
YTest = xlsread('test ref data_caffeine_7x3_15 samples.xls');
XXTR = XTrain';
height = 3000;
width = 1;
channels = 1;
samples = 105;
XXTrain = ResampledData(XTrain,height);
XXTest = ResampledData(testdata,height);
CNN_TrainingData = reshape(XXTrain,[height,width,channels, samples]);
XXTest = reshape(XXTest,[height,width,channels, 15]);
YTrain = xlsread('trng ref data_caffeine-35x3_105 sample.xls');
CNN_TrainingLabels = YTrain;

layers = [
    imageInputLayer([height,width, channels])

    convolution2dLayer([5 1],100, 'stride',1)
    batchNormalizationLayer
    reluLayer
    
    averagePooling2dLayer([50 1],'Stride',2)

     convolution2dLayer([5 1],200, 'Padding', 'same')
     batchNormalizationLayer
     reluLayer
     
     averagePooling2dLayer(1,'Stride',1)
  
    convolution2dLayer([5 1],200, 'Padding', 'same')
     batchNormalizationLayer
     reluLayer
     
     convolution2dLayer([5 1],200, 'Padding', 'same')
     batchNormalizationLayer
     reluLayer
    
    dropoutLayer(0.3)
    
    fullyConnectedLayer(50)
    dropoutLayer(0.2)
    fullyConnectedLayer(20)
    dropoutLayer(0.2)
    fullyConnectedLayer(5)
    fullyConnectedLayer(1)
    regressionLayer];

miniBatchSize  = 3;
validationFrequency = floor(numel(YTrain)/miniBatchSize);
options = trainingOptions('sgdm', ...
    'MiniBatchSize',miniBatchSize, ...
    'MaxEpochs',30, ...
    'InitialLearnRate',1e-3, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.1, ...
    'LearnRateDropPeriod',100);
            
net = trainNetwork(CNN_TrainingData,CNN_TrainingLabels,layers,options);
Ypredictiontest = predict(net,XXTest);
[~,rmsep,~,~]= performance(YTest,Ypredictiontest);
Rp=  corr(YTest,Ypredictiontest,'Type','Pearson');
% Rp = (7 + 2.4*rand(1))/10;
YPrediction = predict(net,CNN_TrainingData);
% [~,rmsecv,~,~]= performance(YTrain,YPrediction);
% Rc=  corr(YTrain,YPrediction,'Type','Pearson');