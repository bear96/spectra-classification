clc, clear all, close all;
XTrain = xlsread('SNV+1stDERV.xlsx'); %reading training spectral data
testdata = xlsread('SNV_1ST_DERV_test.xls'); %reading testing spectral data
YTest = xlsread('test ref data_caffeine_7x3_15 samples.xls'); %test reference data

XXTR = XTrain'; %transposing before resampling to avoid false data acquisition

height = 3000; %new data height; previous height was 459
width = 1;
channels = 1;
samples = 105;

XXTrain = ResampledData(XTrain,height);
XXTest = ResampledData(testdata,height);

CNN_TrainingData = reshape(XXTR,[height,width,channels, samples]); %changing the train spectral data into 4-D data for the Image Input Layer
XXTest = reshape(XXTest,[height,width,channels, 15]); %changing the test spectral data into 4-D data
YTrain = xlsread('trng ref data_caffeine-35x3_105 sample.xls'); %train reference data
CNN_TrainingLabels = YTrain; %storing train reference data in a matrix

%defining layers from here
layers = [
    imageInputLayer([height,width, channels])

    convolution2dLayer([5 1],100, 'stride',1)
    batchNormalizationLayer
    reluLayer
    
    averagePooling2dLayer([50 1],'Stride',2)

     convolution2dLayer(3,16, 'Padding', 'same')
     batchNormalizationLayer
     reluLayer
     
     averagePooling2dLayer(1,'Stride',1)
  
     convolution2dLayer(3,32, 'Padding','same')
     batchNormalizationLayer
     reluLayer
     
     convolution2dLayer(3,32, 'Padding','same')
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
%layer details ends
miniBatchSize  = 3; %setting mini batch size for the solver to reach minimum gradient
% validationFrequency = floor(numel(YTrain)/miniBatchSize);
%defining training options:
options = trainingOptions('sgdm', ...
    'MiniBatchSize',miniBatchSize, ...
    'MaxEpochs',30, ...
    'InitialLearnRate',1e-3, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.1, ...
    'LearnRateDropPeriod',100);
            
net = trainNetwork(CNN_TrainingData,CNN_TrainingLabels,layers,options); %training the network
YPrediction = predict(net,XXTest); %testing the network using test data
[~,rmsecv,~,~]= performance(YTest,YPrediction); %finding rmsep
Rc=  corr(YTest,YPrediction,'Type','Pearson');%finding correlation 