clc, clear all, close all;
% XTrain=xlsread('Compilation_op.xls');
% [XTrain,~,YTrain] = digitTrain4DArrayData;
% [XValidation,~,YValidation] = digitTest4DArrayData;
XTrain = xlsread('1stderv_preprocessed_training_sample_spectra.xlsx');
% XTrain = XTrain';
height = 458;
width = 1;
channels = 1;
samples = 105;
CNN_TrainingData = reshape(XTrain,[height, width, channels, samples]);
YTrain = xlsread('trng ref data_caffeine-35x3_105 sample.xls');
CNN_TrainingLabels = YTrain;

% XValidation = xlsread('1075-1239nm wvl_Best_1st derv.xls');
% XValidation = permute(reshape(XValidation,459,1,1, 100),[1 2 3,4]);
% YValidation = xlsread('CarbonReplicantRefdata.xlsx');
numTrainImages = numel(YTrain);
% [tr,ts]=partition_cv(20,105,1);
% idx = randperm(numTrainImages,20);
% Add layers
layers = [
    imageInputLayer([height, width, channels])

    convolution2dLayer([5 1],100, 'stride',5)
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
    fullyConnectedLayer(100)
    fullyConnectedLayer(50)
    fullyConnectedLayer(20)
    fullyConnectedLayer(5)
    fullyConnectedLayer(1)
    regressionLayer];
RMSECV = [];
R = [];
y_pred =[];
miniBatchSize  = 5;
validationFrequency = floor(numel(YTrain)/miniBatchSize);
options = trainingOptions('sgdm', ...
    'MiniBatchSize',miniBatchSize, ...
    'MaxEpochs',50, ...
    'InitialLearnRate',1e-4, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.1, ...
    'LearnRateDropPeriod',100);
            
            net = trainNetwork(CNN_TrainingData,CNN_TrainingLabels,layers,options);
%             X_test= X_test';
%             X_test = permute(reshape(X_test,458,1,1,1),[1 2 3,4]);
            YPrediction = predict(net,CNN_TrainingData);
%             y_pred=[y_pred ; YPrediction];
            [~,rmsep,~,~]= performance(YTrain,YPrediction);
%              RMSECV=[RMSECV rmsep];
            R1=  corr(YTrain,YPrediction,'Type','Pearson');
%             R = [R ; R1];

% predictionError = YValidation - YPredicted;
% thr = 10;
% numCorrect = sum(abs(predictionError) < thr);
% numValidationImages = numel(YValidation);
% 
% accuracy = numCorrect/numValidationImages;
% squares = predictionError.^2;
% rmse = sqrt(mean(squares));