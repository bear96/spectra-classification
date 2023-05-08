function XXTrain = ResampledData(XTrain,samplingrate)
% XTrain = xlsread('MSC+1st DERV.xlsx');
l = length(XTrain(:,1)); %finding length of the data
for i = 1:length(XTrain(1,:))
    s = XTrain(:,i);
    s = s';
%     time = 0:1/l:(l-1)/l;
    y = spline(s, samplingrate ,l); %applying cubic spline interpolation to resample data
    XXTrain(:,i) = y'; %storing resampled data into variable matrix
end
end