function [map, rmse_n, rmse, mae] = performance(act , pred)
mae=0;
for i=1:length(pred)
    
        mse(i)=((act(i)-pred(i)))^2;
        msen(i)=((act(i)-pred(i))/act(i))^2;
        mae=mae+abs(((act(i)-pred(i))));
end
count=0;
for i=1:length(pred)
    
        count=count+abs((abs(act(i)-pred(i))/act(i)));
end

map=(count/length(act))*100;
rmse=sqrt(sum(mse)/(length(mse)));
rmse_n=sqrt(sum(msen)/(length(msen)));
mae=mae/length(act);


end

