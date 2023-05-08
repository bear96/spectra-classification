function [mnmax,detr,MC]=func_ip1(op1)

B=op1;
[a z]=size(B);
%MEAN MAX PROGRAM

for i=1:z
    r(i)=min(B(:,i));
    h(i)=max(B(:,i));
end

 for j=1:a
   for i=1:z
        mnmax(j,i)=(B(j,i)-r(i))/(h(i)-r(i));     
   end
 end
   
% DETRENDED DATA PLOT

detr=detrend(B);

% % MEAN-CENTERING DATA PLOT

m=mean(B);
%s=std(A);
for i=1:a
    for j=1:z
        MC(i,j)=(B(i,j)-m(j));
       
    end
end

end