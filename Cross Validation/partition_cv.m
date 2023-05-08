%this forms the mask of the partitions eg: X(ts(2,:),:) will give test set of
%2nd partition
function [tr,ts]=partition_cv(m,n,q)
ts=false(m,n);
tr=false(m,n);
c=1;
% nn=n;
for i=1:q:n
if(i==1)
ts(c,i:q)=true;
tr(c,:)=~ts(c,:);
else
ts(c,i:i+q-1)=true;
tr(c,:)=~(ts(c,:));
end
c=c+1;
% n=n+5;
end