function [SNV,MSC]=func_ip2(op1)

B=op1;
[a z]=size(B);
% SNV PROGRAM

msnv=mean(B);
stdv=std(B);
for i=1:a
for j=1:z
        
        SNV(i,j)=(B(i,j)-msnv(j))/stdv(j);                                
end
end

% MSC PROGRAM
[m n]=size(B');
rs=mean(B');
cw=ones(1,n);
mz=[];
mz=[mz ones(1,n)'];
mz=[mz rs'];
[mm,nm]=size(mz);
wmz=mz.*(cw'*ones(1,nm));
wz=B'.*(ones(m,1)*cw);
z=wmz'*wmz;
[u,s,v]=svd(z);
sd=diag(s)'; 
cn=10^12;
ms=sd(1)/sqrt(cn);
cs=max(sd,ms ); 
cz=u*(diag(cs))*v';  
zi=inv(cz);
b=zi*wmz'*wz';
C=b';
x_msc=B; 
p=C(:,1);
x_msc=x_msc-(p*ones(1,mm))';
p=C(:,2);
x_msc=x_msc./(p*ones(mm,1)')';
MSC=x_msc;

end