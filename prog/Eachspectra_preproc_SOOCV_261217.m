close all;
clear all;
clc;

% DECLAIRARATION PART

nirspectra='I:\Steps_Model building\Soil Model Building\org Carbon\data\quartz plate_light from bottom\Calibration with quartz\120417\';    % CHECK--
[Refdata]=xlsread('I:\Steps_Model building\Soil Model Building\org Carbon\ref data\Refdata_curbon.xls');               % CHECK--
filename_string='S%d_S%d_%d.ABS';         % CHECK--3 for loop
   
samp_nocompile=20;
samp_no=size(Refdata,1);                 % 'samp_no' corresponds to no of sample.
pos_list=1;%1:3;                             % CHECK--position er name 'A'-'E' declaire kara hache 'pos_list' name
replicant=1:5;     

startwvl1=924.50;                            % CHECK--   929.75, 980.50, 1020.75, 1080.25, 1209.75, 1398.75, 1500.25, 1580.75
endwvl1=1045.25;                            % CHECK--  970.00, 1020.75,1071.50, 1129.25, 1250.00, 1440.75, 1540.50, 1621.00   

startwvl2=1099.50;                            % CHECK--   929.75, 980.50, 1020.75, 1080.25, 1209.75, 1398.75, 1500.25, 1580.75
endwvl2=1222.00;                            % CHECK--  970.00, 1020.75,1071.50, 1129.25, 1250.00, 1440.75, 1540.50, 1621.00   

startwvl3=1390.00;                            % CHECK--   929.75, 980.50, 1020.75, 1080.25, 1209.75, 1398.75, 1500.25, 1580.75
endwvl3=1600.00;                            % CHECK--  970.00, 1020.75,1071.50, 1129.25, 1250.00, 1440.75, 1540.50, 1621.00   

spectrano_sample=(length(pos_list)*length(replicant)); % CHECK-- case 1
%spectrano_sample_rep=15;
total_spec=samp_no*spectrano_sample;

startrow1=((startwvl1-879)/1.75)+1;
endrow1=((endwvl1-879)/1.75)+1;

startrow2=((startwvl2-879)/1.75)+1;
endrow2=((endwvl2-879)/1.75)+1;

startrow3=((startwvl3-879)/1.75)+1;
endrow3=((endwvl3-879)/1.75)+1;



% clmno_55sample=1:(55*spectrano_sample);
% clmno_61to66sample=(60*spectrano_sample+1):(66*spectrano_sample);

All_abs=[];                            %  All_abs name akti matrix er maddhe jato bar for loop chalbe data jama hache.
All_abs1=[];
Allabs_avg=[];                        % 'sample_spectra' name akta matrix creat kara hache 
Refdata_rep=[];


% COMPILATION OF RAW DATA PART

for sample_id=1:samp_nocompile 
for pos=1:length(pos_list)                                                    % 'pos_list' name list a jato gulo alphabet ache tato bar 'pos' name for loop(for loop no-2) ti chalbe
for repl=1:length(replicant)                                                            % rep = replica // j hetu pratita position a 3 te replica tai 'rep' name akti for loop 3 bar chalano hache (for loop no-3)
            
            filename=sprintf(filename_string,sample_id,pos_list(pos),repl); 
            filepath=strcat(nirspectra,filename);                            % 'strcat' syntax ti anekta excell er concatenation er mato kaj kare. akhane 'stdir' a declaire kara path er sathe filename ta concatenate karche. fale natun file path ta "filepath" name save hache. Ai "filepath" nam ta function a babohar hache.
            import_file=importdata(filepath);
            Data_infile=import_file.data;                                    % storing data in Data_infile
            wavel1=Data_infile(startrow1:endrow1,1);                            % akta 'wavel' name data set tairi kara hache jar 1st row ta holo 'Data_infile' name file tir 1st row. ar shes row ta determined hache 'to' name parameter dara. ar column naoa hayeche 'data_infile' name file tir pratham column take.
            wavel2=Data_infile(startrow2:endrow2,1);
            wavel3=Data_infile(startrow3:endrow3,1);
            
            abs1=Data_infile(startrow1:endrow1,2);                              % akta 'abs' name data set tairi kara hache jar 1st row ta holo 'Data_infile' name file tir 1st row. ar shes row ta determined hache 'to' name parameter dara. ar column naoa hayeche 'data_infile' name file tir "2nd" column take.
            abs2=Data_infile(startrow2:endrow2,2);
            abs3=Data_infile(startrow3:endrow3,2);
            
            All_abs1=[abs1; abs2; abs3];                                           % "All_abs" holo akti matrix jar shesh column ta holo current file tar "abs" value ar ager colun gulo purono "All_abs" file. fale pratibar for loop challe All_abs ta upgrade hay o size a barte thake...Ai vabe 3 bar(jehetu rep=replicant) chalar par natun sample er janna for loop suru hay.
            All_abs=[All_abs All_abs1];
end                                                                          % end of (no-3) for loop. 
end                                                                          % end of (no-2) for loop.
end

% AVERAGE CALCULATN - CHECK- eachsamp_spectra should be= length(replicant)

samp_no=size(Refdata,1);           % 'samp_no' is list of no which corresponds to no of sample.
pos_list=1:3;                % CHECK--position er name 'A'-'E' declaire kara hache 'pos_list' name
replicant=1:5;     
All_abs_avg=[];

avgno=spectrano_sample;

Average_ip=All_abs;

for i=1:size(Refdata,1)
k=((i-1)*avgno+1);  

av_abs=All_abs(:,k:k+(avgno-1));    % pratita sample er 3 te replica r janna for loop shesh haoar pare tar average niye 'av_abs' name store kara hache.
mean_abs=mean(av_abs,2);
All_abs_avg=[All_abs_avg mean_abs];      % ager line a j 'av_abs'calculate kara holo 3 te replicar average niye ta 'Allabs_avg' name matrix er last column a update kara hache. fale pratibar 2nad for loop er run er pare matrix ti size a barche. 

end

% SIZE INCREMENT- REF VALUE //  The below portion consider when each replica is considered

for i=1:samp_no
Refdata_rep=[Refdata_rep;Refdata(i).*ones(spectrano_sample,1)];
end
Refdata=Refdata_rep;
% 
% Average_op=All_abs_avg;

% Compilation_op1=All_abs(:,clmno_55sample);
% Compilation_op2=All_abs(:,clmno_61to66sample);
% Compilation_op=[Compilation_op1 Compilation_op2];     

Compilation_op=All_abs;
  
%Preprocess Operation

op1=Compilation_op;
[mnmax,detr,MC]=func_ip1(op1);
[SNV,MSC]=func_ip2(op1);
[a z]=size(op1);
 A1=mnmax;
 A2=detr;
 A3=MC;
 A4=SNV;
 A5=MSC;
 
 detr_new=[];
 mnmax_new=[];
 MC_NEW=[];
 SNV_NEW=[];
 MSC_NEW=[];

 for i=1:5;
    
     if (i==1)
    B= A1;
    
     end
    
    if i==2;
        B=A2;
    end
        if i==3;
    B= A3;
        end
    
    if i==4;
        B=A4;
    end
        
        if i==5;
           B= A5;
        end
        
            
        %MEAN MAX PROGRAM
    [a z]=size(B);
for i=1:z
    r(i)=min(B(:,i));
    h(i)=max(B(:,i));
end

 for j=1:a
   for i=1:z
        mnmax2(j,i)=(B(j,i)-r(i))/(h(i)-r(i));     
   end
 end
 
 mnmax_new=[mnmax_new mnmax2];
 
  % DETRENDED DATA PLOT
  
         detr2=detrend(B);
        
        detr_new=[detr_new detr2];
 
 % % MEAN-CENTERING DATA PLOT

m=mean(B);
%s=std(A);
for i=1:a
    for j=1:z
        MC2(i,j)=(B(i,j)-m(j));
       
    end
end

MC_NEW=[MC_NEW MC2];

% SNV PROGRAM

msnv=mean(B);
stdv=std(B);
for i=1:a
for j=1:z
        
        SNV2(i,j)=(B(i,j)-msnv(j))/stdv(j);                                
end
end

SNV_NEW=[SNV_NEW SNV2];

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
MSC2=x_msc;

MSC_NEW=[MSC_NEW MSC2];
    
 end

 % Data saved
  [m X]=size(B);
  A = 1:X;
  
 % Mean Max

%    minmax=A1;
%    mnmax_mnmax=mnmax_new(:,A);
    %detr_mnmax=mnmax_new(:,A+X);
%    MC_mnmax=mnmax_new(:,A+2*X);
%    SNV_mnmax=mnmax_new(:,A+3*X);
%    MSC_mnmax=mnmax_new(:,A+4*X);
%  
%  
%    
%    % Detrend
%    
% de_trend=A2;
    %mnmax_detr=detr_new(:,A);
    %detr_detr=detr_new(:,A+X);
    %MC_detr=detr_new(:,A+2*X);
    %SNV_detr=detr_new(:,A+3*X);
    %MSC_detr=detr_new(:,A+4*X);
%    
%    % MC
%    
         %MnCnt=A3;
%      mnmax_MC=MC_NEW(:,A);
      %detr_MC=MC_NEW(:,A+X);
%      MC_MC=MC_NEW(:,A+2*X);
%      SNV_MC=MC_NEW(:,A+3*X);
%      MSC_MC=MC_NEW(:,A+4*X);
%    
%    % SNV
%    
%     SdNV=A4;
%     mnmax_SNV=SNV_NEW(:,A);
%       detr_SNV=SNV_NEW(:,A+X);
%     MC_SNV=SNV_NEW(:,A+2*X);
%     SNV_SNV=SNV_NEW(:,A+3*X);
%     MSC_SNV=SNV_NEW(:,A+4*X);
%    
%    % MSC
%    
      MultSC=A5;
%      mnmax_MSC=MSC_NEW(:,A);
      %detr_MSC=MSC_NEW(:,A+X);
%      MC_MSC=MSC_NEW(:,A+2*X);
%      SNV_MSC=MSC_NEW(:,A+3*X);
%      MSC_MSC=MSC_NEW(:,A+4*X);
        
 Preprocs_op=MultSC;    
 
% PLS PROGRAM 

% Declairation part of PLS program

PLS_ip=Preprocs_op;

rownum_PLS_ip=size(PLS_ip,1);
clmnum_PLS_ip=size(PLS_ip,2); 

nirdata_cdbl=[PLS_ip PLS_ip];
Refdata_rdbl=[Refdata;Refdata];

actualcarbn_matx=[];
y1=[];
pred_traincarbn_matx=[];
pred_testcarbn_matx=[];
CF_matrix=[];
CF_matrix_final=[];
CF_matrix_compt_new=[];
pred_testcarbnmatrix_compt_new=[];
pred_testcarbn_matx_compt_new=[];
actualcarbn_matx=[];
actualcarbn_matx_compt_new=[];
Accuracy_comput_A_new=[];
final_data_computd=[];
final_accurc_computd=[];
total_result=[];
final_pred=[];
final_pred_computd=[];
sqr_error=[];

num_samp=size(Refdata,1);
replic=1;%15;
trngset=1:(clmnum_PLS_ip-replic);
component_number=1:20;%length(trngset)-1;
%(size(Allabs_avg,2)/size(Refdata,1));



% LOOP FORMATION 

%'ncomp' = Component Number and 'j' = Index Number of program run(not sample number.In case of j=1,test sample=24)
for ncomp=component_number 
for j=1:clmnum_PLS_ip 

      k=((j-1)*replic+1);

     %k=((ncomp-1)*replic+1);

% Data selection automatically to run PLS program

train_nirdata= nirdata_cdbl(:,k:(k-1)+(length(trngset)));
trainref_carbnval=Refdata_rdbl(k:(k-1)+(length(trngset)),1);
test_nirdata= nirdata_cdbl(:,k+length(trngset):(k-1)+clmnum_PLS_ip);  
testref_carbnval=Refdata_rdbl(k+length(trngset):(k-1)+(clmnum_PLS_ip),1);
actualcarbn_matx=[actualcarbn_matx; testref_carbnval];

X=[train_nirdata'];                                                    % inverse kare row ar column interchange kara hache o natun matrix ti "X" bale nam daoa hache

y=trainref_carbnval;                                                          % "Refdata" name matrix take abar "y" name declaire kara hache
y1=[y1 y];
y1_t=[y1'];


X_test=[test_nirdata'];                                                % inverse kare row ar column interchange kara hache o natun matrix ti "X_test" bale nam daoa hache


% Main PLS program

[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(X,y,ncomp);           % eta PLS er standard formulla ba command
pred_traincarbn=[ones(size(X,1),1) X]*BETA;
pred_traincarbn_matx=[pred_traincarbn_matx pred_traincarbn];

pred_testcarbn=[ones(size(X_test,1),1) X_test]*BETA;
pred_testcarbn_matx=[pred_testcarbn_matx;pred_testcarbn];

% display Predicted value and Correlation factor

% disp('===================================================');
% ncomp
% predicted_testcarbn=pred_testcarbn
% 
% Actualvalue_testcarbn=mean(testref_carbnval)
% disp('===================================================');

end
end

% CF CALCULATION

for l=1:(samp_no*(length(component_number))) 
%for l=1:(clmnum_PLS_ip*(length(component_number)))       % 1:Number of rows in the matrix of predicted value of training set(e.g - pred_traincarbn_matx)     
   
s=y1(:,l);
r=pred_traincarbn_matx(:,l);
CF= corr(s,r);
CF_matrix=[CF_matrix CF];

end;

Z=CF_matrix;

for j=1:length(component_number)
    l=((j-1)*samp_no+1);
    %l=((j-1)*clmnum_PLS_ip+1);
    
    CF_matrix_compt=Z(1,l:l+(samp_no-1));
    %CF_matrix_compt=Z(1,l:l+(clmnum_PLS_ip-1));
    CF_matrix_compt_new=[CF_matrix_compt_new; CF_matrix_compt];
    
    k=((j-1)*clmnum_PLS_ip+1);
    
    pred_testcarbn_matx_compt=pred_testcarbn_matx(k:k+(clmnum_PLS_ip-1),1);
    pred_testcarbn_matx_compt_new=[pred_testcarbn_matx_compt_new pred_testcarbn_matx_compt];
    
    actualcarbn_matx_compt=actualcarbn_matx(k:k+(clmnum_PLS_ip-1),1);
    actualcarbn_matx_compt_new=[actualcarbn_matx_compt_new actualcarbn_matx_compt];
           
end

CF_matrix_compt_new_t=[CF_matrix_compt_new'];

for j=1:samp_no%clmnum_PLS_ip
    
    CF_final_comp=CF_matrix_compt_new_t(j,:);
    CF_matrix_final=[CF_matrix_final;CF_final_comp];  
    
end



for row=1:clmnum_PLS_ip
for j=1:length(component_number)

if ( pred_testcarbn_matx_compt_new(row,j)> 0 && pred_testcarbn_matx_compt_new(row,j)> actualcarbn_matx_compt_new(row,j))
    
  error(row,j)= (pred_testcarbn_matx_compt_new(row,j)-actualcarbn_matx_compt_new(row,j))/actualcarbn_matx_compt_new(row,j);
 
  
else if (pred_testcarbn_matx_compt_new(row,j)> 0 && pred_testcarbn_matx_compt_new(row,j)< actualcarbn_matx_compt_new(row,j))
      
  error(row,j)= (actualcarbn_matx_compt_new(row,j)-pred_testcarbn_matx_compt_new(row,j))/actualcarbn_matx_compt_new(row,j);
      
    
else        
   error(row,j)= (pred_testcarbn_matx_compt_new(row,j)-actualcarbn_matx_compt_new(row,j))/actualcarbn_matx_compt_new(row,j); 
 
end  
end
         
  sqr_error(row,j)= error(row,j)*error(row,j);%((pred_testcarbn_matx_compt_new(row,j)-actualcarbn_matx_compt_new(row,j))*(pred_testcarbn_matx_compt_new(row,j)-actualcarbn_matx_compt_new(row,j)));
  
  Accuracy_perc(row,j)=(1-error(row,j))*100;
 
end
end


RMSECV=sqrt(sum(sqr_error,1)/length(trngset));

for j=1:length(component_number)
for row=1:samp_no%clmnum_PLS_ip
    k=((row-1)*replic+1);
    
    Avg_Accuracy(row,j)= mean(Accuracy_perc(k:k+(replic-1),j));
end
end


for p=1:length(component_number)
final_data= [actualcarbn_matx_compt_new(:,p) pred_testcarbn_matx_compt_new(:,p) Accuracy_perc(:,p)]; %Avg_Accuracy(:,p) CF_matrix_final(:,p)];
%final_accurc= [Accuracy_perc(:,p) Avg_Accuracy(:,p)];

final_data_computd=[final_data_computd  final_data];
%final_accurc_computd=[final_accurc_computd  final_accurc];

final_pred= [actualcarbn_matx_compt_new(:,p) pred_testcarbn_matx_compt_new(:,p)];
final_pred_computd= [final_pred_computd final_pred];
end

comparison_egty=Accuracy_perc(:,1:length(component_number));
comparison_nty=Accuracy_perc(:,1:length(component_number));


for row=1:size(Accuracy_perc)%size(Avg_Accuracy,1)
for j=1:length(component_number)
if (comparison_egty(row,j)>=80 && comparison_egty(row,j)<=100)
    comparison_egty1(row,j)=1;
           
else 
    comparison_egty1(row,j)=0;
    
end
if (comparison_nty(row,j)>=90 && comparison_nty(row,j)<=100)
   
    comparison_nty1(row,j)=1;       
else 
    
    comparison_nty1(row,j)=0;
end
end
end

% Calculation- over & under fitting

spreading_pos=comparison_nty;
spreading_neg=comparison_nty;

%error_1=[];
for row=1:size(Accuracy_perc)%size(Avg_Accuracy,1)
for j=1:length(component_number)
if (spreading_pos(row,j)>0 && pred_testcarbn_matx_compt_new(row,j)>actualcarbn_matx_compt_new(row,j))
    
    spreading_pos(row,j)=1;
else 
    spreading_pos(row,j)=0;
    
end

if(spreading_neg(row,j)>0 && pred_testcarbn_matx_compt_new(row,j)<actualcarbn_matx_compt_new(row,j))
    
        spreading_neg(row,j)=-1;
else 
    spreading_neg(row,j)=0;
end

end
end

% Results

total_egty=sum(comparison_egty1,1);
total_nty=sum(comparison_nty1,1);

perc_egty=[];
for compno=1:length(component_number)
    
   Acc_perc=total_egty(compno)/size(PLS_ip,2)*100;
   perc_egty=[perc_egty Acc_perc];
end

perc_ninty=[];
for compno=1:length(component_number)
    
   Acc_perc=total_nty(compno)/size(PLS_ip,2)*100;
   perc_ninty=[perc_ninty Acc_perc];
end

total_pos=sum(spreading_pos,1);
total_neg=sum(spreading_neg,1);

total_spread=[total_pos;total_neg;RMSECV];
total_result=[total_egty;total_nty;perc_egty;perc_ninty;total_pos;total_neg;RMSECV];    

% 10% marginal data avg

% Tenperc_ip_actual=actualcarbn_matx;
% Tenperc_ip_pred=pred_testcarbn_matx;
% 
% avg_testcarbn_matx_compt2=[];
% upper_margin_rep=[];
% lower_margin_rep=[];
% 
% for i=1:samp_no*length(component_number)
%     
%         l=((i-1)*spectrano_sample+1);
%         
%         avg_testcarbn_matx_compt=mean(Tenperc_ip_pred(l:l+(spectrano_sample-1),1));
%         
%         avg_testcarbn_matx_compt2=[avg_testcarbn_matx_compt2;avg_testcarbn_matx_compt];
%         
%         upper_margin= avg_testcarbn_matx_compt+(avg_testcarbn_matx_compt*0.1);
%         lower_margin= avg_testcarbn_matx_compt-(avg_testcarbn_matx_compt*0.1);
%         
%         upper_margin_rep=[upper_margin_rep;upper_margin.*ones(spectrano_sample,1)];
%         lower_margin_rep=[lower_margin_rep;lower_margin.*ones(spectrano_sample,1)];
%               
% end
% 
% avg_actualcarbn_final=[];
% upper_margin_compt_new=[];
% lower_margin_compt_new=[];
% 
% 
% 
% for row=1:total_spec*length(component_number)
% 
% if ( Tenperc_ip_pred(row,1)> lower_margin_rep(row,1) && Tenperc_ip_pred(row,1)< upper_margin_rep(row,1))
%     
%   pred_testcarbn_matx_avg(row,1)= Tenperc_ip_pred(row,1);
%   
%   pred_testcarbn_matx_avg_count(row,1)=1;
%  
% else pred_testcarbn_matx_avg(row,1)=0;
%      pred_testcarbn_matx_avg_count(row,1)=0;    
% end
% end
% 
% avg_tenpercent_comp=[];
% avg_tenpercent_comp_pred=[];
% avg_tenpercent_comp_actual=[];
% pred_testcarbn_matx_avg_tenperc_matrix=[];
% Reading_pred_total=[];
% 
% for i=1:samp_no*length(component_number)
%  %for row=1:spectrano_sample
%      
%         l=((i-1)*spectrano_sample+1);
%         
%         pred_testcarbn_matx_avg_tenperc=pred_testcarbn_matx_avg(l:l+(spectrano_sample-1),1);
%         
%         pred_testcarbn_matx_avg_tenperc_count=pred_testcarbn_matx_avg_count(l:l+(spectrano_sample-1),1);
%         
%         total_tenperc=sum(pred_testcarbn_matx_avg_tenperc);
%         total_count=sum(pred_testcarbn_matx_avg_tenperc_count);
%         
%         Reading_pred=total_tenperc/total_count;
%         Reading_pred_total=[Reading_pred_total;Reading_pred];
%        
%         avg_tenpercent_actual=mean(Tenperc_ip_actual(l:(l+spectrano_sample-1),1));
%         
%         avg_tenpercent_comp_actual=[avg_tenpercent_comp_actual;avg_tenpercent_actual];
%         
% 
%  end
%  %end
%         
%               
% %end  
% 
% 
% for row=1:samp_no*length(component_number)
% 
% 
% if ( Reading_pred_total(row,1)> 0 && Reading_pred_total(row,1)> avg_tenpercent_comp_actual(row,1))
%     
%  error_avg(row,1)= (Reading_pred_total(row,1)-avg_tenpercent_comp_actual(row,1))/avg_tenpercent_comp_actual(row,1);
%  
%   
% else if (Reading_pred_total(row,1)> 0 && Reading_pred_total(row,1)< avg_tenpercent_comp_actual(row,1))
%       
%   error_avg(row,1)= (avg_tenpercent_comp_actual(row,1)-Reading_pred_total(row,1))/avg_tenpercent_comp_actual(row,1);
%       
%     
% else        
%    error_avg(row,1)= (Reading_pred_total(row,1)-avg_tenpercent_comp_actual(row,1))/avg_tenpercent_comp_actual(row,1); 
%  
% end  
% end
%          
%   %sqr_error(row,j)= error(row,j)*error(row,j);%((pred_testcarbn_matx_compt_new(row,j)-actualcarbn_matx_compt_new(row,j))*(pred_testcarbn_matx_compt_new(row,j)-actualcarbn_matx_compt_new(row,j)));
%   
%   Accuracy_perc_avg(row,1)=(1-error_avg(row,1))*100;
%  
% end
% 
% 
% avg_tenpercent_comp_actual_compt_new=[];
% avg_tenpercent_comp_pred_compt_new=[];
%  Accuracy_perc_avg_compt_new=[];
%  
% for i=1:length(component_number)
%     
%         l=(i-1)*samp_no+1;
%         
%         avg_tenpercent_comp_actual;
%         avg_tenpercent_comp_pred;
%         Accuracy_perc_avg;
%    
%     avg_tenpercent_comp_actual_compt=avg_tenpercent_comp_actual(l:l+(samp_no-1),1);
%     avg_tenpercent_comp_actual_compt_new=[avg_tenpercent_comp_actual_compt_new avg_tenpercent_comp_actual_compt];
%     
%     avg_tenpercent_comp_pred_compt=Reading_pred_total(l:l+(samp_no-1),1);
%     avg_tenpercent_comp_pred_compt_new=[avg_tenpercent_comp_pred_compt_new avg_tenpercent_comp_pred_compt];
%  
%     Accuracy_perc_avg_compt=Accuracy_perc_avg(l:l+(samp_no-1),1);
%     Accuracy_perc_avg_compt_new=[Accuracy_perc_avg_compt_new Accuracy_perc_avg_compt];
% 
%       
% end

% consolidated_avg=[];
% Total_result_avg=[];
% 
% for i=1:length(component_number)
% consolidated_avg=[avg_tenpercent_comp_actual_compt_new(:,i) avg_tenpercent_comp_pred_compt_new(:,i) Accuracy_perc_avg_compt_new(:,i)];
% Total_result_avg=[Total_result_avg consolidated_avg];
% 
% end
% 
% disp('===================================================');
% 
% disp('Actual_carbn  Pred_carbn  Accuracy_percentage');
% 
% Result=[Total_result_avg]
% 
% %Actual_carbn=avg_tenpercent_comp_actual 
% 
% 
% %Actualvalue_testcarbn=avg_tenpercent_comp_actual
% disp('===================================================');
        



