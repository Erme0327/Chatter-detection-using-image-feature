% =====Feature assessment=====
% Developed by Yun Chen, Xiamen University
% Please cite the paper: [1] Chen Y, Li H, Hou L, Bu X. Feature extraction 
% using dominant frequency bands and time-frequency image analysis for 
% chatter detection in milling. Precision Engineering. 2018.

clc;clear;close all
load('FstSt2013_section831_Ver.mat');
load( 'label_Acc2013_m1Ver1.mat');
Y=[l0903;l0904; l0913;];% Generate the classe label array
Y(Y==0)=1;%Unclear tests to stable tests

X=features';
NumFeas=size(features,2);
J12=zeros(NumFeas,1);
for cnt1=1:NumFeas
    class1=X(cnt1,Y==1);
    class2=X(cnt1,Y==-1);
%     figure(1);
%     plotHist_yun(class1',class2');
    classlabels = [1*ones(length(class1),1); -1*ones(length(class2),1)];
%     figure;
    J12(cnt1)=ROC([class1 class2]',classlabels,0);
end
[~,inJ]=sort(J12, 'descend');% Re-order J12

feasPbnd=12;
temp=zeros(feasPbnd,5);
for cnt=1:5
    temp(:,cnt)=J12(1+(cnt-1)*feasPbnd:cnt*feasPbnd);
%     axis tight;
end
figure;
bar(temp);
xlabel('Feature number','fontsize',14);
ylabel('AUC');
set(gca,'FontSize', 14,'Fontname','Times new roman');
legend('1^s^t band','2^n^d band','3^r^d band','4^t^h band','5^t^h band');
%%
figure;
subplot(121);
plot(features(Y==1,inJ(1)),features(Y==1,inJ(2)),'bo', ...,
    features(Y==-1,inJ(1)),features(Y==-1,inJ(2)),'r*');
xlabel('Mean Correlation - 5th Band','fontsize',14);
ylabel('Correlation range - 5th Band');
set(gca,'FontSize', 14,'Fontname','Times new roman');
subplot(122)
plot(features(Y==1,inJ(end-1)),features(Y==1,inJ(end)),'bo', ...,
    features(Y==-1,inJ(end-1)),features(Y==-1,inJ(end)),'r*');
ylabel('Energy range - 3rd Band','fontsize',14);
xlabel('Homogeneity range  - 1st Band');
set(gca,'FontSize', 14,'Fontname','Times new roman');
inJ(end),inJ(end-1)