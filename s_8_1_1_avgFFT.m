% =====Dominant frequency bands=====
% This code is used to examine the effect of different tests selected
% on the sythethized FFT

% Developed by Yun Chen, Xiamen University

% Refer to the paper: [1] Chen Y, Li H, Hou L, Bu X. Feature extraction 
% using dominant frequency bands and time-frequency image analysis for chatter detection in milling. 
% Precision Engineering, 2019, 56: 235-245.


clc;clear;close all
% Exp data 20130903
dir_fm1='...\Exp data 20130903\depth 8\';
Y1=[-1	-1	-1	1	1	1	1	1	1	1	1	1	1	1];%Test label; -1 stands for unstable tests; 1 stands for stable tests
dir_fm2='...\Exp data 20130903\depth 12\';
Y2=[-1	-1	1	-1	-1	-1	1	1	1	1	1	1	1	1];

% Exp data 20130904
dir_fm3='...\Experimental data 20130904\std data\';
Y3=[1 	1 	1 	1 	1 	1 	-1	-1 	-1 	1 	-1	-1 	1 	1 	-1 	-1 	-1 	1	-1 	-1 	-1 	1	1	1	1	1	-1 	-1 	1];

% % Exp data 20130913
dir_fm4='...\Exp data 20130913\std\';
Y4=[1	-1	-1	1	1	1	1	1	1	1	1	-1	1	1	-1	-1	-1	1	-1	1	1	1	-1	-1	-1];

dir_fm_a={dir_fm1,dir_fm2,dir_fm3,dir_fm4};
Y=[Y1 Y2 Y3 Y4];

%% Feature evaluation
fs=2e4;
T=1/fs; %sample T s
seg=1.5; % length of segment (unit s)
n_seg=round(seg*fs);

rng('default');
TestOD=1:82;
SdTest=TestOD(Y==1);
UsdTest=TestOD(Y==-1);
%%
for cntR=1:3
    p1 = randperm(length(SdTest),5);
    in1=SdTest(p1);
    p2 = randperm(length(UsdTest),5);
    in2=UsdTest(p2);
    p=[in1 in2];
    Pre_a=zeros(length(p),2);
    for cnt_p=1:length(p)
        if p(cnt_p)<=14
            Pre_a(cnt_p,1)=1;
            Pre_a(cnt_p,2)=p(cnt_p);
        elseif  p(cnt_p)<=28
            Pre_a(cnt_p,1)=2;
            Pre_a(cnt_p,2)=p(cnt_p)-14;
        elseif p(cnt_p)<=57
            Pre_a(cnt_p,1)=3;
            Pre_a(cnt_p,2)=p(cnt_p)-28;
        else
            Pre_a(cnt_p,1)=4;
            Pre_a(cnt_p,2)=p(cnt_p)-57;
        end
    end
    % Pre_a=[2 1;2 2;2 6;
    %     1 1;2 4;1 2;1 3;
    % %     1 8;1 9;4 7;4 8;
    %     1 12; 1 13; 1 14;];% Pre_a is a matrix for the selected tests; the first
    % column is the number of the experimental subset;
    %                        the second column is the test number in the
    %                        correpsonding experimental subset.
    Num_Pre=size(Pre_a,1);% Total number of tests selected for average FFT
    
    % n = 2^nextpow2(n_seg);
    % sum_P1=zeros(n/2+1,1);
    x_constr=[];
    f_res=2; %Hz
    k=1;
    for cnt_k=1:Num_Pre
        cnt_data=Pre_a(cnt_k,1);% Date order for the subset of the tests
        dir_fm=dir_fm_a{cnt_data};
        
        clear filename;
        filename=dir([dir_fm '*.mat']);   %load all mat files
        
        Num_file=length(filename);
        % Sort the filename order according to the experiment order
        clear C C1 C2;
        C(1:Num_file)={'.mat'};C1(1:Num_file)={''};
        C2=cellfun(@strrep, {filename.name},C,C1,'UniformOutput',false);
        C(:)={'whole'};
        Order_or=str2double(cellfun(@strrep, C2,C,C1,'UniformOutput',false));
        [~,In_or]=sort(Order_or);
        
        for cnt_file=1%:Num_file
            cnt0=Pre_a(cnt_k,2);% test order
            Exp_or=In_or(cnt0);%obtain the real experimental order
            name=filename(Exp_or).name; %Double check the selected test names
            
            data=load([dir_fm name]);    %Load data
            %extract part of data
            [l_sgn,~]=size(data.Accelerationstdy);
            temp=round((l_sgn-n_seg)/2);temp1=temp+n_seg-1;
            if (temp>0) && (temp1<=l_sgn)
                in_seg=temp:temp1; % Extract the data from the middle of data.Accelerationstdy
            elseif (temp<=0) && (temp1<=l_sgn)
                in_seg=1:temp1;
            elseif (temp>0) && (temp1>l_sgn)
                in_seg=temp:l_sgn;
            else
                in_seg=1:l_sgn;
            end
            
            clear x P_org F_org;
            x=data.Accelerationstdy(in_seg,3);  % Acc data in y direction
            [P_org, deltaf_act]=fft_yun(x,f_res, 6e3, fs); % a simple fft function can be used here instead
            n_data_need= length(P_org); F_org=(1:n_data_need)*deltaf_act;
            if k==1
                x_constr=P_org;
            else
                x_constr=x_constr+P_org;
            end
            k=k+1;
            %     figure(cnt_k);
            % %     subplot(211);
            %     plot(time,x,'k','LineWidth',1);
            %     xlabel('Time (s)','fontsize',14);
            %     ylabel('Vibration (m/s^2)');
            %     set(gca,'FontSize', 14);
            %     axis tight;
        end
    end %end of cnt_omega

%     %
%     % FFT
%     figure;
%     plot(F_org ,x_constr,'k','LineWidth',1);
%     xlabel('Frequency (Hz)','fontsize',14);
%     ylabel('Magnitude (m/s^2)');
%     set(gca,'FontSize', 14,'Fontname','Times new roman');
%     title(['Freq Res=' num2str(f_res) ' Hz']);
    %%
    X0=[50 200  400  600 800  1000];
    for cnt_wL=1:length(X0)
        temp0=X0(cnt_wL);
        w = hanning(temp0);
        winL(cnt_wL)=temp0;
        temp=conv(x_constr.^2,w.^2,'same');
        SF(:,cnt_wL)=temp;
        F_sf=(0:length(SF)-1)*deltaf_act;
        y=temp0*ones(1,length(SF));
        figure(3);plot3(F_sf,y,(temp/max(temp)),'linewidth',1);hold on;
    %     title(['Win Size=' num2str(winL)]);
    end
    box on;
    xlabel('Frequency (Hz)','fontsize',14);
    ylabel('Window szie');zlabel('Normlized squared energy');
    set(gca,'FontSize', 14,'Fontname','Times new roman');

    %%
    w = hanning(400);
    temp=conv(x_constr.^2,w.^2,'same');
    temp=temp/max(temp);
    figure(4);plot(F_org,temp,'linewidth',1.5);hold on;
    xlabel('Frequency (Hz)','fontsize',14);
    ylabel('Normlized squared energy');
    set(gca,'FontSize', 14,'Fontname','Times new roman');    
    
    freArray=[0 500;670 1470;1550 2100;2700 3350;4600 5600];    
    num_fre=size(freArray,1);
    peak=zeros(num_fre,1);
    F_res=F_org(2)-F_org(1);
    line_lgth=round(300/F_res);
    for cnt_fre=1:num_fre
        Fre_low=freArray(cnt_fre,1);% frequency lower limit Hz
        Fre_upper=freArray(cnt_fre,2);% frequency upper limit Hz
        ind=find((F_org>Fre_low)&(F_org<Fre_upper));
        [peak(cnt_fre),indmax]=max(temp(ind));
        in_st=ind(indmax)-line_lgth;
        if in_st<=0
            in_st=1;
        end
        in_array=in_st:(ind(indmax)+line_lgth);
        line_x=F_org(in_array);
        line_y=peak(cnt_fre)/sqrt(2)*ones(size(line_x));
        plot(line_x,line_y,'r-.','linewidth',1)
    end
end