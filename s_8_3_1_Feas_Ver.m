% =====Feature extraction=====
% Developed by Yun Chen, Xiamen University
% Please cite the paper: [1] Chen Y, Li H, Hou L, Bu X. Feature extraction 
% using dominant frequency bands and time-frequency image analysis for 
% chatter detection in milling. Precision Engineering. 2018.

clc;clear;close all
% Exp data 20130903
dir_fm1='....\Exp data 20130903\depth 8\';

dir_fm2='...\Exp data 20130903\depth 12\';

% Exp data 20130904
dir_fm3='...\Experimental data 20130904\std data\';

% % Exp data 20130913
dir_fm4='...\Exp data 20130913\std\';

dir_fm_a={dir_fm1,dir_fm2,dir_fm3,dir_fm4};
% FigFeas_a={FigFeas1,FigFeas2,FigFeas3,FigFeas4};

%% Feature evaluation
fs=2e4;
seg=4; % length of segment (unit s)
n_seg=round(seg*fs);

windowsize =512;%round(60/omega*fs*cnt_wd);
window = hanning(windowsize);
nfft= windowsize;
noverlap= round(windowsize/4);

ToNumTest=82;% number of total tests
k=1;% order of tests

freArray=[180 350;1046 1454;1654 1885;2800 3050;4800 5200];
num_fre=size(freArray,1);
perFeas=12;
features=zeros(ToNumTest,perFeas*num_fre);
N_gray=2^8;

for cnt_data=1:4
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
    
    for cnt0=1:Num_file %cnt_test=1:length(test_num)%
        Exp_or=In_or(cnt0);
        name=filename(Exp_or).name;
        data=load([dir_fm name]);    %Load data
        
        %extract a segmentaion of data defined by the variable seg
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
        x=data.Accelerationstdy(in_seg,3);  % Acc data in y direction
        
        % Short-time frequency transform
        [S,F,T] = spectrogram(x,window,noverlap,nfft,fs);
        for cnt_fre=1:num_fre
            Fre_low=freArray(cnt_fre,1);% frequency lower limit Hz
            Fre_upper=freArray(cnt_fre,2);% frequency upper limit Hz
            
            tr_low=round(Fre_low/max(F)*length(F));
            tr_up=round(Fre_upper/max(F)*length(F));
            Ftr=F(tr_low:tr_up);
            Str=S(tr_low:tr_up,:);
            Str=abs(Str);            

            %First-order statistics feature extraction
            features(k,1+(cnt_fre-1)*perFeas:4+(cnt_fre-1)*perFeas)= FstFeas(Str,N_gray );
            % Second-order features
            features(k,5+(cnt_fre-1)*perFeas:12+(cnt_fre-1)*perFeas) = SndFeas(Str,N_gray );
            
        end % for cnt_fre=1:num_fre
        k=k+1;
    end% for cnt0=1:Num_file 
end%end of cnt_data

save('FstSt2013_section831_Ver.mat','features');
