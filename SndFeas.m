function [features] = SndFeas( A,N_gray )
%computes second-order statistics of an image
%Normalization of the pixels intensity in [0, N_gray-1]
A=round((N_gray-1)*((A-min(A(:)))/(max(A(:))-min(A(:)))));
[glcm,~]=graycomatrix(A,'GrayLimits',[0,N_gray-1],...
    'NumLevels',N_gray,'Offset',[0 1;-1 0;-1 1; -1 -1],'Symmetric',true);
stats=graycoprops(glcm,{'Contrast','Correlation',...
    'Energy','Homogeneity'});
features(1)=mean(stats.Contrast);
features(2)=mean(stats.Correlation);
features(3)=mean(stats.Energy);
features(4)=mean(stats.Homogeneity);
features(5)=range(stats.Contrast);
features(6)=range(stats.Correlation);
features(7)=range(stats.Energy);
features(8)=range(stats.Homogeneity);
end