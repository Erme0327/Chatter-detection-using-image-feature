function [ features] = FstFeas( A ,N_gray)
%Normalization of the pixels intensity in [0,N_gray-1]
A=round((N_gray-1)*((A-min(A(:)))/(max(A(:))-min(A(:)))));
features(1)=mean2(A);
features(2)=std2(A);
features(3)=skewness(A(:));
features(4)=kurtosis(A(:));
end


