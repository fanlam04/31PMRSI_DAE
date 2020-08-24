function [out] = contag(input,size)

if nargin == 1
    Data_C = imag(input);
    Data_R = real(input);
    out = [Data_R Data_C];
else
    Data_C = imag(input);
    Data_R = real(input);
    out = [Data_R(:,1:size/2) Data_C(:,1:size/2)];
end


