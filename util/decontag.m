function [RI] = decontag(input)
% Reforming the concatenated real and imag parts to complex FIDs (RI)

RI = complex(input(:,1:size(input,2)/2),input(:,size(input,2)/2+1:end));