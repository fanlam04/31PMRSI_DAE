function [ output_args] = my_forward_Op_1layer( input, para )
%FORWARD_OP Summary of this function goes here
%   Detailed explanation goes here
% input: N_sample by num_spectra
Ns = size(input,1);

%x = [real(input) imag(input)];
x  = input;
y1 = x*para.weightsl1 + repmat(para.biasl1, [Ns, 1]);
y1 = relu(y1);

y2 = y1*para.weightsl2 + repmat(para.biasl2, [Ns, 1]);

output_args = y2;

end


function [output] = relu(input)
    index = find(input>0);
    output = zeros(size(input));
    output(index) = input(index);
end
