function [ output_args, J] = my_forward_Op_leakyRelu( input, para, alpha )
%FORWARD_OP Summary of this function goes here
%   Detailed explanation goes here
% input: N_sample by num_spectra
Ns = size(input,1);

%x = [real(input) imag(input)];
x = input;
y1 = x*para.weightsl1 + repmat(para.biasl1, [Ns, 1]);
y1 = leakyrelu(y1,alpha);

y2 = y1*para.weightsl2 + repmat(para.biasl2, [Ns, 1]);
y2 = leakyrelu(y2,alpha);

y3 = y2*para.weightsl3 + repmat(para.biasl3, [Ns, 1]);
y3 = leakyrelu(y3,alpha);

y4 = y3*para.weightsl4 + repmat(para.biasl4, [Ns, 1]);

y5 = y4*para.weightsl5 + repmat(para.biasl5, [Ns, 1]);
y5 = leakyrelu(y5,alpha);

y6 = y5*para.weightsl6 + repmat(para.biasl6, [Ns, 1]);
y6 = leakyrelu(y6,alpha);

y7 = y6*para.weightsl7 + repmat(para.biasl7, [Ns, 1]);
y7 = leakyrelu(y7,alpha);

y8 = y7*para.weightsl8 + repmat(para.biasl8, [Ns, 1]);

output_args = y8;%(:,1:N) + 1i*y8(:,N+1:end);

if nargout > 1
    % calculate Jacobian
    
    U7 = leakyGradient(y6*para.weightsl7+para.biasl7,alpha);

    U6 = leakyGradient(y5*para.weightsl6+para.biasl6,alpha);

    U5 = leakyGradient(y4*para.weightsl5+para.biasl5,alpha);

    U3 = leakyGradient(y2*para.weightsl3+para.biasl3,alpha);

    U2 = leakyGradient(y1*para.weightsl2+para.biasl2,alpha);

    U1 = leakyGradient(x*para.weightsl1+para.biasl1,alpha);


    J = para.weightsl8.'*...
        U7*para.weightsl7.'*...
        U6*para.weightsl6.'*...
        U5*para.weightsl5.'*...
        para.weightsl4.'*...
        U3*para.weightsl3.'*...
        U2*para.weightsl2.'*...
        U1*para.weightsl1.';
end

end

function [output] = leakyrelu(input,alpha)
    index = find(input>=0);
    index2 = find(input<0);
    output = zeros(size(input));
    output(index) = input(index);
    output(index2) = alpha.*(input(index2));
end

function [output] = leakyGradient(input,alpha)
    output = ones(size(input));
    output(input<0) = alpha;
    output = diag(output);
end
