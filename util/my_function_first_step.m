function x = my_function_first_step(x0,B,S,Y,para,lambda1,mu)

[L1,L2,L3,L4] = size(x0);


if L4 == 1 && L3 ~= 1
    x0 = reshape(x0,[L1*L2,L3]);
    B  = reshape(B,[L1*L2,L3]);
    S  = reshape(S,[L1*L2,L3]);
    Y  = reshape(Y,[L1*L2,L3]);
elseif L4 ~= 1
    x0 = reshape(x0,[L1*L2*L3,L4]);
    B  = reshape(B,[L1*L2*L3,L4]);
    S  = reshape(S,[L1*L2*L3,L4]);
    Y  = reshape(Y,[L1*L2*L3,L4]);
end


options  =  optimoptions('fminunc','Algorithm','quasi-newton','HessUpdate','bfgs','SpecifyObjectiveGradient',true,'MaxIterations',50,'StepTolerance',1e-5,'Display','none');

x0       =  contag(x0);
x        =  zeros(size(x0));


parfor i = 1:size(x0,1)
    x(i,:) = fminunc(@(x)my_function_with_gradient_first_step(x,B(i,:),S(i,:),Y(i,:),para,lambda1,mu),x0(i,:),options);
end

x        =  decontag(x);

if L4 == 1 && L3 ~= 1
    x = reshape(x,[L1,L2,1,L3]);
elseif L4 ~= 1
    x = reshape(x,[L1,L2,L3,L4]);
end
  
end


function [f,g] = my_function_with_gradient_first_step(x,B,S,Y,para,lambda1,mu)
% the input x is in a matrix (L1,L2,M) form
 

[x_out,J]    =     my_forward_Op(x,para);% retrun the jacobian matrix and gradient


%function value
f            =     lambda1*norm(x_out - x)^2 + (mu/2)*norm(B.*decontag(x) - S + Y/mu)^2;

%function gradient
if nargout > 1
    g = gradient1(x,x_out,J,lambda1) + contag(mu*(B.*(decontag(x)) - S + Y/mu).*conj(B));
end
    
    
end


function g = gradient1(x_in,x_out,J,lambda1)


g          = 2*lambda1*(x_out-x_in)*(J-eye(size(J,1)));


end