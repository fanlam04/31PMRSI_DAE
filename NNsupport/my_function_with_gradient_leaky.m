function [f,g] = my_function_with_gradient_leaky(x,d,para,lambda,alpha)

    
[x_out,J] = my_forward_Op_leakyRelu(x,para,alpha);

f = norm(d-x)^2 + lambda*norm(x_out-x)^2;

if nargout > 1
    g = 2*(x-d) + 2*lambda*(x_out-x)*(J-eye(size(J,1)));
end
    
    
end
