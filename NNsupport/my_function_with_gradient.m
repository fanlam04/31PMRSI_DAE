function [f,g] = my_function_with_gradient(x,d,para,lambda)


g = zeros(size(x,1),size(x,2));
x_out = zeros(1,size(x,2));

parfor i = 1:size(x,1)
    
    [x_out(i,:),J] = my_forward_Op(x(i,:),para);

    %if nargout > 1
        g(i,:) = 2*(x(i,:)-d(i,:)) + 2*lambda*(x_out(i,:)-x(i,:))*(J-eye(size(J,1)));
    %end
    
    
end


f = norm(d-x,'fro')^2 + lambda*norm(x_out-x,'fro')^2;