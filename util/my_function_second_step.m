function S_out = my_function_second_step(S0,kMask,X,D,B,Y,dkt,lambda2,mu,cgtol,max_iter)
% to faster the ALG, fftshift the kMask to sample
[L1,L2,L3,M] = size(S0);

dxt       = F3_k2x(kMask.*dkt);
part2     = (mu/2)*(B.*X + Y./mu);
b         = dxt + part2;


S_out     = pcg(@(S)CallAx(S,kMask,B,D,lambda2,mu,L1,L2,L3,M),b(:),cgtol,max_iter, [], [], S0(:));



if M == 1 && L3 ~= 1
    
    S_out = reshape(S_out,[L1,L2,L3]);
    
elseif M ~= 1
    
    S_out = reshape(S_out,[L1,L2,L3,M]);

end

end

function y = CallAx(S,kMask,B,D,lambda2,mu,L1,L2,L3,M)

S        = reshape(S,[L1,L2,L3,M]);

Skt      = F3_x2k(S);
Sxt      = F3_k2x(kMask.*Skt);

BHS      = conj(B).*S;
DHDBHS   = D'*D*reshape(BHS,[L1*L2*L3,M]);

part2    = lambda2*(B.*(reshape(DHDBHS,[L1,L2,L3,M]))) + (mu/2).*S + (lambda2/50).*S;

y        = Sxt + part2;

y        = y(:);
end