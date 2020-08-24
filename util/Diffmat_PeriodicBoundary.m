function [D,Dp] = Diffmat_PeriodicBoundary(N1,N2,N3)
% Generates the sparse three-dimensional finite difference (first-order neighborhoods) matrix D 
% for an image of dimensions N1xN2xN3 (rows x columns x slices).  The optional output
% output argument Dp is the transpose of D.  Also works for 2D images if N3 is 1.
% ------ Fan Lam 15/08/2012 ------ %
% Change-log
% (2015-09-22, Bryan Clifford) Swap out the missing "GetDiffMat" with Diffmat_periodicBoundary2D
%-------------------------------------------------------------------------------

if (not(isreal(N1)&&(N1>0)&&not(N1-floor(N1))&&isreal(N2)&&(N2>0)&&not(N2-floor(N2))))
    error('Inputs must be real positive integers');
end
if ((N1==1)&&(N2==1)&&(N3==1))
    error('Finite difference matrix can''t be generated for a single-pixel image');
end

D1 = [];
D2 = [];
D3 = [];

if (N1 > 1)&&(N2>1)&&(N3>1)    
    e = ones(N1,1);
    if (numel(e)>2)
        T = spdiags([e,-e],[0,1],N1,N1);
        T(N1,1)=-1;
        E = speye(N2);
        E2 = speye(N3);
        D1 = kron(E2,kron(E,T)); % column wise finite difference
    end
    e = ones(N2,1);
    if (numel(e)>2)
        T = spdiags([e,-e],[0,1],N2,N2);
        T(N2,1)=-1;
        E = speye(N1);
        E2 = speye(N3);
        D2 =  kron(E2,kron(T,E)); % row wise finite difference
    end
    e = ones(N3,1);
    if (numel(e)>2)
        T = spdiags([e,-e],[0,1],N3,N3); % slice wise finite difference
        T(N3,1)=-1;
        E = speye(N1);
        E2 = speye(N2);
        D3 =  kron(T,kron(E2,E));
    end
    D = [D1;D2;D3];
elseif (N1 > 1)&&(N2>1)&&(N3 == 1)    % for 2D Image
    [Dv, Dh] = Diffmat_PeriodicBoundary2D(N1,N2,1);
    D = [Dv; Dh];
else
   error('singleton dimensions not supported');
end

if (nargout > 1)
    Dp = D';
end
