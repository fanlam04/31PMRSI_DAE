function edge_weights = process_gre3D(Iref,para_align,gridy,gridx,gridz,ratio,Mask_ref,fignum,flag_disp)
% This function is used calculate weighted-averaged edge weights from a single 3D reference image
% Input:  Iref      : multiple reference images, each of which already normalized
%         para_align: alignment matrix
%         gridy     : grid size along y
%         gridx     : grid size along x
%         opt       :   related parameters
% Output: edge_weights, weights along all directions (gridy*gridx*#directions)
% Created by Fan Lam@Univ of Illinois, 01/29/2014

    % --- prepare parameters --- %
    if isempty(fignum)
        fignum          = 100;
    end
    
    if ~exist('flag_disp','var')
        flag_disp       = false;
    end
    
    Mask_ref            = imresize3d(Mask_ref,[gridy,gridx,gridz])>0.1;
    [N1,N2,N3]          = size(Iref);
    N                   = N1*N2*N3;
    [D,~]               = Diffmat_PeriodicBoundary(N1,N2,N3); % finite difference operators with periodic boundaries
    Y_shift             = para_align(1);
    X_shift             = para_align(2);
    Z_shift             = para_align(3);
    Mask_ref            = circshift(Mask_ref,[Y_shift,X_shift,Z_shift]);

    % --- extract edges --- %
    if N3 == 1
        Ndir            = 2;
    else
        Ndir            = 3;
    end
    edge_weights        = zeros(gridy,gridx,gridz,Ndir);
    
    edges               = abs(D*reshape(Iref,[],1));          % extract the edges
    alpha               = max(edges(:))/ratio;
    w                   = 1.*(edges < alpha) + (alpha./edges).*(edges>=alpha); % a simpler but less general form
    w(isnan(w))         = 1;
    w(isinf(w))         = 1;                                  % 3N*1
    
    % --- alignment --- %
    for nr  = 1:Ndir
        wi              = reshape(w((nr-1)*N+1:nr*N),N1,N2,N3);
        % alignment
        wi              = imresize3d(wi,[gridy,gridx,gridz]);
%         wi              = abs(wi./max(abs(wi(:))));
        wi              = circshift(wi,[Y_shift,X_shift,Z_shift]);
        wi              = wi.*Mask_ref;             % apply mask for relevant regions only

        wi(isnan(wi))   = 1;             % irrelevant region should have maximum weights
        wi(wi<=0)       = 1; 

        edge_weights(:,:,:,nr) = wi;
    end % end-of-loop for directions
        

    

    edgemap         = sqrt(sum(edge_weights.^2,4)); % for visualizing the weights map
    edgemap         = edgemap./max(edgemap(:));
    
    if flag_disp
        figure(fignum), montagesc(edgemap), colormap gray
    end
%     saveFigEps(gcf,'edge_weights3D_dircombined');
%     imagesc(edgemap(:,:,floor(end/2)+1)),colormap gray,axis image off%,title('aligned composite edge weightes of center slice')
    
end