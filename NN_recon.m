%------------------NN reconsturction with ADMM algorithm------------------%
%
%   Description: Perform NN spatiospectral reconstruction of the provided
%                noiy data
%   
%   Input: 1. noisy x-t space data
%          2. B0 map
%          3. anatomical images
%
%   Output: reconstructed x-t data
%
%
%   Parameters: 1. dt 
%               2. TE 
%               3. unippm
%               4. recon_dims:   reconstruction dimension.
%                                Default is the same as the data dimension
%
%               5. lambda1:      NN regularization parameter.
%                                Default: 10
%
%               6. lambda2:      spatial regularization parameter.
%                                Default: 0.1               
%
%               7. mu:           penalty parameter in the ADMM algorithm.
%                                Default: 3
%
%               8. cgtol:        tolerance of the linear solver. 
%                                Default: 1e-6
%
%               9. max_iter:     max iteration of the linear solver. 
%                                Default:50
%
%               10. err_thresh:  threshold for exiting the ADMM algorithm.
%                                Default:1e-3
%
%               11. coef, alpha: coefficient that modify the penalty
%                                parameter mu for different iterations for faster
%                                convergence.
%                                Default: 1.06 and 1
%
%
%   Author: Yahang Li @ UOFI, Lam's Group
%   Date: 2020/08/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear;
clc;

% add path to the utility codes
addpath('./NNsupport');
addpath('./util');

% load NN parameters; Network with middle layer dimension 24 is used here
load('./NNsupport/para_nn24.mat');

% load data
load('./data.mat');
load('./B0.mat');
load('./anat.mat');
xt_data = permute(comb,[2 3 4 1]);

% time grids
dt          = 1/5000;
unippm      = 161.5461;
TE          = 1.36e-3;


% set up ppm range of the metabolites
ppm_PCr     = [-0.30, 0.40];
ppm_ATPa    = [-7.80,-7.00];
ppm_ATPg    = [-3.0, -2.1];
ppm_Pi      = [4.80, 5.20];

%% Data preparation
xt_noisy    = (xt_data(:,:,6:8,:));
ktData      = F3_x2k(xt_noisy);

dims_data   = size(xt_noisy);
L1          = dims_data(1);
L2          = dims_data(2);
L3          = dims_data(3);
M           = dims_data(4);
recon_dims  = [L1,L2,L3,M]; % reconstruction dimension (can be different from data dimension)
xt_noisy    = imresize_ktrunc(xt_noisy,recon_dims(1:3));
Nt          = M;
t_vec       = TE+[0:M-1]*dt;
fvec        = genfvec(Nt, 1/dt);
ppm_vec     = fvec./unippm;

B0Map       = FreqShift; 
fm_hz_31p   = B0Map;
fm_hz_orig  = fm_hz_31p;
Iref        = normRange(abs(anat));
mask_orig   = Iref>0.001;


kMask       = abs(ktData(:,:,:,1)) > 1e-12;

mask        = logical(imresize3d(logical(mask_orig),recon_dims(1:3)));
fm_hz_31p   = real(imresize_ktrunc(real(imresize_ktrunc(fm_hz_31p,size(fm_hz_31p),1)),recon_dims(1:3))).*mask;

    
% edge info
para_align2 = [0,0,0];
edge_ratio  = 20;
edgemap     = process_gre3D(Iref,para_align2,recon_dims(1),recon_dims(2),recon_dims(3),edge_ratio,mask,102);
w           = edgemap(:);
[D,~]       = Diffmat_PeriodicBoundary(recon_dims(1),recon_dims(2),recon_dims(3));
D           = spdiags(sqrt(w),0,size(D,1),size(D,1))*D; % appy weightings to D

% B0 info
B0MapD1     = fm_hz_31p; % 
B           = reshape(exp(1i*2*pi*B0MapD1(:)*t_vec), recon_dims);
B_orig      = reshape(exp(1i*2*pi*fm_hz_orig(:)*t_vec), [size(B0Map),length(t_vec)]);

% B0 correction
xt_noisy    = imresize_ktrunc(conj(B).*xt_noisy,recon_dims(1:3)); % ideal x0 should not contain B0 effects

% normalization
xt_noisy    = reshape(xt_noisy,[L1*L2*L3,M]);
xt_noisy    = decontag(normalize(contag(xt_noisy),[],'global'));
xt_noisy    = reshape(xt_noisy,[L1,L2,L3,M]);

% linear subspace
% Vi          = Vt_all(1:M,1:48);
% [Ut,St,Vt]  = svd(Vi,'econ');
% Vi          = Ut(:,1:48).';
% xt_tmp      = reshape(xt_noisy,[L1*L2*L3,M]); 
% [Un,Sn,Vn]  = svd(xt_tmp,'econ');
% xt_tmp      = reshape(Un(:,1:64)*Sn(1:64,1:64)*Vn(:,1:64)',[L1,L2,L3,M]);
% xt_noisy    = reshape((xt_noisy*Vi')*inv(Vi*Vi')*Vi,[L1,L2,L3,M]); % Projection onto the learned subspace

%% Reconstruction
% setting up initializatoin
delete(gcp('nocreate'));

kMask          = repmat(kMask,[1,1,1,M]);
para           = para_nn;

% ALG parameter setup
lambda1        = 10;
lambda2        = 0.1;
mu             = 3;
iter           = 15;
cgtol          = 1e-6;
max_iter       = 50;

err_thresh     = 1e-3;
coef           = 1.06;
alpha          = 1;


% initialization
x0             = xt_noisy; % ideal x0 should not contain B0 effects
dkt            = F3_x2k(xt_noisy); 

B              = ones(size(x0));
S_last         = B.*x0;
Y_last         = zeros(L1,L2,L3,M);

err            = zeros(iter,1);

x0             = NN_initialization(x0,para_nn); % initialize x0 value
x0             = x0.*mask;

S_last         = my_function_second_step(S_last,kMask,x0,D,B,Y_last,dkt,lambda2,mu,cgtol,max_iter);
X_last         = x0;
for i = 1:iter

    % ALG begins
    disp(['iter:', num2str(i), '/', num2str(iter)]);

    X          = my_function_first_step(X_last,B,S_last,Y_last,para,lambda1,mu);

    S          = my_function_second_step(S_last,kMask,X,D,B,Y_last,dkt,lambda2,mu,cgtol,max_iter);

    Y          = Y_last + mu*alpha*(B.*X - S);

    err(i)     = norm(X_last(:) - X(:))/norm(X_last(:));

    if err(i) < err_thresh
        break;
    end


    %update
    X_last     = X;
    S_last     = S;
    Y_last     = Y;

    alpha      = alpha*coef;


end


% exam results
xt_denoise  = reshape((X),L1,L2,L3,[]);

sxf_D2      = F3_t2f(xt_noisy);
sxf_denoise = F3_t2f(xt_denoise);
locs_PCr    = (ppm_vec>=ppm_PCr(1)) & (ppm_vec<=ppm_PCr(2));
locs_ATPa   = (ppm_vec>=ppm_ATPa(1)) & (ppm_vec<=ppm_ATPa(2));
locs_ATPg   = (ppm_vec>=ppm_ATPg(1)) & (ppm_vec<=ppm_ATPg(2));
locs_Pi     = (ppm_vec>=ppm_Pi(1)) & (ppm_vec<=ppm_Pi(2));

img_PCr_D2  = sqrt(sum(abs(sxf_D2(:,:,:,locs_PCr)).^2,4));
img_PCr_de  = sqrt(sum(abs(sxf_denoise(:,:,:,locs_PCr)).^2,4));

img_ATPa_D2 = sqrt(sum(abs(sxf_D2(:,:,:,locs_ATPa)).^2,4));
img_ATPa_de = sqrt(sum(abs(sxf_denoise(:,:,:,locs_ATPa)).^2,4));

img_ATPg_D2 = sqrt(sum(abs(sxf_D2(:,:,:,locs_ATPg)).^2,4));
img_ATPg_de = sqrt(sum(abs(sxf_denoise(:,:,:,locs_ATPg)).^2,4));

img_Pi_D2   = sqrt(sum(abs(sxf_D2(:,:,:,locs_Pi)).^2,4));
img_Pi_de   = sqrt(sum(abs(sxf_denoise(:,:,:,locs_Pi)).^2,4));  

figure,
subplot(3,2,1), montagesc(normRange(img_PCr_de)),axis image off,title(sprintf('PCr denoise in lambda1=%.2f, lambda2=%.2f',lambda1,lambda2));
subplot(3,2,3), montagesc(normRange(img_ATPa_de)),axis image off,title(sprintf('ATPa denoise in lambda1=%.2f, lambda2=%.2f',lambda1,lambda2));
subplot(3,2,5), montagesc(normRange(img_Pi_de)),axis image off,title(sprintf('Pi denoise in lambda1=%.2f, lambda2=%.2f',lambda1,lambda2));
subplot(3,2,2), montagesc(normRange(img_PCr_D2)),axis image off,title(sprintf('PCr noisy in lambda1=%.2f, lambda2=%.2f',lambda1,lambda2));
subplot(3,2,4), montagesc(normRange(img_ATPa_D2)),axis image off,title(sprintf('ATPa noisy in lambda1=%.2f, lambda2=%.2f',lambda1,lambda2));
subplot(3,2,6), montagesc(normRange(img_Pi_D2)),axis image off,title(sprintf('Pi noisy in lambda1=%.2f, lambda2=%.2f',lambda1,lambda2));


figure, hold on
plot(ppm_vec,real(fftshift(fft(squeeze(xt_noisy(17,6,2,:)))))),
set(gca,'XDir','Reverse','XLim',[-15,15],'YLim',[-10,15],'Box','off'),
plot(ppm_vec,squeeze(real(fftshift(fft(X(17,6,2,:)))))),
set(gca,'XDir','Reverse','XLim',[-15,15],'YLim',[-10,15],'Box','off'),
title(sprintf('ALG denoised in lambda1=%.2f, lambda2=%.2f',lambda1,lambda2)),box off,

% save data
name  = sprintf('ReALG_EXdenoise_data2_NNinit_alpha%.2f_lambda1_%d_lambda2_%.2f.mat',coef,lambda1,lambda2);
save(name,'X','err','t_end');
    


% sendgmail('researchShuai@gmail.com','Finished!','Process finished!');


