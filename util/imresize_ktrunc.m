function img = imresize_ktrunc(img,newdim,doFilt,doNorm)
% IMRESIZE_KTRUNC creates a resized image through zero-padding or truncation in
% k-space. If doFilt is true then a Hamming window is applied.
% ---- INPUTS ----
% img    : original image
% newdim : new dimensions
% doFilt : flag for hamming window (defualts to false)
% doNorm : flag for normalization (defaults to true)
%
% ---- OUTPUTS ----
% img    : new image (old image over-written)
% Modified by Ronny @ 02/04/2018
%     for: change doFilt from global to for each dimension 
%     for: datatype
% ------------------------------------------------------------------------------

    if nargin < 4
        doNorm = true;
        if nargin < 3
            doFilt = false(1,ndims(img));
        end
    end
    if(length(doFilt)==1)
        doFilt = doFilt*ones(1,ndims(img));
    end
    
    dim = size(img);
    numdim = length(dim);
    for n=1:length(dim)
        if n<=length(newdim)
            img = reshape(img,dim(1),[]);
            img = resize1d(img,newdim(n),doFilt(n),doNorm);
            dim(1) = newdim(n);
            img = reshape(img,dim);
        end
        
        img = permute(img,[2:numdim,1]);
        dim = circshift(dim,[0,-1]);
    end


end


%% Subfunction: perform the operation along 1dim
function img = resize1d(img,newM,doFilt,doNorm)
    M = size(img,1);
    N = size(img,2);
    
    img = fftshift(img,1);
    img = fft(img,[],1)/sqrt(M);
    img = ifftshift(img,1);
    
    if newM > M
        tmp = zeros(newM,N,'like',img);
        tmp(floor((newM - M)/2) + 1 : floor((newM + M)/2),:) = img;
        img = tmp;
        clear tmp;
    else
        img = img(floor((M - newM)/2) + 1 : floor((M + newM)/2),:);
    end
    
    if doFilt
        img = img.*repmat(hamming(newM),[1,N]);
    end
    
    img = fftshift(img,1);
    img = ifft(img,[],1)*sqrt(newM);
    img = ifftshift(img,1);
    
    if doNorm
        img = img*sqrt(newM/M);
    end
    
end
