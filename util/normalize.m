function [out] = normalize(input,dim,mode)

switch mode
    case 'voxel-wise'
        if nargin>2
            out = input./max(abs(input),[],dim);
        else
            out = input./max(abs(input),[],2);
        end
    case 'global'
        out = input./max(abs(input(:)));
    otherwise
        error('Wrong mode for normalization.');
end
