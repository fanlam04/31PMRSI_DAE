function NNout3D = NN_initialization(input,para)

[L1,L2,L3,M]        = size(input);

input               = reshape(input,[L1*L2*L3,M]);
input               = contag(input);

sam                 = L1*L2*L3;
NNout               = zeros(L1*L2*L3,M*2);
for i=1:sam
    NNout(i,:) = my_forward_Op(input(i,:),para);
end

NNout3D        = reshape(decontag(NNout),[L1,L2,L3,M]);

end