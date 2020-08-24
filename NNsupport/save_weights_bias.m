%% save NN parameters
close all;
clear;
clc;

load('./GZ_0ordPhase_gaussian_10HzB0_weights2000_250_100_24_dropout05.mat');
load('./GZ_0ordPhase_gaussian_10HzB0_bias2000_250_100_24_dropout05.mat');

para_nn.weightsl1 = double(weightsl1);
para_nn.weightsl2 = double(weightsl2);
para_nn.weightsl3 = double(weightsl3);
para_nn.weightsl4 = double(weightsl4);
para_nn.weightsl5 = double(weightsl5);
para_nn.weightsl6 = double(weightsl6);
para_nn.weightsl7 = double(weightsl7);
para_nn.weightsl8 = double(weightsl8);

para_nn.biasl1 = double(biasl1);
para_nn.biasl2 = double(biasl2);
para_nn.biasl3 = double(biasl3);
para_nn.biasl4 = double(biasl4);
para_nn.biasl5 = double(biasl5);
para_nn.biasl6 = double(biasl6);
para_nn.biasl7 = double(biasl7);
para_nn.biasl8 = double(biasl8);


save('./para_nn24.mat','para_nn');