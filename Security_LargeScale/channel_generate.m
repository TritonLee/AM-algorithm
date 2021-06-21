%% channel generation for beamforming design 
clc
clear;

 

M = 10;    % the number of cluster in system 
N_vec = [64 128 256 512 1024]; % the number of antennas in MISO systems
sigma_e = 10^(5/10);  % parameter related to sigmal_E
alp_e = 5;            % large scale path loss for eavesdropper 
P_t = 10^(10/10);     % maximum total transmit power 
epsilon = 0.1;

Monte = 2; % for testing with Monte Carlo Simulation
ka = zeros(1,M);
R = zeros(1,M);
for s = 1:length(N_vec)
    N = N_vec(s);
    BF = zeros(N,M,Monte);
    for m = 1: Monte
  
        W = gener_Channel(epsilon,P_t,alp_e,sigma_e,N,M)
        BF(1:N,:,m) = W

    end

end
% save file
save('AMchannel.mat','BF');



