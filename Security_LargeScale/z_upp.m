function [ value ] = z_upp(x,tao,delta,sigma_b,P_t,K)
% Bisect get minmax value of Transmission Power

%   Detailed explanation goes here
    left = exp(-x*sigma_b/P_t)/((1+x*tao^(-1))^(K-1));
    value = left - (1-delta);
end