function [ value ] = xi_upp(x,tao,delta,alp_k,sigma_b,P_t,M)
% Bisect get minmax value of Transmission Power

%   Detailed explanation goes here
    left = exp(-x*sigma_b/(2*P_t*alp_k))/((1+x*tao)^(M-1));
    value = left - (1-delta);
end