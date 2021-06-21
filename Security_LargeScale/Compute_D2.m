function [fvalue] = Compute_D2(theta,xi,K,ka,sigma_b,P_t,alp_k,tao,M)
for k = 1:K
    theta_temp = 0;
    for ii = 1:K
        if ii ~= k
            theta_temp = theta_temp + theta(ii);
        end
    end
    thta_tem = sum(theta(1:k-1));
    R_m(k) = log2(1+xi(k)*theta(k)/(1+xi(k)*thta_tem));
    D_m(k) = log2(1+theta(k)/(ka+theta_temp));
    right(k) = max(R_m(k)-D_m(k),0);  % expression of A_k
    left(k) = exp(-xi(k)*sigma_b/(2*P_t*alp_k(k)))/(1+xi(k)*tao)^(M-1);
    rate_temp(k)= left(k)*right(k);
end
fvalue = sum(rate_temp);
end