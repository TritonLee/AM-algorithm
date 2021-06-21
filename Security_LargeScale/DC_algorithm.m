%% Outer iteratively optimize the sum rate 
% output the converged objective function value, variable Xi and Theta 

function [xi,theta_temp,fvalue] = DC_algorithm(ka,B,alp_k,P_t,tao,sigma_b,delta,M,N,K)

%% initialize the updated variables (xi, theta)
% directly analytical solution for upper bound value of Xi 
% clc;
% clear all;
% ka = 10;
% M = 10;    % the number of cluster in system 
% K = 10;    % the number of users 
% N = 100;   % the number of antennas in MISO systems
% B = 16;    % the feedback bits due to quantization and imperfect CSI 
% sigma_b = 10^(0/10);  % parameter related to sigmal_B
% tao = 2^(-B/(N-1))/(2*M);
% P_t = 10^(0/10);     % maximum total transmit power 
% alp_k = [10,7,6,6,4,3,3,1,1,1]
% delta = 0.5;

ITER_MAX = 50; % max number of outer iterations (usually converge within 5 iterations)
%% determine the upper bound of Xi with analytical solution
% A warm start of updating variable xi
Del = (1-delta)^(1/(M-1));
theta_tp = ones(K,1)*1/(M*K);
tep = zeros(1,K);
xi_analy_ub = zeros(1,K);
for k = 1:K
    tep(k) = sigma_b*M*2^(B/(N-1))/(P_t*(M-1)*alp_k(k));
    xi_analy_ub(k) = 2*P_t*alp_k(k)*(M-1)/sigma_b*lambertw(tep(k)*exp(tep(k))/Del) - 2*M*2^(B/(N-1));

    xi_s(k) = sum(theta_tp(k+1:K));
    xi_ini(k) = 1/(ka+xi_s(k));
end
xi0 = ones(1,K).*xi_ini;
Tol = 10^(-4);                   % search tolerance for iterative algorithm
%% iteratively optimize the based on AM framework
% update xi
[xi, fvalue_l] = update_Xi(xi0,theta_tp,ka,tao,alp_k,P_t,xi_analy_ub,sigma_b,M,K);
% update theta
[theta_te, fvalue_temp] = update_theta(theta_tp,xi,ka,sigma_b,P_t,alp_k,tao,M,K);

[theta_temp, fvalue_r] = CCP_Theta(xi,theta_te,ka,ITER_MAX,K,sigma_b,P_t,alp_k,tao,M);
Ra_left = fvalue_l;
Ra_right = fvalue_r;

% k = 1;
% while (abs(Ra_right-Ra_left)>Tol) && k<100
% 
%     [xi_t, fvalue_l_t] = update_Xi(xi,theta_temp,ka,tao,alp_k,P_t,xi_analy_ub,sigma_b,M,K);
%     [theta_r, fvalue_r_t] = CCP_Theta(xi_t,theta_temp,ka,ITER_MAX,K,sigma_b,P_t,alp_k,tao,M);
%     theta_temp = theta_r;
%     xi = xi_t;
%     Ra_left = fvalue_l_t;
%     Ra_right = fvalue_r_t;
%     k = k+1;
% 
% end

history.dual(1)=0;
iter = 0:ITER_MAX;
for i =1:length(iter)
    %stopping criterion
    Rate_opt = (Ra_right+Ra_left)/2;
    dual_obj = Rate_opt;
    history.dual(i+1)=dual_obj; 
    if i>=ITER_MAX
        if abs(history.dual(i+1)-history.dual(i))<Tol % line 15 of Algorithm 2
            break;
        end
    end
    [xi_t, fvalue_l_t] = update_Xi(xi,theta_temp,ka,tao,alp_k,P_t,xi_analy_ub,sigma_b,M,K);
    [theta_r, fvalue_r_t] = CCP_Theta(xi_t,theta_temp,ka,ITER_MAX,K,sigma_b,P_t,alp_k,tao,M);

    theta_temp = theta_r;
    xi = xi_t;
    Ra_left = fvalue_l_t;
    Ra_right = fvalue_r_t;
end
Rate_opt = (Ra_right+Ra_left)/2;
fvalue = Rate_opt;
end


% plot(iter,history.dual,'b-', 'LineWidth',2);hold on; 
% %% text of figure
% xlabel('Number of Outer Iterations');
% ylabel('Objective function value of (26)');
% legend('M = 2','M = 5','M = 10');
% grid on;
