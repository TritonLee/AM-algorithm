%% updating variable Xi based on FP method with a local optimal solution 

function [xi, fvalue] = update_Xi(xi0,theta_tp,ka,tao,alp_k,P_t,xi_analy_ub,sigma_b,M,K)

% clc;
% clear;
% ka = 100
% delta = 0.5;
% M = 10;    % the number of cluster in system 
% K = 10;    % the number of users 
% N = 100;   % the number of antennas in MISO systems
% B = 16;    % the feedback bits due to quantization and imperfect CSI 
% sigma_b = 10^(0/10);  % parameter related to sigmal_B
% tao = 2^(-B/(N-1))/(2*M);
% P_t = 10^(10/10);     % maximum total transmit power 
% alp_k = [10,7,6,6,4,3,3,1,1,1]
% Del = (1-delta)^(1/(M-1));
% theta_tp = ones(K,1)*1/(M*K);
% for k = 1:K
%     xi_s(k) = sum(theta_tp(k+1:K));
%     xi_ini(k) = 1/(ka+xi_s(k));
%     tep(k) = sigma_b*M*2^(B/(N-1))/(P_t*(M-1)*alp_k(k));
%     xi_analy_ub(k) = 2*P_t*alp_k(k)*(M-1)/sigma_b*lambertw(tep(k)*exp(tep(k))/Del) - 2*M*2^(B/(N-1));
% end
% %xi0 = ones(1,K).*xi_ini;
% xi0 = ones(1,K);


xi = zeros(1,K);
rate = zeros(1,K);
right = zeros(1,K);  
D_m = zeros(1,K);  
R_m = zeros(1,K);  
opt_y = zeros(1,K);  
thta_tem = zeros(1,K);
%% determine the left hand side of updated Xi 
for k = 1:K
    xi_s(k) = sum(theta_tp(k+1:K));
    xi_ini(k) = 1/(ka+xi_s(k));
end
xi_left = ones(1,K).*xi_ini;

for k = 1:K
    theta_temp = 0;
    for ii = 1:K
        if ii ~= k
            theta_temp = theta_temp + theta_tp(ii);
        end
    end
    thta_tem(k) = sum(theta_tp(1:k-1));
    D_m(k) = log2(1+theta_tp(k)/(ka+theta_temp));   
    %% update xi based on FP framework 
    ITER_MAX=20; % max number of iterations
    chi=1e-4; % target accuracy
    history.dual(1)=0;
    iter = 0:ITER_MAX;
    for i=1:length(iter)

        R_m(k) = log2(1+xi0(k)*theta_tp(k)/(1+xi0(k)*thta_tem(k)));
        right(k) = max(R_m(k)-D_m(k),0);  % expression of A_k

        A(k) = right(k);
        B(k) = exp(xi0(k)*sigma_b/(2*P_t*alp_k(k)))*(1+xi0(k)*tao).^(M-1);
        opt_y(k) = sqrt(A(k))/B(k);  % the update y_k

        g_gm = @(xi)(2*opt_y(k)*sqrt((log2(1+xi*theta_tp(k)/(1+xi*thta_tem(k)))-D_m(k)))-(opt_y(k))^2*exp(xi*sigma_b/(2*P_t*alp_k(k)))*(1+xi*tao)^(M-1));
        [xi_ud,fval] = goldmax(g_gm,xi_left(k),xi_analy_ub(k),1e-8); % update xi 
        xi0(k)=xi_ud; % update xi0

        dual_obj=real(fval);
        history.dual(i+1)=dual_obj; 

        % stopping criterion
        if i>=ITER_MAX
            if abs(history.dual(i+1)-history.dual(i))<chi 
                break;
            end
        end  
    end
    xi(k) = xi0(k);
    %rate(k) = fval;
    R_m(k) = log2(1+xi(k)*theta_tp(k)/(1+xi(k)*thta_tem(k)));
    right(k) = max(R_m(k)-D_m(k),0);  % expression of A_k
    A(k) = right(k);
    B(k) = exp(xi(k)*sigma_b/(2*P_t*alp_k(k)))*(1+xi(k)*tao).^(M-1);
    rate(k) = A(k)/B(k);
    
end
SumRate = sum(rate);
fvalue = SumRate;
end


% plot(iter,history.dual,'b^-', 'LineWidth',2);hold on;
% 
% %% text of figure
% xlabel('Number of FP Iterations');
% ylabel({'$g(\Xi_{m,k},y_k)$ in Q2'},'Interpreter','latex');
% legend('M=2','M=5','M=10');
% grid on;


