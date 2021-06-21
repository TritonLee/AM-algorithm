%% updating variable Theta based on projected gradient method 

%% considering the accelaration method for first-order method 
function [theta, fvalue] = update_theta(theta_tp,xi,ka,sigma_b,P_t,alp_k,tao,M,K)

% clc;
% clear;
% ka = 10;
% M = 10;    % the number of cluster in system 
% K = 10;    % the number of users 
% N = 100;   % the number of antennas in MISO systems
% B = 16;    % the feedback bits due to quantization and imperfect CSI 
% sigma_b = 10^(0/10);  % parameter related to sigmal_B
% tao = 2^(-B/(N-1))/(2*M);
% P_t = 10^(10/10);     % maximum total transmit power 
% alp_k = [10,7,6,6,4,3,3,1,1,1];
% theta_tp = ones(K,1)*1/(M*K);
% for k = 1:K
%     xi_s(k) = sum(theta_tp(k+1:K));
%     xi_ini(k) = 1/(ka+xi_s(k));
% end
% xi = ones(1,K).*xi_ini;

%% update theta based on projected-gradient method 
ITER_MAX=300;   % max number of outer iterations (usually converge within 30 iterations)
chi=1e-4;       % target accuracy
L = 1e-4;         % constant step size chosen by Armijo rule 
c = 1;
history.dual(1)=0;
theta = theta_tp;
theta1 = zeros(1,K);
iter = 0:ITER_MAX;
for i = 1:length(iter)
    
    c0 = c;
    c = (1+sqrt(1+4*c0^2))./2;                    % tuned parameter 
    rho = real(theta+(c0-1)./c*(theta-theta_tp)); % accelerate procedure
    
    % update theta with gradient step 
    gradient = zeros(1,K);
    for k = 1:K
%         temp_tep = zeros(1,k-1);
%         for ii = 1:k-1
%            temp_tep(ii) = theta_tp(ii)/(ka+1/M)/(ka+1/M-theta_tp(ii)); 
%         end
%         gradient(k) = 1/log(2)*(xi(k)/(1+xi(k)*sum(theta_tp(1:k)))-1/(ka+1/M)+sum(temp_tep));
%       
        
%% derive the gradient approx. 
        theta_temp = 0;
        for ii = 1:K
            if ii ~= k
                theta_temp = theta_temp + theta_tp(ii);
            end
        end
        
        temp_gra = zeros(1,K);
        for jj = 1:k-1
            temp_gra(jj) = xi(k)/(1+xi(k)*sum(theta_tp(1:k)))-xi(k)/(1+xi(k)*sum(theta_tp(1:k-1)))+1/(ka+theta_temp); 
        end
        temp_gra(k) = xi(k)/(1+xi(k)*sum(theta_tp(1:k)));
        for jj = k+1:K
            temp_gra(jj) = 1/(ka+theta_temp);
        end
        gradient(k) = 1/log(2)*sum(temp_gra);
    end
       
 
    L = 1e-3*1/max(gradient);
    theta1 = rho' + L*gradient;

    % update with projected step 
    mu0 = 2/K*(sum(theta1)-1/M);
    for k = 1:K 
        theta(k) = max(real(theta1(k) - mu0/2),0);
    end
    theta_tp = theta; % update theta0
    
    % compute objective value 
    Obj_D2 = Compute_D2(theta,xi,K,ka,sigma_b,P_t,alp_k,tao,M);
    dual_obj = Obj_D2;
    history.dual(i+1)=dual_obj; 
    % stopping criterion
    if i>=ITER_MAX
        if abs(history.dual(i+1)-history.dual(i))<chi % line 15 of Algorithm 2
            break;
        end
    end
end
fvalue = dual_obj;
end


%% text of figure
% plot(iter,history.dual,'b-', 'LineWidth',2);hold on;
% 
% xlabel('Number of PG Iterations');
% ylabel('Objective Function value of problem (26)');
% ylabel({'$f(\theta_m)$ in D2'},'Interpreter','latex');
% legend('M = 2','M = 5','M = 10');
% grid on;


