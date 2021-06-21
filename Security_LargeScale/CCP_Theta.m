function [Theta,fvalue] = CCP_Theta(xi,theta_ini,ka,ITER_MAX,K,sigma_b,P_t,alp_k,tao,M)

% f1 = zeros(1,K);
% f2 = zeros(1,K);
grad = zeros(1,K);
iter = 0:ITER_MAX;
cvx_prev = inf;
tol=1e-4;       % target accuracy

% for k = 1:K
%     xi_s(k) = sum(theta0(k+1:K));
%     xi_ini(k) = 1/(ka+xi_s(k));
% end
% xi = ones(1,K).*xi_ini;
% theta0 = ones(K,1)*1/(M*K);
theta0 = theta_ini;
for i = 1:length(iter)

    cvx_begin 
    cvx_solver sdpt3
    cvx_precision best
    %cvx_expert true
    variable theta(K,1) nonnegative
    variable t(K,1) 

    maximize sum(t)
    subject to   
    for k=1:K
        t(k)>=0;
        theta(k)>=0;
    end
    %% derive the gradient approx. 
    for k = 1:K
        temp_gra = zeros(1,K);
        for jj = 1:k-1
            temp_gra(jj) = xi(k)/(1+xi(k)*sum(theta0(1:k-1))); 
        end
        for jj = k:K
            temp_gra(jj) = 0;
        end
        grad(k) = sum(temp_gra);
      
        theta_temp = 0;
        for ii = 1:K
            if ii ~= k
                theta_temp = theta_temp + theta(ii);
            end 
        end

        f1 = (log(1+xi(k)*sum(theta(1:k)))+log(ka+theta_temp))/log(2);
        f2 = log(1+xi(k)*sum(theta0(1:k-1)))/log(2);         
        log(1+xi(k)*sum(theta(1:k)))+log(ka+theta_temp) - log(1+xi(k)*sum(theta0(1:k-1))) - grad(k)*(theta(k)-theta0(k))>=log(2)*(log(ka+1/M)+t(k));
        %Rate_temp >= t(k)      
    end

    sum(theta) == 1/M;
    %sum(f1) - sum(f2)-grad*(theta-theta0)>=t   
    cvx_end
    
    theta0=theta; % update variables for next iteration
    
    if (norm(cvx_prev-sum(t))<tol) % stopping criterion
        break;
    end
    cvx_prev = sum(t);
    
    %Rate_history(i+1)=real(sum(t));
end
Theta = theta;
Obj_D2 = Compute_D2(theta,xi,K,ka,sigma_b,P_t,alp_k,tao,M);
fvalue = Obj_D2;

end