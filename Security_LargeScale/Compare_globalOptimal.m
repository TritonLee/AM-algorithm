%% channel generation for beamforming design 
clc
clear all;

M = 2;    % the number of cluster in system 
K = 2;    % the number of users 
N = 100;   % the number of antennas in MISO systems
B = 2;    % the feedback bits due to quantization and imperfect CSI 
sigma_b = 10^(0/10);  % parameter related to sigmal_B
sigma_e = 10^(5/10);  % parameter related to sigmal_E
tao = 2^(-B/(N-1))/(2*M);
a_temp = ceil(rand(1,K)*10);
%alp_k = sort(a_temp,'descend');
alp_k = [10,3];
alp_e = 2;           % large scale path loss for eavesdropper 
P_t = 10^(10/10);     % maximum total transmit power 
epsilon = 0.5;
delta = 0.3;
eps = log(epsilon^(-1));

R = zeros(1,M);

load('channel.mat','W');
%load('AMchannel.mat','BF');

G = 1/sqrt(2)*(randn(M,N)+j*randn(M,N));  % random generated channel vector after channel estimation 
W = (G'*G)^(-1)*G';          % N times M
Del = (1-delta)^(1/(M-1));
tep = zeros(1,length(alp_k));
xi_analy_ub = zeros(1,length(alp_k));
for m = 1:M   % for each cluster 
    uni_wm = W(:,m)/norm(W(:,m))^2;   % normalized beamforming vector 
    S_m = uni_wm*uni_wm'; 
    Jtemp = zeros(N,N);
    for s = 1:M
        if s ~= m
            Jtemp = Jtemp + uni_wm*uni_wm';  
        end
    end
end

for a = 1:length(alp_k)
    tep(a) = sigma_b*M*2^(B/(N-1))/(P_t*(M-1)*alp_k(a));
    xi_analy_ub(a) = 2*P_t*alp_k(a)*(M-1)/sigma_b*lambertw(tep(a)*exp(tep(a))/Del) - 2*M*2^(B/(N-1));
    
end
xi_1 = [0:0.1:xi_analy_ub(1)];
xi_2 = [0:0.1:xi_analy_ub(2)];
theta_1 = [0.1:0.1:1];

for t = 1:length(theta_1)
    
    Teta_2(t) = 1-theta_1(t);
    for i = 1:length(xi_1)   
        Dm_temp = [0:0.1:10]
        for ii = 1:length(Dm_temp)
            Gama = (P_t*alp_e)/sigma_e*(theta_1(t).*S_m-Dm_temp(ii).*S_m*Teta_2(t)-Dm_temp(ii)/M.*Jtemp);
            V =  real(eigs(Gama,3));
            P_so(ii) = exp(-Dm_temp(ii)/V);
            if P_so(ii)<= epsilon
                R_m(i,t) = log2(1+xi_1(i)*theta_1(t)/(1+xi_1(i)*0))
                right(i,t,ii) = max(R_m(i,t)- Dm_temp(ii),0)
                B(i) = exp(xi_1(i)*sigma_b/(2*P_t*alp_k(1)))*(1+xi_1(i)*tao).^(M-1);
                Rate(i,ii,t) = right(i,t,ii)/B(i);
            else
                Rate(i,ii,t) = 0;
            end
        end
    end
    Sum_ra1(t) = max(max(Rate));
    
    for j = 1:length(xi_2)
        Dm_temp2 = [0:0.1:10]
        for ii = 1:length(Dm_temp2)

            Gama2 = (P_t*alp_e)/sigma_e*(Teta_2(t).*S_m-Dm_temp2(ii).*S_m*theta_1(t)-Dm_temp2(ii)/M.*Jtemp);
            V2 =  real(eigs(Gama2,3));
            P_so2(ii) = exp(-Dm_temp2(ii)/V2);

            if P_so2(ii)<= epsilon

                R_m2(j,t) = log2(1+xi_2(j)*Teta_2(t)/(1+xi_2(i)*theta_1(t)))
                right(j,t,ii) = max(R_m2(j,t)- Dm_temp2(ii),0)
                B(j) = exp(xi_2(j)*sigma_b/(2*P_t*alp_k(2)))*(1+xi_2(j)*tao).^(M-1);
                Rate2(j,ii,t) = right(j,t,ii)/B(j);
            else
                Rate2(j,ii,t) = 0;
            end     

        end 
    end
    Sum_ra2(t) = max(max(Rate2));
    
    to_sum_ra(t) = Sum_ra1(t) + Sum_ra2(t);
end
Obj_rate = max(to_sum_ra)





plot(epsilon,Ac_Pso,'r', 'LineWidth',2);hold on;
%plot(zk,q,'r', 'LineWidth',2);hold on;
xlabel('\epsilon');
ylabel({'Actual SOP, $p_{so}^{m,k}$'},'Interpreter','latex');
grid on;

