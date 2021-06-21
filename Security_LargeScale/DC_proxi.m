%% solving optimization problem D1 with d.c structure
%% updating variable Theta based on DC programming


clc;
clear;
ka = 1;
M = 10;    % the number of cluster in system 
K = [2:2:10];    % the number of users 
N = 100;   % the number of antennas in MISO systems
B = 4;    % the feedback bits due to quantization and imperfect CSI 
sigma_b = 10^(0/10);  % parameter related to sigmal_B
tao = 2^(-B/(N-1))/(2*M);
P_t = 10^(10/10);     % maximum total transmit power 

%a_temp = ceil(rand(1,K)*10);
%alp_k = sort(a_temp,'descend');
alp_k = [10,7,6,6,4,3,3,1,1,1];
delta = 0.5;
t_CCP = zeros(1,length(K));
t_PG = zeros(1,length(K));
R_am = zeros(1,length(K));
R_dc = zeros(1,length(K));
for s = 1:length(K)
    
    %% solving using the first-order method
    tic;
    [xi,theta_temp,R_am(s)] = AM(ka,B,alp_k,P_t,tao,sigma_b,delta,M,N,K(s))
    t_PG(s) = toc;

    
    %% solving use the DC programming method 
    tic;
    %[xi1,theta_temp1,R_dc(s)] = DC_algorithm(ka,B,alp_k,P_t,tao,sigma_b,delta,M,N,K(s))
    
    %t_CCP(s) = toc;

end
No_u = M.*K;   % total number of users in transmission networks (No. of clusters * No. of users in each cluster)
%Obj_D2 = Compute_D2(theta,xi,K,ka,sigma_b,P_t,alp_k,tao,M);
%fvalue = Obj_D2;
semilogy(No_u,t_PG,'ms-', 'LineWidth',2); hold on;
%semilogy(No_u,t_CCP,'gs--', 'LineWidth',2); hold on;
grid on;
xlabel('Number of Users');
ylabel('Average Computation time (s)');
%ylabel('Objective function value of (29)');
%legend('Algorithm 3','DC programming');



% figure
% plot(xi,'b-', 'LineWidth',2);hold on;
% plot(xi1,'r-', 'LineWidth',2);hold on;
% plot(No_u,R_am,'b-', 'LineWidth',2);hold on;
% plot(No_u,R_dc,'r-', 'LineWidth',2);hold on;

%% text of figure
% plot(iter,history.dual,'b-', 'LineWidth',2);hold on;
% 
% xlabel('Number of PG Iterations');
% ylabel('Objective Function value of problem (26)');
% ylabel({'$f(\theta_m)$ in D2'},'Interpreter','latex');
% legend('M = 2','M = 5','M = 10');
% grid on;


