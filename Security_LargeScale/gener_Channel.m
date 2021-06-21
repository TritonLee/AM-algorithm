%% channel generation for beamforming design 
function [BF] = gener_Channel(epsilon,P_t,alp_e,sigma_e,N,M)
    ka = zeros(1,M);
    eps = log(epsilon^(-1));

    G = 1/sqrt(2)*(randn(M,N)+j*randn(M,N));  % random generated channel vector after channel estimation 
    %g_e = randn(1,N)+j*randn(1,N);
    W = (G'*G)^(-1)*G';          % N times M
    for m = 1:M   % for each cluster 
        uni_wm = W(:,m)/norm(W(:,m))^2;   % normalized beamforming vector 
        S_m = uni_wm*uni_wm';  
        Jtemp = zeros(N,N);
        for jj = 1:M
            if jj ~= m
                Jtemp = Jtemp + uni_wm*uni_wm';  
            end
        end
        J(:,:,m) = Jtemp;   % the obtained complementary beamformer for cluster-m
        deno = (1+eps+sqrt(2*eps))*trace(S_m);
        nomi = sigma_e/(P_t*alp_e) + (trace(Jtemp) - sqrt(2*eps)*norm(Jtemp,'fro'))/M;
        ka(m) = real(nomi/deno);   
    end
    BF = W;
    
end




