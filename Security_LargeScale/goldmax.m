function [x,fval] = goldmax(f,z1,z2,eps)
%%  黄金分割法求解函数的最值
    % 输入
    % f 待优化函数
    % eps 精度
    %  输出
    % x 最优解 beta 
    % fval 最优解对应的最大值 Rs
    k = 0;
    g = (3-sqrt(5))/2;
    Z1 = z1+g*(z2-z1);
    Z2 = z2-g*(z2-z1);
    F1 = f(Z1);
    F2 = f(Z2);

    while abs(z2-z1)>=eps
        if  F1>F2
            z2 = Z2;
            Z2 = Z1;
            F2 = F1;
            Z1 = z1+g*(z2-z1);
            F1 = f(Z1);

        elseif F1<F2
            z1 = Z1;
            Z1 = Z2;
            F1 = F2;
            Z2 = z2-g*(z2-z1);
            F2 = f(Z2);
        else
            z1 = Z1;
            z2 = Z2;
            Z1 = z1+g*(z2-z1);
            Z2 = z2-g*(z2-z1);
            F1 = f(Z1);
            F2 = f(Z2);
        end
    end
    x = (z1+z2)/2;
    fval = f(x);


end



