function [T, U, V] = ODE45(f, g, a, b, n, u0, v0)
   
   alpha = [0, 1/4, 3/8, 12/13, 1, 1/2];
    beta = [ ...
        0           0           0           0         0       0;
        1/4         0           0           0         0       0;
        3/32        9/32        0           0         0       0;
        1932/2197  -7200/2197   7296/2197   0         0       0;
        439/216    -8           3680/513   -845/4104  0       0;
        -8/27       2          -3544/2565   1859/4104 -11/40  0 ];

    c = [16/135  0  6656/12825  28561/56430  -9/50  2/55]; 


    h = (b - a) / n;
    T = linspace(a, b, n+1);
    U = zeros(1, n+1);  
    V = zeros(1, n+1);  
    U(1) = u0;
    V(1) = v0;

    for i = 1:n
        t = T(i);
        u = U(i);
        v = V(i);

        ku1 = f(t, u, v);
        ku2 = f(t + alpha(2)*h, u + h*beta(2,1)*ku1, v + h*beta(2,1)*g(t, u, v));
        ku3 = f(t + alpha(3)*h, u + h*(beta(3,1)*ku1 + beta(3,2)*ku2), ...
                                 v + h*(beta(3,1)*g(t,u,v) + beta(3,2)*g(t,u,v)));
        ku4 = f(t + alpha(4)*h, u + h*(beta(4,1)*ku1 + beta(4,2)*ku2 + beta(4,3)*ku3), ...
                                 v + h*(beta(4,1)*g(t,u,v) + beta(4,2)*g(t,u,v) + beta(4,3)*g(t,u,v)));
        ku5 = f(t + alpha(5)*h, u + h*(beta(5,1)*ku1 + beta(5,2)*ku2 + beta(5,3)*ku3 + beta(5,4)*ku4), ...
                                 v + h*(beta(5,1)*g(t,u,v) + beta(5,2)*g(t,u,v) + beta(5,3)*g(t,u,v) + beta(5,4)*g(t,u,v)));
        ku6 = f(t + alpha(6)*h, u + h*(beta(6,1)*ku1 + beta(6,2)*ku2 + beta(6,3)*ku3 + beta(6,4)*ku4 + beta(6,5)*ku5), ...
                                 v + h*(beta(6,1)*g(t,u,v) + beta(6,2)*g(t,u,v) + beta(6,3)*g(t,u,v) + beta(6,4)*g(t,u,v) + beta(6,5)*g(t,u,v)));

        kv1 = g(t, u, v);
        kv2 = g(t + alpha(2)*h, u + h*beta(2,1)*v, v + h*beta(2,1)*kv1);
        kv3 = g(t + alpha(3)*h, u + h*(beta(3,1)*v + beta(3,2)*kv2), v + h*(beta(3,1)*kv1 + beta(3,2)*kv2));
        kv4 = g(t + alpha(4)*h, u + h*(beta(4,1)*v + beta(4,2)*kv2 + beta(4,3)*kv3), v + h*(beta(4,1)*kv1 + beta(4,2)*kv2 + beta(4,3)*kv3));
        kv5 = g(t + alpha(5)*h, u + h*(beta(5,1)*v + beta(5,2)*kv2 + beta(5,3)*kv3 + beta(5,4)*kv4), v + h*(beta(5,1)*kv1 + beta(5,2)*kv2 + beta(5,3)*kv3 + beta(5,4)*kv4));
        kv6 = g(t + alpha(6)*h, u + h*(beta(6,1)*v + beta(6,2)*kv2 + beta(6,3)*kv3 + beta(6,4)*kv4 + beta(6,5)*kv5), v + h*(beta(6,1)*kv1 + beta(6,2)*kv2 + beta(6,3)*kv3 + beta(6,4)*kv4 + beta(6,5)*kv5));

        U(i+1) = u + h*(c(1)*ku1 + c(3)*ku3 + c(4)*ku4 + c(5)*ku5 + c(6)*ku6);
        V(i+1) = v + h*(c(1)*kv1 + c(3)*kv3 + c(4)*kv4 + c(5)*kv5 + c(6)*kv6);
    end
end