function [t, u, v] = NRK4SED(f, g, a, b, n, u0, v0)
%NRK4SED Método de Runge-Kutta de 4ª ordem para um Sistema de EDOs/PVI
%   Resolve um sistema de duas equações diferenciais de 1ª ordem do tipo:
%       u' = f(t, u, v)
%       v' = g(t, u, v)
%   com condições iniciais: u(a) = u0, v(a) = v0
%
%   Fórmulas de iteração:
%     k1u = h * f(t(i), u(i), v(i))
%     k1v = h * g(t(i), u(i), v(i))
%
%     k2u = h * f(t(i) + h, u(i) + k1u, v(i) + k1v)
%     k2v = h * g(t(i) + h, u(i) + k1u, v(i) + k1v)
%
%     k3u = h * f(t(i) + h/2, u(i) + k2u/2, v(i) + k2v/2)
%     k3v = h * g(t(i) + h/2, u(i) + k2u/2, v(i) + k2v/2)
%
%     k4u = h * f(t(i) + h, u(i) + k3u, v(i) + k3v)
%     k4v = h * g(t(i) + h, u(i) + k3u, v(i) + k3v)
%
%     u(i+1) = u(i) + (k1u + 2*k2u + 2*k3u + k4u) / 6
%     v(i+1) = v(i) + (k1v + 2*k2v + 2*k3v + k4v) / 6
%
% INPUT:
%   f  - função do lado direito da equação u' = f(t, u, v)
%   g  - função do lado direito da equação v' = g(t, u, v)
%   a  - limite inferior do intervalo
%   b  - limite superior do intervalo
%   n  - número de subintervalos
%   u0 - condição inicial para u
%   v0 - condição inicial para v
%
% OUTPUT:
%   t  - vetor dos tempos (1 x n+1)
%   u  - vetor das soluções aproximadas de u(t)
%   v  - vetor das soluções aproximadas de v(t)
%
%   28/03/2025  Afonso Mariz Luís         a2022127026@isec.pt
%   28/03/2025  Diogo Bento Santos        a2022108969@isec.pt
%   28/03/2025  João Manuel Almeida Nunes a2022122159@isec.pt

    h = (b - a) / n;
    t = a:h:b;
    u = zeros(1, n+1);
    v = zeros(1, n+1);
    u(1) = u0;
    v(1) = v0;

    for i = 1:n
        % k1
        k1u = h * f(t(i), u(i), v(i));
        k1v = h * g(t(i), u(i), v(i));

        % k2
        k2u = h * f(t(i) + h, u(i) + k1u, v(i) + k1v);
        k2v = h * g(t(i) + h, u(i) + k1u, v(i) + k1v);

        % k3
        k3u = h * f(t(i) + h/2, u(i) + k2u/2, v(i) + k2v/2);
        k3v = h * g(t(i) + h/2, u(i) + k2u/2, v(i) + k2v/2);

        % k4
        k4u = h * f(t(i) + h, u(i) + k3u, v(i) + k3v);
        k4v = h * g(t(i) + h, u(i) + k3u, v(i) + k3v);

        % Atualização das soluções
        u(i+1) = u(i) + (k1u + 2*k2u + 2*k3u + k4u)/6;
        v(i+1) = v(i) + (k1v + 2*k2v + 2*k3v + k4v)/6;
    end
end
