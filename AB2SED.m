function [t, u, v] = AB2SED(f, g, a, b, n, u0, v0)
% AB2SED Método de Adams-Bashforth de 2ª ordem para um Sistema de EDOs/PVI
%   Resolve um sistema de duas equações diferenciais de 1ª ordem do tipo:
%       u' = f(t, u, v)
%       v' = g(t, u, v)
%   com condições iniciais: u(a) = u0, v(a) = v0
%
%   Fórmulas de iteração:
%     u(i+1) = u(i) + h/2 * (3*f(t(i), u(i), v(i)) - f(t(i-1), u(i-1), v(i-1)))
%     v(i+1) = v(i) + h/2 * (3*g(t(i), u(i), v(i)) - g(t(i-1), u(i-1), v(i-1)))
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

    % Primeira iteração usando o método de Euler (ou outro método inicial)
    u(2) = u(1) + h * f(t(1), u(1), v(1));
    v(2) = v(1) + h * g(t(1), u(1), v(1));

    for i = 2:n
        % Aplicando a fórmula de Adams-Bashforth de 2ª ordem
        u(i+1) = u(i) + h/2 * (3*f(t(i), u(i), v(i)) - f(t(i-1), u(i-1), v(i-1)));
        v(i+1) = v(i) + h/2 * (3*g(t(i), u(i), v(i)) - g(t(i-1), u(i-1), v(i-1)));
    end
end
