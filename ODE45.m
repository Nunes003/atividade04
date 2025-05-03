function [t, u, v] = ODE45(f, g, a, b, n, u0, v0)
%ODE45 Método de resolução numérica de EDOs usando o solver embutido do MATLAB
%   Resolve um sistema de duas equações diferenciais de 1ª ordem do tipo:
%       u' = f(t, u, v)
%       v' = g(t, u, v)
%   com condições iniciais: u(a) = u0, v(a) = v0
%
% INPUT:
%   f  - função do lado direito da equação u' = f(t, u, v)
%   g  - função do lado direito da equação v' = g(t, u, v)
%   a  - limite inferior do intervalo
%   b  - limite superior do intervalo
%   n  - número de subintervalos (define os pontos em tspan)
%   u0 - condição inicial para u
%   v0 - condição inicial para v
%
% OUTPUT:
%   t  - vetor dos tempos (1 x n+1)
%   u  - vetor das soluções aproximadas de u(t)
%   v  - vetor das soluções aproximadas de v(t)
%
% NOTA:
%   Usa o solver adaptativo ode45 do MATLAB com o vetor de estado Y = [u; v].
%   A saída u e v são transpostas para manter a consistência (vetores linha).
%
%   28/03/2025  Afonso Mariz Luís         a2022127026@isec.pt
%   28/03/2025  Diogo Bento Santos        a2022108969@isec.pt
%   28/03/2025  João Manuel Almeida Nunes a2022122159@isec.pt

    % Define o sistema vetorial como uma função anónima
    sistema = @(t, Y) [f(t, Y(1), Y(2)); g(t, Y(1), Y(2))];

    % Cria vetor de tempos com n+1 pontos igualmente espaçados
    tspan = linspace(a, b, n+1);
    Y0 = [u0; v0];  % Condições iniciais

    % Chamada ao ode45
    [t, Y] = ode45(sistema, tspan, Y0);

    % Separação dos resultados
    u = Y(:, 1).';  % transposto para vetor linha
    v = Y(:, 2).';
end
