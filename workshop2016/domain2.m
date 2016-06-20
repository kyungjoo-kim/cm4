%*************************************************************************%
%                                                                         %
%  Problem domain setup for local and nonlocal                            %
%                                                                         %
%  Author: Marta D'Elia                                                   %
%                                                                         %
%  Modified: 06-01-2016                                                   %
%                                                                         %
%  NOTE 1: h is THE SAME for nonlocal and local problems                  %
%                                                                         %
%  NOTE 2: nonlocal discretization - piecewise linear Disc Galerkin       %
%          local discretization - piecewise linear Cont Galerkin          %
%                                                                         %
%*************************************************************************%
% N - # of points in nonlocal domain
% epsilon - nonlocal operator neighborhood size
% nonlocal - nonlocal region
% local - local region
function [xN,xL,n,hN,hL] = domain2(NN,NL,epsilon,nonlocal,local)
    a = 0; b = 0; h = 0;

    h = 0;
    N = NN;
    if nonlocal == 0
        printf('domain:: nonlocal domain is not createdn\n');
        xN = 0; n = 0;
    else
        a = nonlocal(1);
        b = nonlocal(2);
        h = (b-a)/N;
        
        % # of points in epsilon
        if (mod(epsilon,h) == 0)
            n = round(epsilon/h);
        else
            n = ceil(epsilon/h);
        end
        
        xN = [a-epsilon [a-(n-1)*h:h:b+(n-1)*h] b+epsilon];
    end
    hN = h;

    h =0;
    N = NL;
    if local == 0
        printf('domain:: local domain is not createdn\n');
        xL = 0;
    else
        a = local(1);
        b = local(2);
        h = (b-a)/N;

        xL = [a:h:b];
    end
    hL = h;