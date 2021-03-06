%*************************************************************************%
%                                                                         %
%  This function computes the local and nonlocal states given the values   %
%  the volume constraint and the boundary condition                        %
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
function value = evaluate_nonlocal_solution(k,v,x,uN)
%%% nonlocal error
    value = ((uN(2*k)-uN(2*k-1))./(x(k+1)-x(k)).*(v'-x(k)) + (uN(2*k-1)))';
