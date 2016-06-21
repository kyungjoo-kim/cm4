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
function value = exact_solution(x,test)
    switch test
      case 0
        %% this is used for creating boundary condition only
        value = x.*0;
      case 1
        value = x.*1;
      case 2
        value = x.^2;
      case 3
        value = x.^3;
      case 4
        value = x.^2 - x.^4;
      case 5
        %% this is used for creating boundary condition only
        value = x.*0;
      case 6
        %% this is used for creating boundary condition only
        value = x.*0;
    end
