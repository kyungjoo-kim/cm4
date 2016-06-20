function [uN,uL,xN_plot,errN,errL] = state(xL,hL,theta_local_left, ...
                                           xN,hN,epsilon,n,theta_nonlocal_right,...
                                           reuse,print,accuracy,test)
%*************************************************************************%
%                                                                         %
%  This function computes the local and nonlocal states given the values  %
%  the volume constraint and the boundary condition                       %
%                                                                         %
%  Author: Marta D'Elia                                                   %
%                                                                         %
%  Modified: 06-01-2016                                                   %
%                                                                         %
%*************************************************************************%

    theta_local_right = exact_solution(xL(end),test);
    [AL,bL,uL,errL] = local_problem(xL,hL,epsilon,...
                                    theta_local_left,theta_local_right,...
                                    reuse,print,accuracy,test); 

    theta_nonlocal_left = zeros(2*n,1);
    theta_nonlocal_left_tmp = exact_solution(xN(1:n+1),test); 
    theta_nonlocal_left(1:2:end) = theta_nonlocal_left_tmp(1:end-1);
    theta_nonlocal_left(2:2:end) = theta_nonlocal_left_tmp(2:end); 

    [AN,bN,uN,xN_plot,errN] = nonlocal_problem(xN,hN,epsilon,n,...
                                               theta_nonlocal_left,theta_nonlocal_right,...
                                               reuse,print,accuracy,test);
    
