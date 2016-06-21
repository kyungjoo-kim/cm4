% The code solves a local 1D poisson problem. 
%
% [errL, AL, bL, uL] = run_local_problem(N, epsilon, test)
% * Input
%   - N: # of elements
%   - epsilon: interaction radius, dummy for test 1 - 4
%   - test: problem id
%     0 -> a source function with discontinuity
%     1 -> problem with a manufactured solution u = x
%     2 -> problem with a manufactured solution u = x^2
%     3 -> problem with a manufactured solution u = x^3
%     4 -> problem with a manufactured solution u = x^2 - x^4
% * Output 
%   - errL: local error with respect to manufactured solutions
%   - AL: local stiffness matrix
%   - bL: load vector
%   - uL: solution vector
%
%  Author: Marta D'Elia
function [errL,AL,bL,uL] = run_local_problem(N,epsilon,test)
    close all;
    more off;

    %%% domains ---------------------------------------------------------------

    % non local dummy
    nonlocal_domain = 0;

    if (test == 0)
        %% for integration of source term it needs nonzero epsilon
    else
        %% otherwise, epsilon is dummy and should be zero
        epsilon = 0;
    end
    
    % local domain
    local_domain = [ 0 1.75 ];
    [xN,xL,n,h] = domain(N,epsilon,nonlocal_domain,local_domain);

    %%% initializing the control ----------------------------------------------
    accuracy = true;
    print    = false;
    reuse    = false;

    theta_local_left  = exact_solution(xL(1),  test);
    theta_local_right = exact_solution(xL(end),test);

    tic;
    [AL,bL,uL,errL] = local_problem(xL,h,epsilon,...
                                    theta_local_left,theta_local_right,...
                                    reuse,print,accuracy,test);
    toc;
    
    figure;
    plot(xL,uL,'k-+');

    SL = sparse(AL);

    figure;
    spy(SL);

    fprintf('nnz(A) = %d\n', nnz(SL));    