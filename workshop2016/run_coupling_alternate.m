%*************************************************************************%
%                                                                         %
%  This function solves the optimization-based local-nonlocal coupling    %
%  problem using the built-in Matlab function fminunc.m                   %
%                                                                         %
%  Author: Marta D'Elia                                                   %
%                                                                         %
%  Modified: 01-06-2016                                                   %
%                                                                         %
%  NOTE 1: h is THE SAME for nonlocal and local problems                  %
%                                                                         %
%  NOTE 2: nonlocal discretization - piecewise linear Disc Galerkin       %
%          local discretization - piecewise linear Cont Galerkin          %
%                                                                         %
%*************************************************************************%
function [errN,errL] = run_coupling_alternate(NN,NL,epsilon,test)
    close all;
    more off;
    
    %%% domains ---------------------------------------------------------------
    problem_domain  = [ 0    1.75];
    nonlocal_domain = [ 0    1   ];
    local_domain    = [ 0.75 1.75];
    [xN,xL,n,hN,hL] = domain2(NN,NL,epsilon,nonlocal_domain,local_domain);

    printf('epsilon = %f, hN = %f, hL = %f\n', epsilon, hN, hL);
    if (epsilon < hN)
        warning('nonlocal neighborhood (epsilon = %f) is smaller than  mesh size (h = %f)\n', ...
                epsilon, hN);
    end 
    
    %% overlap region in nonlocal model
    xNr = xN(end-n:end);

    %%% initializing the control ----------------------------------------------
    accuracy   = true;
    print      = true;
    reuse      = false;
    use_coarse_solution = false;
    t          = 1.0; %% relaxation parameter (0,1]
    
    % local model boundary
    theta_local_left = 1.0; %% arbitrary value
    theta_local_right = exact_solution(xL(end),test); 

    % nonlocal model boundary
    theta_nonlocal_left = zeros(2*n,1);
    theta_nonlocal_left_tmp = exact_solution(xN(1:n+1),test);
    theta_nonlocal_left(1:2:end) = theta_nonlocal_left_tmp(1:end-1);
    theta_nonlocal_left(2:2:end) = theta_nonlocal_left_tmp(2:end);
    
    theta_nonlocal_right = zeros(2*n,1) + 1.5; %% arbitrary value
    theta_nonlocal_right_tmp = zeros(n+1,1);

    % use coarse solution for initial guess
    if (use_coarse_solution)
        [dummy,xA,dummy,hA] = domain(NL,epsilon,0,problem_domain);

        theta_all_left  = exact_solution(xA(1),  test);
        theta_all_right = exact_solution(xA(end),test);   

        [AA,bA,uA,errA] = local_problem(xA,hA,epsilon, ...
                                        theta_all_left,theta_all_right, ...
                                        reuse,print,accuracy,test);

        figure;
        plot(xA,uA,'k-');
        vline=get(gca,'ylim');
        hold on;
        lpos(1:2) = nonlocal_domain(2);        plot(lpos,vline,'r');
        lpos(1:2) = nonlocal_domain(2)+epsilon;plot(lpos,vline,'m');
        lpos(1:2) = local_domain(1);           plot(lpos,vline,'k:');
        hold off;

        lower_bound = find(xA <= xL(1));
        k = lower_bound(end);
        v = xL(1);
        theta_local_left = evaluate_local_solution(k,v,xA,uA);
        
        for i=1:n+1
            lower_bound = find(xA <= xNr(i));
            k = lower_bound(end);
            v = xNr(i); 

            theta_nonlocal_right_tmp(i) = evaluate_local_solution(k,v,xA,uA);
        end
        theta_nonlocal_right(1:2:end) = theta_nonlocal_right_tmp(1:end-1);
        theta_nonlocal_right(2:2:end) = theta_nonlocal_right_tmp(2:end);
    end
    
    [AL,bL,uL,errL] = local_problem(xL,hL,epsilon,...
                                    theta_local_left,theta_local_right,...
                                    reuse,print,accuracy,test);

    [AN,bN,uN,xN_plot,errN] = nonlocal_problem(xN,hN,epsilon,n,...
                                               theta_nonlocal_left,theta_nonlocal_right,...
                                               reuse,print,accuracy,test);

    % initial condition plot
    figure;
    plot(xL,uL,'k-',xN_plot,uN,'r-','Linewidth',4); 
    vline=get(gca,'ylim');
    hold on;
    lpos(1:2) = nonlocal_domain(2);        plot(lpos,vline,'r');
    lpos(1:2) = nonlocal_domain(2)+epsilon;plot(lpos,vline,'m');
    lpos(1:2) = local_domain(1);           plot(lpos,vline,'k:');
    hold off;
    
    %%% loop in solving
    reuse = true;
    converged = false;
    rel_thres = 1.0e-4;
    iter = 0; max_iter = 100;

    figure(10);

    tic;
    while ~converged

        theta_local_left_prev = theta_local_left;
        
        % adjust boundary condition
        lower_bound = find(xN <= xL(1));
        k = lower_bound(end);
        v = xL(1);
        theta_local_left = ...
            (1-t)*theta_local_left + ...
            t*    evaluate_nonlocal_solution(k,v,xN,uN);
        
        theta_nonlocal_right_tmp(1:end-1) = theta_nonlocal_right(1:2:end);
        theta_nonlocal_right_tmp(2:end)   = theta_nonlocal_right(2:2:end);
        for i=1:n+1
            lower_bound = find(xL <= xNr(i));
            k = lower_bound(end);
            v = xNr(i); 

            theta_nonlocal_right_tmp(i) = ...
                (1-t)*theta_nonlocal_right_tmp(i) + ...
                t    *evaluate_local_solution(k,v,xL,uL);
        end
        theta_nonlocal_right(1:2:end) = theta_nonlocal_right_tmp(1:end-1);
        theta_nonlocal_right(2:2:end) = theta_nonlocal_right_tmp(2:end);

        % run the problem with updated boundary conditions
        [AL,bL,uL,errL] = local_problem(xL,hL,epsilon,...
                                        theta_local_left,theta_local_right,...
                                        reuse,print,accuracy,test);

        [AN,bN,uN,xN_plot,errN] = nonlocal_problem(xN,hN,epsilon,n,...
                                                   theta_nonlocal_left,theta_nonlocal_right,...
                                                   reuse,print,accuracy,test);

        if (print)
            figure(10);
            plot(xL,uL,'k-',xN_plot,uN,'r-','Linewidth',4);         
            hold on;
            vline=get(gca,'ylim');
            lpos(1:2) = nonlocal_domain(2);        plot(lpos,vline,'r');
            lpos(1:2) = nonlocal_domain(2)+epsilon;plot(lpos,vline,'m');
            lpos(1:2) = local_domain(1);           plot(lpos,vline,'k:');
            hold off;
            printf('[Enter] to continue ...\n');
            pause;
        end
        
        % simple
        rate = abs(theta_local_left - theta_local_left_prev)/ ...
               abs(theta_local_left_prev);
        if rate < rel_thres
            converged = true;
        end

        % count iteration
        iter = iter + 1;
        if (iter > max_iter)
            break;
        end

        printf('iter = %d, rate = %f, err(local,nonlocal) = %f, %f\n' ...
               , iter, rate, errL, errN);
    end
    toc;
    
    if converged
        printf('converged\n');
    else
        printf('reached max number of iteration %d\n', max_iter);        
    end
    
    figure;
    plot(xL,uL,'k-',xN_plot,uN,'r-','Linewidth',4);         
    vline=get(gca,'ylim');
    hold on;
    lpos(1:2) = nonlocal_domain(2);        plot(lpos,vline,'r');
    lpos(1:2) = nonlocal_domain(2)+epsilon;plot(lpos,vline,'m');
    lpos(1:2) = local_domain(1);           plot(lpos,vline,'k:');
    hold off;
end