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
function [errN,errL] = run_coupling_optimize(NN,NL,epsilon,test)
    close all;
    more off;
    
    %%% domains ---------------------------------------------------------------
    problem_domain  = [ 0    1.75];
    nonlocal_domain = [ 0    1   ];
    local_domain    = [ 0.75 1.75];
    overlap_domain  = [ local_domain(1) nonlocal_domain(2)+epsilon];
    [xN,xL,n,hN,hL] = domain2(NN,NL,epsilon,nonlocal_domain,local_domain); 

    fprintf('epsilon = %f, hN = %f, hL = %f\n', epsilon, hN, hL);
    if (epsilon < hN)
        warning('nonlocal neighborhood (epsilon = %f) is smaller than  mesh size (h = %f)\n', ...
                epsilon, hN);
    end

    %% overlap region in nonlocal model
    xNr = xN(end-n:end);

    %%% initializing the control ----------------------------------------------
    accuracy = true;
    print    = false;
    reuse    = false;
    use_coarse_solution = true;
    
    theta_local_left = 1.0; %% arbitrary value
    theta_nonlocal_right = zeros(2*n,1) + 1.5; %% arbitrary value  

    theta  = [theta_local_left;theta_nonlocal_right]; % control vector
    theta_initial = theta;

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
    
    %%% solving --------------------------------------------------------------- 
    % computing the local and nonlocal matrix once an for all
    [uN,uL,xN_plot]= state(xL,hL,theta_local_left,...
                           xN,hN,epsilon,n,theta_nonlocal_right,...
                           reuse,print,accuracy,test);

    figure;
    plot(xL,uL,'k-',xN_plot,uN,'r-','Linewidth',4);
    vline=get(gca,'ylim');
    hold on;
    lpos(1:2) = nonlocal_domain(2);        plot(lpos,vline,'r');
    lpos(1:2) = nonlocal_domain(2)+epsilon;plot(lpos,vline,'m');
    lpos(1:2) = local_domain(1);           plot(lpos,vline,'k:');
    hold off;

    reuse = true;
    tic;
    options = optimset('Display','iter', ...
                       'TolFun',1.e-10, ... 
                       'TolX',1.e-10, ...
                       'MaxIter',1000, ...
                       'MaxFunEvals',1000, ...
                       'GradObj','off');
    [theta,fval,info,output] = ...
        fminunc(@(t)functional(t,...
                               xL,hL,...
                               xN,hN,epsilon,n,...
                               overlap_domain,...
                               reuse,test),...
                theta_initial,options);

    info
    output
    
    theta_local_left = theta(1);
    theta_nonlocal_right = theta(2:end);

    [uN,uL,xN_plot,errN,errL]= state(xL,hL,theta_local_left,...
                                     xN,hN,epsilon,n,theta_nonlocal_right,...
                                     reuse,print,accuracy,test);
    toc
    
    figure;
    plot(xL,uL,'k-',xN_plot,uN,'r-','Linewidth',4);
    vline=get(gca,'ylim');
    hold on;
    lpos(1:2) = nonlocal_domain(2);        plot(lpos,vline,'r');
    lpos(1:2) = nonlocal_domain(2)+epsilon;plot(lpos,vline,'m');
    lpos(1:2) = local_domain(1);           plot(lpos,vline,'k:');
    hold off;

%----- FUNCTIONS -----%

function [Jvalue] = functional(t,...
                               xL,hL,...
                               xN,hN,epsilon,n,...
                               overlap_domain,...
                               reuse,test)
    [uN,uL,xN_plot] = state(xL,hL,t(1),...
                            xN,hN,epsilon,n,t(2:end),...
                            reuse,false,false,test);

    idx = find(xN_plot >= overlap_domain(1) & ...
               xN_plot <= overlap_domain(2));

    xN_tmp = xN_plot(idx);
    xN_overlap = xN_tmp(2:end);

    uN_tmp = uN(idx);
    uN_overlap = uN_tmp(2:end);

    % evaluate uL on xN_overlap    
    n_overlap = length(xN_overlap);
    uL_overlap = zeros(n_overlap,1);
    for i=1:n_overlap
        xN_overlap(i);
        lower_bound = find(xL <= xN_overlap(i));
        k = lower_bound(end);
        v = xN_overlap(i);
        uL_overlap(i) = evaluate_local_solution(k,v,xL,uL);
    end
    
    % computation of J
    Jvalue = 0;
    for i=1:n_overlap/2
        Jvalue = Jvalue + quadgk(@(v)gaussJ(i,v,xN_overlap,uL_overlap,uN_overlap),...
                                 xN_overlap(2*k-1),xN_overlap(2*k),'RelTol',1e-10);
    end
    Jvalue = 0.5*Jvalue;
    
function value = gaussJ(k,v,x,uL,uN)  
%%% evaluation of the L2 norm |un-ul|^2
    k1 = 2*k-1; k2 = 2*k;

    ul = ((uL(k2)-uL(k1))./(x(k2)-x(k1)).*(v'-x(k1)) + (uL(k1)))';
    un = ((uN(k2)-uN(k1))./(x(k2)-x(k1)).*(v'-x(k1)) + (uN(k1)))';
    value = (un-ul).^2;