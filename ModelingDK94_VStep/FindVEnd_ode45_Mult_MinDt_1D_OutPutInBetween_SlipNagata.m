% find V for the end of the time step

%%%%%%%%%%%%% Inputs %%%%%%%%%%%%%%%%%%

% All input vector must be 1*N

% tau_ext_start, external shear stress at each element at
% the begining of this timestep, vectors

% Theta_start, state variable at each element at each element at
% the begining of this timestep, vectors

% sigma_start, normal stress at each element at each element at
% the begining of this timestep, vectors

% v_start, slip rate at each element at each element at
% the begining of this timestep, vectors

% slip_start, slip at each element at each element at
% the begining of this timestep, vectors

% Kshear_matrix, Stress transfer matrix for shear stress
% Knormal_matrix, Stress transfer matrix for normal stress

% step_control_coef, the dt = step_control_coef*Dc/[max slip rate], scalar
% Dc, rate-state Dc, vector

% a, rate-state a, vector
% b, rate-state b, vector
% mu0, rate-state mu0, vector
% tau_rate, shear stressing rate at each element, vector
% sigma_rate, normal stressing rate at each element, vector

% c, nagata stress weakening parameter Bhattacharya et al. 2017 JGR


%%%%%%%%%%%%%%%% Outputs %%%%%%%%%%%%%%%%%%%

%tau_ext, external shear stress of all calculated steps

%sigma, normal stress of all calculated steps

%v, slip rate of all calculated steps

%theta, state variable of all calculated steps

%dt_, pre determined timestep length of this time step

function [tau_ext,sigma,v,slip,theta,t_out] = FindVEnd_ode45_Mult_MinDt_1D_OutPutInBetween_SlipNagata(tau_ext_start,sigma_start,...
        v_start,slip_start,theta_start,...
        Kshear_matrix,Knormal_matrix,...
        dt_now,Dc,...
        G,Cs,a,b,...
        mu0,V0,Theta0,...
        v_load,...
        c) 
%%%% Preprocessing
    size_tau_start = size(tau_ext_start);
    size_Kshear = size(Kshear_matrix);
    size_Knormal = size(Knormal_matrix);
    
    if ~isequaln(size(tau_ext_start),size(v_start),size(theta_start),size(sigma_start),size(slip_start)) ...
            || size_tau_start(1) ~= 1
        error(' one or more of the start array have dimensions that are not correct')
    end
    if ~isequaln(size_Kshear,size_Knormal) || length(size(Kshear_matrix)) ~= 2 ...
        || size_Kshear(1) ~= size_Kshear(2)
        error('Green function dimension not correct.')
    end
    if size_Kshear(2) ~= length(v_start)
        error('Green function size does not match the slip rate profile.')
    end

%%%%%%% finding the solution    
    % setup intial v_trial
    
    % building the y0 vection for ode45
    % if there is N element, the total number of row(1 colomn) of y is 5N
    % y(1) - y(N): external shear stress 
    % y(N+1) - y(2N): normal stress
    % y(2N+1) - y(3N): slip rate
    % y(3N+1) - y(4N): slip
    % y(4N+1) - y(5N): theta
    
    N = size_tau_start(2);
    y0 = nan(5*N,1);
    
    y0(1:N,1) = tau_ext_start';
    y0(N+1:2*N,1) = sigma_start';
    y0(2*N+1:3*N,1) = log(v_start'./V0');
    y0(3*N+1:4*N,1) = slip_start';
    y0(4*N+1:5*N,1) = log(theta_start'./Theta0');
    
%     v_trial = v_start;
%     dt_now = min(abs(step_control_coef*Dc./v_trial));%giving an approximate time step
    
    tspan = [0 dt_now];
    
    %options = odeset('RelTol',1e-5,'Stats','on','OutputFcn',@odeplot);
    %options = odeset('RelTol',1e-5);
    options=odeset('Refine',1,'RelTol',3e-7,'InitialStep',1e-3,'MaxStep',3e6); 
    
    % Quasi-dynamic
%     [t_out,y_out] = ode45(@(t,y) RSmarching_RD_mult_1D(t,y,Kshear_matrix,Knormal_matrix,Dc,G,Cs,a,b,mu0,V0,Theta0,v_load),...
%         tspan, y0,options);% add option, absolute tol, relative tol, max time step
    [t_out,y_out] = ode45(@(t,y) RSmarching_RD_mult_1D_SlipNagata(t,y,Kshear_matrix,Knormal_matrix,Dc,G,Cs,a,b,mu0,V0,Theta0,v_load,c),...
        tspan, y0,options);% add option, absolute tol, relative tol, max time step

    
    % Quasi-static
%     [t_out,y_out] = ode45(@(t,y) RSmarching_QS(t,y,Kshear_matrix,Knormal_matrix,Dc,G,Cs,a,b,mu0,V0,Theta0,tau_rate,sigma_rate),...
%         tspan, y0,options);% add option, absolute tol, relative tol, max time step    
    
    
%     tau_ext_end = y_out(end,1:N);
%     sigma_end = y_out(end,N+1:2*N);
%     v_end = V0.*exp(y_out(end,2*N+1:3*N));
%     slip_end = y_out(end,3*N+1:4*N);
%     theta_end = Theta0.*exp(y_out(end,4*N+1:5*N));  
    
    tau_ext = y_out(:,1:N);
    sigma = y_out(:,N+1:2*N);
    v = V0.*exp(y_out(:,2*N+1:3*N));
    slip = y_out(:,3*N+1:4*N);
    theta = Theta0.*exp(y_out(:,4*N+1:5*N));      
    
    %MinDt = min(diff(t_out));
    MinDt = t_out(end) - t_out(end-1);
    
    %dimension is specified with 1:N, so no need for transpose
%     tau_ext_end = tau_ext_end';
%     sigma_end = sigma_end';
%     v_end = v_end';
%     slip_end = slip_end';
%     theta_end = theta_end'; 
    
end