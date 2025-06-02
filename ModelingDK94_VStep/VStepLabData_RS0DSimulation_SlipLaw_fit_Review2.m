clear all
clc

% simulate RS response for a 0d model

%% set up rate and state parameters

Simulation_TimeStep = 5000;

%% material setup
ct = 3000;
cl = 5600;
rho = 2800;
mu = rho*ct*ct;
lambda = rho*cl*cl-2*mu;
nu = (cl^2-2*ct^2)/(2*(cl^2-ct^2));
G =mu;
Cs = ct;

savedata_flag = 1;


%% Friction parameter, all in SI, my test
% Dc = 1e-6; %m
% G = mu; % Pa

% a = 0.01;
% b = 0.015;
% mu0 = 0.6;
% V0 = 1e-6;
% Theta0 = Dc/V0;

Dc = 6.8e-6; %m
a = 0.0546;
b = 0.0558;
mu0 = 0.765;
V0 = 1e-6;
Theta0 = Dc/V0;
sig = 2.5e6;%Pa
chi = 502e6;%Pa, indentation strength, Acrylic Plastic from D&K1994, pg291
%chi = 402e6;
Ar_Factor = 1;

%c factor in Nagata Law
%c = 0.7;
%c = 0;
c_test = 0;
c = c_test;

% loading parameter
%Kshear = -0.5*(sig/Dc); 
%Kshear = -(0.5)*(sig/Dc); 
Kshear = -0.1*1e12;
%Kshear = -0.28*1e12;
Knormal = 0;
%v_load = 1e-6;

% tau_rate = -Kshear*Vloading;
% sigma_rate = 0;

%% Sylvain's parameters

% n = 12;
% m = 3;
% p1 = m*0.98;
% mu0 = 0.765;
% a = mu0/n;
% %a = 0.05;
% b = m/p1*a;
% 
% lambda1=1.2e2;         % reciprocal of characteristic strain
% W=5e-4; 
% 
% Dc = 2*W/(lambda1*p1); %m
% 
% V0 = 1e-6;
% Theta0 = Dc/V0;
% sig = 2.5e6;%Pa
% chi = 530e6;%Pa, indentation strength, Acrylic Plastic from D&K1994, pg291
% Ar_Factor = 1;
% 
% % c factor in Nagata Law
% c = 0.7;
% %c = 0;
% 
% 
% % loading parameter
% %Kshear = -0.5*(sig/Dc); 
% Kshear = -1e14; 
% Knormal = 0;


%% Initial condition


sigma_init = sig;
V_init = V0; % set an very small number to make the initial estimate of time step possible
D_init = 0;
Theta_init = Theta0;

tau_init = 1.0*mu0*sigma_init*(V_init/V0).^(a/mu0)*(Theta_init/Theta0).^(b/mu0);
%% Threhold

TimeStep_Threhold = 0.01;
% Estimated Timestep =  TimeStep_Threhold*Dc/Vmax 
%just as an estimation, would not affect the accuracy which is controlled by ode45;

%% Initialize array

% results are at the end of the time step
%dt_everyStep = [];
t_end_results = [];

tau_results = [];
sigma_results = [];
V_results = [];
D_results = [];
Theta_results = [];
t_results = [];
LoadPointD_results = [];

% assign initial condition
tau_ext_start = tau_init;
sigma_start = sigma_init;
v_start = V_init;
slip_start = D_init;
theta_start = Theta_init;
LoadPointD_start = 0;

t_now = 0;

t_start = t_now;

%% max_dt
max_dt = 100*Dc/v_start;%s

%% setup different stage
% Stg1_Dend = 1e-4 + 2e-6;
% Stg2_Dend = 2e-4;
% Stg3_Dend = 3e-4 + 1e-6;
% Stg4_Dend = 4e-4;

Stg1_Dend = 1e-4;
Stg2_Dend = 2e-4;
Stg3_Dend = 3e-4;
Stg4_Dend = 4e-4;



%for ind_step = 1:Simulation_TimeStep  
step_count = 0;
while(slip_start<=Stg4_Dend)
    %ind_step
    step_count = step_count + 1;
    ind_step = step_count;
    
    % control loading speed
    if (slip_start >= 0 && slip_start <= Stg1_Dend) || (slip_start > Stg2_Dend && slip_start <= Stg3_Dend)
        v_load = 1e-6;
    elseif (slip_start > Stg1_Dend && slip_start <= Stg2_Dend) || (slip_start > Stg3_Dend && slip_start <= Stg4_Dend)
        v_load = 0.1e-6;
    else
        error('something wrong')
    end
    
    dt_now = TimeStep_Threhold*Dc/v_start;
    dt_now = min([dt_now max_dt]);
    

%     [tau_ext_end,sigma_end,v_end,slip_end,theta_end,MinDt] = FindVEnd_ode45_Mult_MinDt_1D(tau_ext_start,sigma_start,...
%             v_start,slip_start,theta_start,...
%             Kshear,Knormal,...
%             dt_now,Dc,...
%             G,Cs,a,b,...
%             mu0,V0,Theta0,...
%             v_load);

      [tau_ext_withinStep,sigma_withinStep,v_withinStep,slip_withinStep,theta_withinStep,t_withinStep]...
            = FindVEnd_ode45_Mult_MinDt_1D_OutPutInBetween_SlipNagata(tau_ext_start,sigma_start,...
            v_start,slip_start,theta_start,...
            Kshear,Knormal,...
            dt_now,Dc,...
            G,Cs,a,b,...
            mu0,V0,Theta0,...
            v_load,c);
        
     tau_results = cat(1,tau_results,tau_ext_withinStep);
     sigma_results = cat(1,sigma_results,sigma_withinStep);
     V_results = cat(1,V_results,v_withinStep);
     D_results = cat(1,D_results,slip_withinStep);
     Theta_results = cat(1,Theta_results,theta_withinStep);
     t_results = cat(1,t_results,t_withinStep+t_start);
     LoadPointD_results = cat(1,LoadPointD_results,LoadPointD_start +...
         t_withinStep.*v_load);
     
     slip_lastEnd = slip_start;
     
     tau_ext_start = tau_ext_withinStep(end);
     sigma_start = sigma_withinStep(end);
     v_start = v_withinStep(end);
     slip_start = slip_withinStep(end);
     theta_start = theta_withinStep(end); 
     t_start = t_start + t_withinStep(end);
     LoadPointD_start = LoadPointD_results(end);
    
end    

%% plot results just friction
colors = load('ColorScheme_MatlabOrder.mat');
mu_results = tau_results./sigma_results;

figure(1)

plot(D_results*1e6,mu_results,...
    'linewidth',3,'Linestyle','-','Color',colors.c(1,:))

xlabel('Displacement, \mu m')
ylabel('COF')
ylim([0.5 0.9])

set(gca,'Fontsize',20,'Fontweight','bold')
set(gca, 'FontName', 'Helvetica')
set(gcf, 'Renderer', 'Painters');% make eps clear


%% plot model against data


D_series = D_results;
Theta_series = Theta_results;
 

%% set the applied V condition in the 4 stages, m

V_series = nan(size(D_series));
V_series(D_series<=Stg1_Dend | (D_series>Stg2_Dend & D_series<=Stg3_Dend)) = 1e-6;
V_series((D_series>Stg1_Dend & D_series<=Stg2_Dend) | (D_series>Stg3_Dend & D_series<=Stg4_Dend)) = 0.1e-6; 


tau_series = mu0*sig.*(V_results./V0).^(a/mu0).*(Theta_series./Theta0).^(b/mu0);
f_series = tau_series./sig;
A_series = Ar_Factor*(mu0*sig/chi)*(Theta_series./Theta0).^(b/mu0);


%% plot model against DK94 fig8 data
colors = load('ColorScheme_MatlabOrder.mat');
ImportDK94_friction;
ImportDK94_PA

D_DK94Friction = DK94Fig8DataDigitizedFriction.D_DK94Friction;
f_DK94Friction = DK94Fig8DataDigitizedFriction.f_DK94Friction;
D_DK94PA = DK94Fig8DataDigitizedPercentArea.D_DK94PA;
PA_DK94PA = DK94Fig8DataDigitizedPercentArea.PA_DK94PA;

figure(3)

yyaxis left
plot(D_DK94Friction,f_DK94Friction,...
    'linewidth',3,'Linestyle','-','Color','k')
hold on
plot(D_series.*1e6 + 400,f_series,'linewidth',3,'Linestyle','-','Color',colors.c(1,:))
hold off

xlabel('Block Displacement, \mu m')
ylabel('COF')
ylim([0.5 0.9])

yyaxis right
plot(D_DK94PA,PA_DK94PA,...
    'linewidth',3,'Linestyle','-','Color','k')
hold on
plot(D_series.*1e6 + 400,A_series*100,'linewidth',3,'Linestyle','-','Color',colors.c(2,:))
hold off

%xlabel('Displacement, \mu m')
ylabel('area percentage')
ylim([0.3 0.8])

xlim([400 800])
set(gca,'Fontsize',20,'Fontweight','bold')
set(gca, 'FontName', 'Helvetica')
set(gcf, 'Renderer', 'Painters');% make eps clear


%% plot model against DK94 fig8 data, using loadpoint velocity

figure(4)

yyaxis left
plot(D_DK94Friction,f_DK94Friction,...
    'linewidth',3,'Linestyle','-','Color','k')
hold on
plot(LoadPointD_results.*1e6 + 400,f_series,'linewidth',3,'Linestyle','-','Color',colors.c(1,:))
hold off

xlabel('Loadpoint Displacement, \mu m')
ylabel('COF')
ylim([0.5 0.9])

yyaxis right
plot(D_DK94PA,PA_DK94PA,...
    'linewidth',3,'Linestyle','-','Color','k')
hold on
plot(LoadPointD_results.*1e6 + 400,A_series*100,'linewidth',3,'Linestyle','-','Color',colors.c(2,:))
hold off

%xlabel('Displacement, \mu m')
ylabel('area percentage')
ylim([0.3 0.8])

xlim([400 800])
set(gca,'Fontsize',20,'Fontweight','bold')
set(gca, 'FontName', 'Helvetica')
set(gcf, 'Renderer', 'Painters');% make eps clear


%% save the whole workspace

if savedata_flag == 1
    save(['vstep_workspace_c',num2str(c,'%.2f'),'_SlipLaw_Review2.mat'])
end