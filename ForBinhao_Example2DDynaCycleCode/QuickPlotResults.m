%% plot Slip rate evoltion of the whole simulation

clear all
clc

%% load plot colors

% load default colors
DClr = load('ColorScheme_MatlabOrder.mat');

%% load simulation data
% slow simulation
%Sim1 = load('MatchSvet_SaveState_Ru6.76_Rb0.11_sig4.50_b0.0180_a0.0160_L2.6e-07_QDFlag0_SimStep120000.mat');
%Sim1 = load('MatchSvet_SaveState_Ru6.76_Rb0.11_sig4.50_b0.0180_a0.0160_L2.6e-07_QDFlag0_SimStep120000_FSym.mat');
%Sim1 = load('MatchSvet_SaveState_Ru6.74_Rb0.11_sig5.00_b0.0180_a0.0160_L2.9e-07_QDFlag0_SimStep120000_FSym.mat');
%Sim1 = load('MatchSvet_SaveState_Ru5.46_Rb0.28_sig5.00_b0.0180_a0.0130_L1.0e-06_QDFlag0_SimStep120000_aging_OpenEnd_FSym.mat');
%Sim1 = load('MatchSvet_SaveState_Ru16.39_Rb0.33_sig5.00_b0.0180_a0.0120_L4.0e-07_QDFlag0_SimStep120000_aging_OpenEnd_FSym.mat');
%Sim1 = load('MatchSvet_SaveState_Ru10.93_Rb0.33_sig5.00_b0.0180_a0.0120_L6.0e-07_QDFlag0_SimStep120000_aging_OpenEnd_FSym.mat');
Sim1 = load('MatchSvet_SaveState_Ru6.56_Rb0.30_sig5.00_b0.0200_a0.0140_L1.0e-06_QDFlag0_SimStep60000_aging_OpenEnd_FSym.mat');


% Fast simulation
%Sim2 = load('MatchSvet_SaveState_Ru7.81_Rb0.11_sig5.20_b0.0180_a0.0160_L2.6e-07_QDFlag0_SimStep120000_FSym.mat');
%Sim2 = load('MatchSvet_SaveState_Ru10.14_Rb0.17_sig4.50_b0.0180_a0.0150_L2.6e-07_QDFlag0_SimStep120000_FSym.mat');
%Sim2 = load('MatchSvet_SaveState_Ru10.28_Rb0.11_sig5.00_b0.0180_a0.0160_L1.9e-07_QDFlag0_SimStep120000_FSym.mat');
Sim2 = load('MatchSvet_SaveState_Ru7.81_Rb0.28_sig5.00_b0.0180_a0.0130_L7.0e-07_QDFlag0_SimStep120000_aging_OpenEnd_FSym.mat');


%% import simulation variables

%%%%%%%%%%%%% Parameters that are different between the two simulations 
%Sim1
Sim1.element_location_x1 = Sim1.SaveStateData.InitDistribution.element_location_x1;
Sim1.Fault = Sim1.SaveStateData.Fault;
Sim1.MAX_TIMESTEP = Sim1.SaveStateData.NumOfTimeStep_perSave;
Sim1.t = Sim1.SaveStateData.t;
Sim1.Ru = Sim1.Fault(1).Length/Sim1.Fault(1).hstar;
Sim1.Rb = (Sim1.Fault(1).b - Sim1.Fault(1).a)/Sim1.Fault(1).b;

%Sim2
Sim2.element_location_x1 = Sim2.SaveStateData.InitDistribution.element_location_x1;
Sim2.Fault = Sim2.SaveStateData.Fault;
Sim2.MAX_TIMESTEP = Sim2.SaveStateData.NumOfTimeStep_perSave;
Sim2.t = Sim2.SaveStateData.t;
Sim2.Ru = Sim2.Fault(1).Length/Sim2.Fault(1).hstar;
Sim2.Rb = (Sim2.Fault(1).b - Sim2.Fault(1).a)/Sim2.Fault(1).b;

%%%%%%%%%%%%% Parameters that are shared by the two simulations (Kernel info, or Use Sim1)

% Kernel Info
mfile_knInfo = matfile(['KernelInfo_N',num2str(length(Sim1.element_location_x1)),'W',num2str(Sim1.Fault(1).Length),'_SingleFault.mat']);
KernelInfo = mfile_knInfo.KernelInfo(1,:);
KernelTlength = KernelInfo.MaxEndIndex_AllElem;
MaxEndIndex_AllElem = KernelInfo.MaxEndIndex_AllElem;

% Use Sim1
dt_dyna = Sim1.SaveStateData.GlobalSetup.dt_dyna;
dx = Sim1.SaveStateData.GlobalSetup.dx;
ct = Sim1.SaveStateData.GlobalSetup.ct;
cl = Sim1.SaveStateData.GlobalSetup.cl;
mu = Sim1.SaveStateData.GlobalSetup.mu;

%% PostProcess some output for plots


%% plot whole simulation cyle in timestep, log velocity

%Simulation 1, slow
figure(1)

%only plot half of the simulation
PlotRangeCond = ...
    Sim1.element_location_x1 > Sim1.element_location_x1(length(Sim1.element_location_x1)/2);

imagesc((Sim1.element_location_x1(PlotRangeCond)-Sim1.Fault.Length/2).*1e3,...
    1:Sim1.MAX_TIMESTEP,...
    log10(Sim1.SaveStateData.SlipRate(:,PlotRangeCond)))
%imagesc(element_location_x1,1:600,log10(SaveStateData.SlipRate(1:600,:)))
h=colorbar();
colormap(cycles);
ylabel(h,'log(slip rate), m/s')
ylabel('Time step (integer)')
xlabel('Along fault distance (mm)')
%title([num2str(Sim1.MAX_TIMESTEP),' step simulation, Slow'])
title('Simulation 1, Slow')

pbaspect([0.7 1 1])
xlim([0 200])
xticks(0:50:200)
yticks([0:20000:120000])
title('Sim1')

set(gca,'Ydir','normal')
set(gca,'clim',[-10 1])
set(gca,'Fontsize',20,'Fontweight','bold')
set(gca, 'FontName', 'Helvetica')
set(h,'Fontsize',20,'Fontweight','bold')
set(h, 'FontName', 'Helvetica')
set(gcf, 'Renderer', 'Painters');% make eps clear

%% Plot Timestep length with step Sim1

figure(2)

x = Sim1.t';
y = [1:Sim1.MAX_TIMESTEP]';
plot(x,y,...
     'linewidth',5,'LineStyle','-','Color',DClr.c(1,:));

xlabel('Time (s)')
ylabel('Time step (integer)')
xlim([0 300])
xticks([0:100:1200])
yticks([0:20000:120000])

pbaspect([0.4 1 1])
set(gca,'Fontsize',20,'Fontweight','bold')
set(gca, 'FontName', 'Helvetica')

set(gcf, 'Renderer', 'Painters');% make eps clear


%%
%Simulation 2, fast
figure(3)

PlotRangeCond = ...
    Sim2.element_location_x1 > Sim2.element_location_x1(length(Sim2.element_location_x1)/2);

imagesc((Sim2.element_location_x1(PlotRangeCond)-Sim2.Fault.Length/2).*1e3,...
    1:Sim2.MAX_TIMESTEP,...
    log10(Sim2.SaveStateData.SlipRate(:,PlotRangeCond)))
%imagesc(element_location_x1,1:600,log10(SaveStateData.SlipRate(1:600,:)))
h=colorbar();
colormap(cycles);
ylabel(h,'log(slip rate), m/s')
ylabel('Time step (integer)')
xlabel('Along fault distance (mm)')
%title([num2str(Sim2.MAX_TIMESTEP),' step simulation, Fast'])
title('Sim2')



yticks([0:20000:120000])
pbaspect([0.7 1 1])

xlim([0 200])
xticks(0:50:200)

set(gca,'Ydir','normal')
set(gca,'clim',[-10 1])
set(gca,'Fontsize',20,'Fontweight','bold')
set(gca, 'FontName', 'Helvetica')
set(h,'Fontsize',20,'Fontweight','bold')
set(h, 'FontName', 'Helvetica')
set(gcf, 'Renderer', 'Painters');% make eps clear


%% Plot Timestep length with step Sim2
figure(4)

x = Sim2.t';
y = [1:Sim2.MAX_TIMESTEP]';
plot(x,y,...
     'linewidth',5,'LineStyle','-','Color',DClr.c(1,:));

xlabel('Time (s)')
ylabel('Time step (integer)')
xlim([0 300])
xticks([0:100:1200])
yticks([0:20000:120000])

pbaspect([0.4 1 1])

title('Sim2')
set(gca,'Fontsize',20,'Fontweight','bold')
set(gca, 'FontName', 'Helvetica')

set(gcf, 'Renderer', 'Painters');% make eps clear

%% plot whole simulation cyle in timestep, dyna flag

%Simulation 1, slow
figure(5)

%only plot half of the simulation
PlotRangeCond = ...
    Sim1.element_location_x1 > Sim1.element_location_x1(length(Sim1.element_location_x1)/2);

imagesc((Sim1.element_location_x1(PlotRangeCond)-Sim1.Fault.Length/2).*1e3,...
    1:Sim1.MAX_TIMESTEP,...
    log10(Sim1.SaveStateData.SimulationStateFlag(:,PlotRangeCond)))
%imagesc(element_location_x1,1:600,log10(SaveStateData.SlipRate(1:600,:)))
%h=colorbar();
%colormap(cycles);
%ylabel(h,'log(slip rate), m/s')
ylabel('Time step (integer)')
xlabel('Along fault distance (mm)')
%title([num2str(Sim1.MAX_TIMESTEP),' step simulation, Slow'])
title('Simulation 1, Slow')

pbaspect([0.7 1 1])
xlim([0 200])
xticks(0:50:200)
yticks([0:20000:120000])
title('Sim1')

set(gca,'Ydir','normal')
set(gca,'clim',[-10 1])
set(gca,'Fontsize',20,'Fontweight','bold')
set(gca, 'FontName', 'Helvetica')
set(h,'Fontsize',20,'Fontweight','bold')
set(h, 'FontName', 'Helvetica')
set(gcf, 'Renderer', 'Painters');% make eps clear

