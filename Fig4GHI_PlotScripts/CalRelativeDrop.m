% Obtain information for fracture mechanics analysis, such as, rupture
% speed, stress, and slip weakening distance for one event
clear all
close all
clc

addpath('../') 

DClr = load('ColorScheme_MatlabOrder.mat');

%% because I only use one savestate file, I can also save it outside the loop

%load('MatchSvet_SaveState_Ru6.56_Rb0.30_sig5.00_b0.0200_a0.0140_L1.0e-06_QDFlag0_SimStep60000_aging_OpenEnd_FSym.mat');
%load('MatchSvet_SaveState_Ru8.20_Rb0.30_sig5.00_b0.0200_a0.0140_L8.0e-07_QDFlag0_SimStep120000_aging_OpenEnd_FSym.mat')
%load('MatchSvet_SaveState_Ru9.37_Rb0.25_sig5.00_b0.0240_a0.0180_L7.0e-07_QDFlag0_SimStep60000_aging_OpenEnd_FSym.mat')
%load('MatchSvet_SaveState_Ru5.46_Rb0.25_sig5.00_b0.0120_a0.0090_L6.0e-07_QDFlag0_SimStep60000_aging_OpenEnd_FSym.mat');
%load('MatchSvet_SaveState_Ru8.20_Rb0.25_sig5.00_b0.0120_a0.0090_L4.0e-07_QDFlag0_SimStep60000_aging_OpenEnd_FSym.mat');
%load('MatchSvet_SaveState_Ru10.93_Rb0.33_sig5.00_b0.0180_a0.0120_L6.0e-07_QDFlag0_SimStep120000_aging_OpenEnd_FSym.mat')
%load('MatchSvet_SaveState_Ru5.46_Rb0.28_sig5.00_b0.0180_a0.0130_L1.0e-06_QDFlag0_SimStep120000_aging_OpenEnd_FSym.mat')
%load('MatchSvet_SaveState_Ru4.62_Rb0.33_sig5.00_b0.0180_a0.0120_L1.4e-06_QDFlag0_SimStep120000_aging_OpenEnd_FSym.mat');

%load('MatchSvet_SaveState_Ru8.20_Rb0.30_sig5.00_b0.0200_a0.0140_L8.0e-07_QDFlag0_SimStep120000_aging_OpenEnd_FSym.mat');
%load('MatchSvet_SaveState_Ru6.56_Rb0.30_sig5.00_b0.0200_a0.0140_L1.0e-06_QDFlag0_SimStep60000_aging_OpenEnd_FSym.mat');

%load('MatchSvet_SaveState_Ru10.93_Rb0.27_sig5.00_b0.0220_a0.0160_L6.0e-07_QDFlag0_SimStep60000_aging_OpenEnd_FSym.mat')
%load('MatchSvet_SaveState_Ru8.20_Rb0.27_sig5.00_b0.0220_a0.0160_L8.0e-07_QDFlag0_SimStep60000_aging_OpenEnd_FSym.mat')
%load('MatchSvet_SaveState_Ru6.56_Rb0.27_sig5.00_b0.0220_a0.0160_L1.0e-06_QDFlag0_SimStep60000_aging_OpenEnd_FSym.mat')
%load('MatchSvet_SaveState_Ru5.46_Rb0.27_sig5.00_b0.0220_a0.0160_L1.2e-06_QDFlag0_SimStep60000_aging_OpenEnd_FSym.mat');

%load('MatchSvet_SaveState_Ru10.93_Rb0.25_sig5.00_b0.0240_a0.0180_L6.0e-07_QDFlag0_SimStep60000_aging_OpenEnd_FSym.mat');
%load('MatchSvet_SaveState_Ru8.20_Rb0.25_sig5.00_b0.0240_a0.0180_L8.0e-07_QDFlag0_SimStep60000_aging_OpenEnd_FSym.mat');
%load('MatchSvet_SaveState_Ru5.96_Rb0.25_sig5.00_b0.0240_a0.0180_L1.1e-06_QDFlag0_SimStep60000_aging_OpenEnd_FSym.mat');

%load('MatchSvet_SaveState_Ru8.20_Rb0.23_sig5.00_b0.0260_a0.0200_L8.0e-07_QDFlag0_SimStep60000_aging_OpenEnd_FSym.mat');
%load('MatchSvet_SaveState_Ru6.56_Rb0.23_sig5.00_b0.0260_a0.0200_L1.0e-06_QDFlag0_SimStep60000_aging_OpenEnd_FSym.mat');

%load('MatchSvet_SaveState_Ru6.56_Rb0.21_sig5.00_b0.0280_a0.0220_L1.0e-06_QDFlag0_SimStep60000_aging_OpenEnd_FSym.mat')

%load('MatchSvet_SaveState_Ru10.93_Rb0.17_sig5.00_b0.0240_a0.0200_L4.0e-07_QDFlag0_SimStep60000_aging_OpenEnd_FSym.mat');
%load('MatchSvet_SaveState_Ru21.86_Rb0.17_sig5.00_b0.0240_a0.0200_L2.0e-07_QDFlag0_SimStep60000_aging_OpenEnd_FSym.mat');

%load('MatchSvet_SaveState_Ru6.83_Rb0.21_sig5.00_b0.0240_a0.0190_L8.0e-07_QDFlag0_SimStep60000_aging_OpenEnd_FSym.mat');

%load('MatchSvet_SaveState_Ru6.56_Rb0.38_sig5.00_b0.0160_a0.0100_L1.0e-06_QDFlag0_SimStep60000_aging_OpenEnd_FSym.mat');

%load('MatchSvet_SaveState_Ru10.93_Rb0.25_sig5.00_b0.0160_a0.0120_L4.0e-07_QDFlag0_SimStep60000_aging_OpenEnd_FSym.mat');
%load('MatchSvet_SaveState_Ru7.29_Rb0.25_sig5.00_b0.0160_a0.0120_L6.0e-07_QDFlag0_SimStep60000_aging_OpenEnd_FSym.mat');
%load('MatchSvet_SaveState_Ru5.46_Rb0.25_sig5.00_b0.0160_a0.0120_L8.0e-07_QDFlag0_SimStep60000_aging_OpenEnd_FSym.mat');

%load('MatchSvet_SaveState_Ru4.37_Rb0.23_sig5.00_b0.0350_a0.0270_L2.0e-06_QDFlag0_SimStep60000_aging_OpenEnd_FSym.mat')
%load('MatchSvet_SaveState_Ru5.83_Rb0.23_sig5.00_b0.0350_a0.0270_L1.5e-06_QDFlag0_SimStep60000_aging_OpenEnd_FSym.mat')
%load('MatchSvet_SaveState_Ru8.74_Rb0.23_sig5.00_b0.0350_a0.0270_L1.0e-06_QDFlag0_SimStep60000_aging_OpenEnd_FSym.mat')

%load('MatchSvet_SaveState_Ru9.84_Rb0.25_sig5.00_b0.0360_a0.0270_L1.0e-06_QDFlag0_SimStep60000_aging_OpenEnd_FSym.mat')

load('MatchSvet_SaveState_Ru6.56_Rb0.12_sig5.00_b0.0240_a0.0210_L5.0e-07_QDFlag0_SimStep60000_aging_OpenEnd_FSym.mat');

%% import some variables

element_location_x1 = SaveStateData.InitDistribution.element_location_x1;
Fault = SaveStateData.Fault;

mfile_knInfo = matfile(['KernelInfo_N',num2str(length(element_location_x1)),'W',num2str(Fault(1).Length),'_SingleFault.mat']);
KernelInfo = mfile_knInfo.KernelInfo(1,:);
KernelTlength = KernelInfo.MaxEndIndex_AllElem;
MaxEndIndex_AllElem = KernelInfo.MaxEndIndex_AllElem;

MAX_TIMESTEP = SaveStateData.NumOfTimeStep_perSave;
dt_dyna = SaveStateData.GlobalSetup.dt_dyna;
dx = SaveStateData.GlobalSetup.dx;
ct = SaveStateData.GlobalSetup.ct;
cl = SaveStateData.GlobalSetup.cl;
mu = SaveStateData.GlobalSetup.mu;


t = SaveStateData.t;
Ru = Fault(1).Length/Fault(1).hstar;
Rb = (Fault(1).b - Fault(1).a)/Fault(1).b;

%% Calculate Global Normalized Real Area of Contact

AScaled = (SaveStateData.StateVariable.*Fault(1).V0./Fault(1).Dc).^(Fault(1).b./Fault(1).mu0);
AScaled = AScaled./max(max(AScaled));

%% plot cycle in timestep, log velocity

figure(51)

%only plot half of the simulation
PlotRangeCond = ...
    element_location_x1 > element_location_x1(length(element_location_x1)/2);

imagesc((element_location_x1(PlotRangeCond)-Fault.Length/2).*1e3,...
    1:MAX_TIMESTEP,...
    log10(SaveStateData.SlipRate(:,PlotRangeCond)))
h=colorbar();
colormap(cycles);
ylabel(h,'log(slip rate), m/s')
ylabel('Time step (integer)')
xlabel('Along fault distance (mm)')
%title([num2str(Sim1.MAX_TIMESTEP),' step simulation, Slow'])
title('Simulation 1, Slow')

axis tight
xlim([0 200])
xticks(0:50:200)
title('Sim1')

set(gca,'Ydir','normal')
set(gca,'clim',[-10 1])
set(gca,'Fontsize',20,'Fontweight','bold')
set(gca, 'FontName', 'Helvetica')
set(h,'Fontsize',20,'Fontweight','bold')
set(h, 'FontName', 'Helvetica')
set(gcf, 'Renderer', 'Painters');% make eps clear

%% one event, contact area

% make a mesh grid for space time

[x1_oridata_mat,t_oridata_mat] = meshgrid(element_location_x1,t);
AScaled_oridata_mat = AScaled;
ShearStress_oridata_mat = SaveStateData.ShearStress; 
Slip_oridata_mat = SaveStateData.Slip;
SlipRate_oridata_mat = SaveStateData.SlipRate;
% [x1_oridata_mat,t_oridata_mat] = meshgrid(element_location_x1,t(1:1800));
% AScaled_oridata_mat = AScaled(1:1800,:);
% ShearStress_oridata_mat = SaveStateData.ShearStress(1:1800,:); 
% Slip_oridata_mat = SaveStateData.Slip(1:1800,:);

% interpolate


%RupEndInd = 96345;
%RupEndInd = 82611;
%RupEndInd = 9900;
%RupEndInd = 108531;
%RupEndInd = 46476;
%RupEndInd = 58328;
RupEndInd = 55923;

interp_tstart = t(RupEndInd)-48000*dt_dyna;
interp_tend = t(RupEndInd);
interp_tarr = linspace(interp_tstart,interp_tend,20000);
[x1_interp_mat,t_interp_mat] = meshgrid(element_location_x1,interp_tarr);

plot_starttime = 0.65*(interp_tend-interp_tstart) + interp_tstart;
plot_endtime = 1*(interp_tend-interp_tstart) + interp_tstart;

% speed sign setting
CR_sign_Loc = 110;%ms
%CR_sign_Time = 1.2;%ms
CR_sign_Time = (plot_endtime-plot_starttime)*1e3*0.3;%ms
VrSignLength = 20;
VrtextFontSize = 15;


AScaled_interp_mat...
    = interp2(x1_oridata_mat,t_oridata_mat,AScaled_oridata_mat,x1_interp_mat,t_interp_mat);
ShearStress_interp_mat...
    = interp2(x1_oridata_mat,t_oridata_mat,ShearStress_oridata_mat,x1_interp_mat,t_interp_mat);
Slip_interp_mat...
    = interp2(x1_oridata_mat,t_oridata_mat,Slip_oridata_mat,x1_interp_mat,t_interp_mat);
SlipRate_interp_mat...
    = interp2(x1_oridata_mat,t_oridata_mat,SlipRate_oridata_mat,x1_interp_mat,t_interp_mat);

figure(81)
%imagesc(element_location_x1,interp_tarr-interp_tarr(1),AScaled_interp_mat);
%imagesc(element_location_x1,interp_tarr-interp_tarr(1),AScaled_interp_mat./max(max(AScaled_interp_mat)));
%imagesc(element_location_x1*1e3,(interp_tarr-interp_tarr(1))*1e3,AScaled_interp_mat./max(max(AScaled_interp_mat)));

% use interp_start time as t=0
%imagesc((element_location_x1(PlotRangeCond)-Fault.Length/2)*1e3,(interp_tarr-interp_tarr(1))*1e3,AScaled_interp_mat(:,PlotRangeCond)./AScaled_interp_mat(1,PlotRangeCond));
% use plotstart time as t=0
imagesc((element_location_x1(PlotRangeCond)-Fault.Length/2)*1e3,(interp_tarr-plot_starttime)*1e3,...
    AScaled_interp_mat(:,PlotRangeCond)./AScaled_interp_mat(1,PlotRangeCond));


hold on
%draw Rayleigh wave speed
start_point_R = [CR_sign_Loc CR_sign_Time];
cR = 0.92*ct;
LineLx = VrSignLength;
end_point_R = [start_point_R(1)+LineLx start_point_R(2)+LineLx/cR];

plot([start_point_R(1) end_point_R(1)],[start_point_R(2) end_point_R(2)],'-',...
'linewidth',5,'Color','k')
text(end_point_R(1),end_point_R(2)-0.001,' C_R','Color','k','FontSize',VrtextFontSize)

%draw 0.1 Rayleigh wave speed
start_point_01R = [CR_sign_Loc CR_sign_Time];
cR = 0.92*ct;
LineLx = VrSignLength;
end_point_01R = [start_point_R(1)+LineLx start_point_R(2)+LineLx/(0.1*cR)];

plot([start_point_01R(1) end_point_01R(1)],[start_point_01R(2) end_point_01R(2)],'-',...
'linewidth',5,'Color','k')
text(end_point_01R(1),end_point_01R(2)-0.001,' 0.1C_R','Color','k','FontSize',VrtextFontSize)


hold off

axis tight
xlim([0 200])
xticks(0:50:200)

%y plot range
ylim([0 plot_endtime-plot_starttime]*1e3)

set(gca,'Ydir','normal')
h = colorbar;
colormap('jet')
set(gca,'clim',[0.5 1.1])
xlabel('Along fault distance (mm)')
ylabel('Time (ms)')
title('First event, contact area, FD')
ylabel(h,'contact area, normalized')

set(gca,'Fontsize',20,'Fontweight','bold')
set(gca, 'FontName', 'Helvetica')
set(h,'Fontsize',20,'Fontweight','bold')
set(h, 'FontName', 'Helvetica')
set(gcf, 'Renderer', 'Painters');% make eps clear


%% Calculate Drop

MeanStartInd = 1;
MeanStart_CorrespondT = (interp_tarr(MeanStartInd) - plot_starttime)*1e3;
%MeanFinalInd = 17800;
MeanFinalInd = 20000;
MeanFinal_CorrespondT = (interp_tarr(MeanFinalInd) - plot_starttime)*1e3;


MeanAr_initial = mean(AScaled_interp_mat(MeanStartInd,PlotRangeCond)./AScaled_interp_mat(1,PlotRangeCond));
MeanAr_final = mean(AScaled_interp_mat(MeanFinalInd,PlotRangeCond)./AScaled_interp_mat(1,PlotRangeCond));
% MeanAr_initial = mean(AScaled_interp_mat(MeanStartInd,PlotRangeCond));
% MeanAr_final = mean(AScaled_interp_mat(MeanFinalInd,PlotRangeCond));

ArDrop_perc = (MeanAr_initial - MeanAr_final)/MeanAr_initial

MeanTau_initial = mean(ShearStress_interp_mat(MeanStartInd,PlotRangeCond));
MeanTau_final = mean(ShearStress_interp_mat(MeanFinalInd,PlotRangeCond));

TauDrop_perc = (MeanTau_initial-MeanTau_final)/MeanTau_initial