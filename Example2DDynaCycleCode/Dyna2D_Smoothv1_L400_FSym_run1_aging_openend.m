% 2D dynamic rupture code for slip weakening law, a first very simple model

% 07-07-2023: make an open end model

clear all
close all
clc

%% material setup
ct = 1361;
%cl = 2680;%Plane Strain
cl = 2345;%Plane Stress
rho = 1170;
mu = rho*ct*ct;
lambda = rho*cl*cl-2*mu;
nu = (cl^2-2*ct^2)/(2*(cl^2-ct^2));
Simulation_TimeStep = 60000;

ForceQuasiDynamic_Flag = 0;

savedata_flag = 1;
ForceSysmetric_Flag = 1;

OpenEnd_Flag = 1;

%% Discretisaztion setup, need to be consistent with kernel calculation
dx = 1e-3; % from G&F 2021
discretization_size = dx;
delta_s = dx;
CFL = 0.45;
dt_dyna = CFL*dx/cl;
delta_t = dt_dyna;

%% save global model setup
GlobalSetup.ct = ct;
GlobalSetup.cl = cl;
GlobalSetup.rho = rho;
GlobalSetup.mu = mu;
GlobalSetup.lambda = lambda;
GlobalSetup.nu = nu;
GlobalSetup.Simulation_TimeStep = Simulation_TimeStep;
GlobalSetup.dx = dx;
GlobalSetup.CFL = CFL;
GlobalSetup.dt_dyna = dt_dyna;
GlobalSetup.delta_t = delta_t;
GlobalSetup.delta_s = delta_s;

%% Fault geometry setup a simple step over, need to be consistent with kernel calculation
% making geometry, all fault along x3
%fix parameter

% load geometry from the setup of Kernel calculation
load('Fault_N400W0.4_SingleFault.mat')

% Fault 1 frictional and stress setting
ind = 1;
ElementOnesArray = Fault(ind).Element_No./Fault(ind).Element_No;
%Fault(ind).a = 0.014; %MPa
Fault(ind).a = 0.018; %MPa
%Fault(ind).a = 0.01; %MPa
Fault(ind).b = 0.024;
Fault(ind).mu0 = 0.45;
%Fault(ind).Dc = 1e-5; %m
%Fault(ind).sig0 = 10e6;
%Fault(ind).Dc = 1e-6; %m
Fault(ind).Dc = 0.4e-6; %m
Fault(ind).sig0 = 5.0e6;
%Fault(ind).sig0 = 2.5e6;
%Fault(ind).sig0 = 45e6;
Fault(ind).V0 = 1e-6;%m/s
Fault(ind).Theta0 = Fault(ind).Dc/Fault(ind).V0;
%Fault(ind).tau0 = (Fault(ind).mu0+0.1)*Fault(ind).sig0; 
Fault(ind).tau0 = (Fault(ind).mu0)*Fault(ind).sig0; 
%Fault(ind).Theta_init = 10^9;
Fault(ind).Theta_init = Fault(ind).Dc/Fault(ind).V0;

%% set up heterogeneous stressing rate

% stressing rate
%gaussian pulse
% StressingPulseWidth = 60*dx;
% WithinPulseRelAmp = 1;
% OutsidePulseRelAmp = 3;
% ElementStressingShape = (WithinPulseRelAmp-OutsidePulseRelAmp)*exp(-0.5*((Fault(ind).Element_x1-Fault(ind).Element_x1(1))./StressingPulseWidth).^2)+OutsidePulseRelAmp;

%linear
%StressingPulseWidth = 60*dx;
% LeftEndRelAmp = 1;
% RightEndRelAmp = 1;
% ElementStressingShape = (RightEndRelAmp - LeftEndRelAmp)*(Fault(ind).Element_x1-Fault(ind).Element_x1(1))/(Fault(ind).Element_x1(end)-Fault(ind).Element_x1(1))+LeftEndRelAmp;

%1/r^n decay
% StressingPulseWidth=20*dx;
% StressingDecayWidth=20*dx;
% OutsidePulseRelAmp = 0.3;
% ElementStressingShape = ElementOnesArray;
% InsidePulseCond = Fault(ind).Element_x1-Fault(ind).Element_x1(1) < StressingPulseWidth;
% OutsidePulseCond = Fault(ind).Element_x1-Fault(ind).Element_x1(1) >= StressingPulseWidth;
% ElementStressingShape(InsidePulseCond)=1;
% ElementStressingShape(OutsidePulseCond) = ...
%     (StressingPulseWidth./(Fault(ind).Element_x1(OutsidePulseCond)-Fault(ind).Element_x1(1))).^5.*(1-OutsidePulseRelAmp)+OutsidePulseRelAmp;

%ElementStressingShape = 1-((Fault(ind).Element_x1-Fault(ind).Element_x1(1))./Fault(ind).Length).^(2);
%

% piecewise
% ElementStressingShape = ElementOnesArray;
% StressingPulseWidthInd = 30;
% OutsidePulseRelAmp = 0.4;
% ElementStressingShape(Fault(ind).Element_x1-Fault(ind).Element_x1(StressingPulseWidthInd) > 0) = OutsidePulseRelAmp ;
% ElementStressingShape(Fault(ind).Element_x1-Fault(ind).Element_x1(StressingPulseWidthInd) <= 0) = 1;

% uniform 
ElementStressingShape = ElementOnesArray;

figure(5)
plot([1:length(ElementOnesArray)]*dx*1e3,ElementStressingShape,'linewidth',3)
xlabel('Along fault distance (mm)')
ylabel('Stressing rate scale factor')
%ylim([0 1])
set(gca,'Fontsize',20,'Fontweight','bold')
set(gca, 'FontName', 'Helvetica')
set(gcf, 'Renderer', 'Painters');% make eps clear
drawnow

%% set up heterogeneous normal stress, initial stress need to be adjusted accordingly
%gaussian pulse
SigPulseWidth_Left = 60*dx;
SigPulseWidth_Right= 20*dx;
WithinPulse_SigRelAmp_Left = 1;
WithinPulse_SigRelAmp_Right = 1;
OutsidePulse_SigRelAmp_Left = 1;
OutsidePulse_SigRelAmp_Right = 1;
ElementSigShape_Left = (WithinPulse_SigRelAmp_Left-OutsidePulse_SigRelAmp_Left)*exp(-0.5*(abs(Fault(ind).Element_x1-Fault(ind).Element_x1(1))./SigPulseWidth_Left).^2)+OutsidePulse_SigRelAmp_Left;
ElementSigShape_Right = (WithinPulse_SigRelAmp_Right-OutsidePulse_SigRelAmp_Right)*exp(-0.5*(abs(Fault(ind).Element_x1-Fault(ind).Element_x1(end))./SigPulseWidth_Right).^2)+OutsidePulse_SigRelAmp_Right;
ElementSigShape = ElementSigShape_Left .* ElementSigShape_Right;

% uniform 
%ElementSigShape = ElementOnesArray;

% need to adjust initial shear stress too
ElementTau0Shape = ElementSigShape;

figure(6)
plot([1:length(ElementOnesArray)]*dx*1e3,ElementSigShape,'linewidth',3)
xlabel('Along fault distance (mm)')
ylabel('Normal stress scale factor')
%ylim([0 1])
ylim([0 max(ElementSigShape)])
set(gca,'Fontsize',20,'Fontweight','bold')
set(gca, 'FontName', 'Helvetica')
set(gcf, 'Renderer', 'Painters');% make eps clear
drawnow

%%

Fault(ind).Bulk_backslip_rate = 2e-6;%m/s
%Fault(ind).Bulk_backslip_rate = 0;%m/s
Fault(ind).tau_rate = mu*Fault(ind).Bulk_backslip_rate/Fault(ind).Length; 
Fault(ind).sigma_rate = 0; 


% for testing
hstar = 2/pi*(mu/(1-nu)*Fault(ind).Dc)/((Fault(ind).b-Fault(ind).a)*Fault(ind).sig0)

L_b = 2/pi*(mu/(1-nu)*Fault(ind).Dc)/((Fault(ind).b)*Fault(ind).sig0)

hstar_min = 2/pi*(mu/(1-nu)*Fault(ind).Dc)/((Fault(ind).b-Fault(ind).a)*Fault(ind).sig0*max(ElementSigShape))

L_b_min = 2/pi*(mu/(1-nu)*Fault(ind).Dc)/((Fault(ind).b)*Fault(ind).sig0*max(ElementSigShape))

hstar_max = 2/pi*(mu/(1-nu)*Fault(ind).Dc)/((Fault(ind).b-Fault(ind).a)*Fault(ind).sig0*min(ElementSigShape))

L_b_max = 2/pi*(mu/(1-nu)*Fault(ind).Dc)/((Fault(ind).b)*Fault(ind).sig0*min(ElementSigShape))

Fault(ind).hstar=hstar;
Fault(ind).L_b=L_b;

Ru = Fault(1).Length/Fault(1).hstar;
Rb = (Fault(1).b - Fault(1).a)/Fault(1).b;
RuDotRb = Ru*Rb




Fault(ind).Element_a = Fault(ind).a * ElementOnesArray;%
Fault(ind).Element_b = Fault(ind).b * ElementOnesArray;%
Fault(ind).Element_mu0 = Fault(ind).mu0 * ElementOnesArray;%
Fault(ind).Element_Dc = Fault(ind).Dc * ElementOnesArray;%
Fault(ind).Element_sig0 = Fault(ind).sig0 * ElementOnesArray .* ElementSigShape;%
Fault(ind).Element_V0 = Fault(ind).V0 * ElementOnesArray;%
Fault(ind).Element_Theta0 = Fault(ind).Theta0 * ElementOnesArray;%
Fault(ind).Element_tau0 = Fault(ind).tau0 * ElementOnesArray .* ElementTau0Shape;%
Fault(ind).Element_Theta_init = Fault(ind).Theta_init * ElementOnesArray;%
Fault(ind).Element_sigma_rate = Fault(ind).sigma_rate * ElementOnesArray;%

%%%%%%%%%% added on 07-07-2023 to make available open end simulation
if OpenEnd_Flag == 0
    Fault(ind).Element_tau_rate = Fault(ind).tau_rate .* ElementOnesArray .* ElementStressingShape;%
elseif OpenEnd_Flag == 1
    StaticK = load('StaticKernel_PlaneStressW400.mat');
    Fault(ind).Element_tau_rate =  -StaticK.ElementKernel_static_K12 * (Fault(ind).Bulk_backslip_rate*ElementOnesArray);%
else
    error('something wrong');
end


%% set up nucleation zone
% ForSafeRatio = 10;
% RoughNucLength = Fault(1).D0*mu/Fault(1).tau0*ForSafeRatio;
% HalfNucLengthElementNum = ceil(0.5*RoughNucLength/discretization_size);
% NucCentralEleNo = floor(Fault(1).NumOfElement/4);
% Fault(1).Element_mu0(NucCentralEleNo-HalfNucLengthElementNum:NucCentralEleNo+HalfNucLengthElementNum) = 0.39;

%% Making total array for rupture calculation, need to be the same as in the kernel calculation

% element_array
element_location_x1 = [];
element_location_x2 = [];
element_location_n = [];

% Fault(ind).Element_a = Fault(ind).a * ElementOnesArray;%
% Fault(ind).Element_b = Fault(ind).b * ElementOnesArray;%
% Fault(ind).Element_mu0 = Fault(ind).mu0 * ElementOnesArray;%
% Fault(ind).Element_Dc = Fault(ind).Dc * ElementOnesArray;%
% Fault(ind).Element_sig0 = Fault(ind).sig0 * ElementOnesArray;%
% Fault(ind).Element_V0 = Fault(ind).V0 * ElementOnesArray;%
% Fault(ind).Element_Theta0 = Fault(ind).Theta0 * ElementOnesArray;%
% Fault(ind).Element_tau0 = Fault(ind).tau0 * ElementOnesArray;%

element_location_a = [];
element_location_b = [];
element_location_mu0 = [];
element_location_Dc = [];
element_location_sig0 = [];
element_location_V0 = [];
element_location_Theta0 = [];
element_location_tau0 = [];
element_location_Theta_init = [];
element_location_tau_rate = [];
element_location_sigma_rate = [];

for ind = 1:length(Fault)
    %geometry
    element_location_x1 = cat(1,element_location_x1,Fault(ind).Element_x1);
    element_location_x2 = cat(1,element_location_x2,Fault(ind).Element_x2);
    element_location_n = cat(1,element_location_n,Fault(ind).Element_n);
    
    %stress and friction
    element_location_a = cat(1,element_location_a,Fault(ind).Element_a);
    element_location_b = cat(1,element_location_b,Fault(ind).Element_b);
    element_location_mu0 = cat(1,element_location_mu0,Fault(ind).Element_mu0);
    element_location_Dc = cat(1,element_location_Dc,Fault(ind).Element_Dc);
    element_location_sig0 = cat(1,element_location_sig0,Fault(ind).Element_sig0);   
    element_location_V0 = cat(1,element_location_V0,Fault(ind).Element_V0);
    element_location_Theta0 = cat(1,element_location_Theta0,Fault(ind).Element_Theta0);
    element_location_tau0 = cat(1,element_location_tau0,Fault(ind).Element_tau0);
    element_location_Theta_init = cat(1,element_location_Theta_init,Fault(ind).Element_Theta_init);
    
    element_location_tau_rate = cat(1,element_location_tau_rate,Fault(ind).Element_tau_rate);
    element_location_sigma_rate = cat(1,element_location_sigma_rate,Fault(ind).Element_sigma_rate);
end

%% save geometery and initial distribution
InitDistribution.element_location_x1 = element_location_x1;
InitDistribution.element_location_x2 = element_location_x2;
InitDistribution.element_location_n = element_location_n;
InitDistribution.element_location_a = element_location_a;
InitDistribution.element_location_b = element_location_b;
InitDistribution.element_location_mu0 = element_location_mu0;
InitDistribution.element_location_Dc = element_location_Dc;
InitDistribution.element_location_sig0 = element_location_sig0;
InitDistribution.element_location_V0 = element_location_V0;
InitDistribution.element_location_Theta0 = element_location_Theta0;
InitDistribution.element_location_tau0 = element_location_tau0;
InitDistribution.element_location_Theta_init = element_location_Theta_init;

InitDistribution.element_location_tau_rate = element_location_tau_rate;
InitDistribution.element_location_sigma_rate = element_location_sigma_rate;



%% Load dynamic kernel
% for now, load all the kernel at once, for a larger problem later, may need to
% divide the loading in chunks and use "time" to substitude some space
% Note Mar23, 2022: This is a version that not save kernel in a compact
% way. Will use matrix version instead, the following version is more
% intuitively understood.

% mfile_kn = matfile('ElementKernel.mat');
% ElementKernel = mfile_kn.ElementKernel(1,:);
% KernelTlength = length(ElementKernel(1,1).K11(:,1));

% load matrix version kernel

mfile_kn = matfile(['MatrixKernel_N',num2str(length(element_location_x1)),'W',num2str(Fault(1).Length),'_SingleFault.mat']);
Kernel_K11_all = mfile_kn.Kernel_K11_all;
Kernel_K12_all = mfile_kn.Kernel_K12_all;
Kernel_K22_all = mfile_kn.Kernel_K22_all;

%% load other info
mfile_knCut = matfile(['ElementCut_N',num2str(length(element_location_x1)),'W',num2str(Fault(1).Length),'_SingleFault.mat']);
ElementCut = mfile_knCut.ElementCut(1,:);
mfile_knInfo = matfile(['KernelInfo_N',num2str(length(element_location_x1)),'W',num2str(Fault(1).Length),'_SingleFault.mat']);
KernelInfo = mfile_knInfo.KernelInfo(1,:);
KernelTlength = KernelInfo.MaxEndIndex_AllElem;
MaxEndIndex_AllElem = KernelInfo.MaxEndIndex_AllElem;

GlobalSetup.KernelInfo = KernelInfo;


%% Calculate static kernel 
%%%%%%%%%%%%%%%%%%%%%%%%% using segall solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% since it is not big so I just calculate than in everytime
%need to bury the fault very deep, so I set a ref depth at 5e6
RefDepth = 5e6;
% other parameter in Segall code
slip_PS = -1;
opening_PS = 0;
width_PS = delta_s;
dip_PS = 0;

ElementKernel_static_K11 = nan(length(element_location_x1),length(element_location_x2));
ElementKernel_static_K12 = nan(length(element_location_x1),length(element_location_x2));
ElementKernel_static_K22 = nan(length(element_location_x1),length(element_location_x2));
for ind_i = 1:length(element_location_x1)
    for ind_j = 1:length(element_location_x2)
        
        % stress source is ind_j, obs is ind_i
        x1 = element_location_x1(ind_i)-element_location_x1(ind_j);
        x2 = element_location_x2(ind_i)-element_location_x2(ind_j);
        
        [K11s_PS,K12s_PS,K22s_PS] = ...
            computeStressPlaneStrainSegall10...
            (x1,RefDepth+x2,-0.5*width_PS,RefDepth,slip_PS,opening_PS,width_PS,dip_PS,mu,nu);
        
        ElementKernel_static_K11(ind_i,ind_j) = K11s_PS;
        ElementKernel_static_K12(ind_i,ind_j) = K12s_PS;
        ElementKernel_static_K22(ind_i,ind_j) = K22s_PS;
        
    end
end


%for comparison
% S = zeros(length(element_location_x1),length(element_location_x1));
% for i=1:length(element_location_x1)
%     for j=1:i % using the property of symetry (Sij=Sji)
%         S(i,j) =  mu/(2*pi*dx*(1-nu))*S_part(i,j);
%         if i~=j
%             S(j,i) = S(i,j);
%         end
%     end
% end
% K_matrix = S;
% ElementKernel_static_K12 = K_matrix;


%% Rupture process setup

% simulation time
MAX_TIMESTEP = Simulation_TimeStep;
t = nan(1,MAX_TIMESTEP);
Dt_arr = nan(1,MAX_TIMESTEP);

% making rupture or not flag array
element_RupFlag = isnan(element_location_x1);
% making array to temporally save shear stress along faults
ShearStress_now = nan(1,length(element_location_x1));
% making array to temporally save normal stress along faults
NormalStress_now = nan(1,length(element_location_x1));

% making array to temporally save the slip rate at the timestep start along faults
SR_StepStart = nan(1,length(element_location_x1));
% At the middle of time step, for corrector predictor
SR_StepMid = nan(1,length(element_location_x1));
% Note, slip rate should be the constant for through out the time step, the
% differences of "Start" and "Mid" here is only for the use of solving equation

% making array to temporally save slip rate along faults
SR_now = nan(1,length(element_location_x1));
% making array to temporally save slip rate along faults of the previous
% time step
SR_prev = nan(1,length(element_location_x1));

% When say, start, mid, or end, The following variables really literally
% mean them. Because these variable can change within a time step.

% making array to temporally save the external shear stress (no radiation damping component) along faults at the
% beginning of timestep, initial value is the step 0 initial stress
ExtShearStress_StepStart = element_location_tau0';
% including the radiation component
ShearStress_StepStart = element_location_tau0';
% making array to temporally save normal stress along faults at the
% beginning of timestep, initial value is the step 0 initial stress

% positive normal stress is compress here !!!!!!!!!!!!!!!!!!!! 

NormalStress_StepStart = element_location_sig0';
% making array to temporally save the state variable along faults at the
% beginning of timestep, initial value is the step 0 initial stress
Theta_StepStart = element_location_Theta_init';

% Making the arrays at the middle of timestep
ExtShearStress_StepMid = nan(1,length(element_location_x1));
ShearStress_StepMid = nan(1,length(element_location_x1));
NormalStress_StepMid = nan(1,length(element_location_x1));
Theta_StepMid = nan(1,length(element_location_x1));

% Making the arrays at the end of timestep
ExtShearStress_StepEnd = nan(1,length(element_location_x1));
ShearStress_StepEnd = nan(1,length(element_location_x1));
NormalStress_StepEnd = nan(1,length(element_location_x1));
Theta_StepEnd = nan(1,length(element_location_x1));



% making array to temporally save slip along faults
Slip_start = nan(1,length(element_location_x1));

% doesn't need to save slip mid, can be calculated by v*0.5dt
%Slip_mid = nan(1,length(element_location_x1));

% making array to temporally save slip along faults at the previous time
% step
Slip_end = nan(1,length(element_location_x1));
% ruptured or not flag for slip weakening law
% RFlag = ~isnan(Slip_start);



%% Set up quasidynamic simulation parameter
SimulationStateFlag = 0;
% state == 0, quasi-dyanmic
% state == 1, fully-dynamic

IfInDynaTransitionFlag = 0;% If in the transition from dynamic to quasidynamic
%TransitionStepAfterDynamic_Thre = 250;
TransitionStepAfterDynamic_Thre = 0;
TransitionStepAfterDynamic_count = 0;

%simulation parameters
%delta_d = 0.1*max(element_location_Dc);%displacement step
%step_control_coef = 0.1;

step_control_coef = 0.1;
% control the time length of each step.
%dt = step_control_coef*Dc/[max slip rate]

FricOnly_threRatio_forTransition = 0.001;

FricOnly_threRatio = 1e-5; % if RD stress for FricOnly is lower than external stress by this ratio, RD is ignored
RDOnly_threRatio = 1e-5; % if Fric stress for RDOnly is lower than external stress by this ratio, Fric stress is ignored

GradientMethodExit_DiffThreRatio = 1e-8;
iter_Threhold_GradientMethod = 1e6; % number of iter allowed in one search

%PredCorr_ThreRatio = 1e-5;
PredCorr_ThreRatio = 1e-4;
iter_Threhold_PredCorr = 1e4;


% June-2022, threshold of stress perturbation diff between dynamic and static for getting out of the dynamic phase
Relative_DynaStaticDiff_thre = 1e-3;

%% Initialize save state matrix

% except for slip rate that is constant through out a time-step, other
% parameter is recorded at the end of the time step
%NumOfTimeStep_perSave = 1000;
NumOfTimeStep_perSave = length(t);
SaveStateData.NumOfTimeStep_perSave = NumOfTimeStep_perSave;
SaveStateData.StartStep = 1;
SaveStateData.EndStep = SaveStateData.StartStep+NumOfTimeStep_perSave-1;
SaveStateData.SlipRate = nan(NumOfTimeStep_perSave,length(element_location_x1));
SaveStateData.Slip = nan(NumOfTimeStep_perSave,length(element_location_x1));
SaveStateData.ExtShearStress = nan(NumOfTimeStep_perSave,length(element_location_x1));
SaveStateData.ShearStress = nan(NumOfTimeStep_perSave,length(element_location_x1));
SaveStateData.NormalStress = nan(NumOfTimeStep_perSave,length(element_location_x1));
SaveStateData.StateVariable = nan(NumOfTimeStep_perSave,length(element_location_x1));
SaveStateData.dt_now = nan(NumOfTimeStep_perSave,length(element_location_x1));

SaveStateData.SimulationStateFlag = nan(NumOfTimeStep_perSave,length(element_location_x1));

%% Save setup in savestate (2022.11.5)
SaveStateData.GlobalSetup = GlobalSetup;
SaveStateData.Fault = Fault;
SaveStateData.InitDistribution = InitDistribution;



%% Allocating Recent Slip Rate Profile matrix for the purpose of calculating dynamic stress
% Will allocate again every time change from quasi-dynamic simulation to
% dynamic simulation
RecentSRProfile = zeros(MaxEndIndex_AllElem,length(element_location_x1));%starting from recent to old

RecentSlipProfile_1StepEarlier = zeros(MaxEndIndex_AllElem,length(element_location_x1));
%starting from recent to old, one timestep offset (earlier than)RecentSRProfile
% The use of is to obtain StaticSlip_BeforeRecent;

StaticSlip_BeforeRecent = nan(1,length(element_location_x1));
StaticSlip_BeforeEnterDynamic = nan(1,length(element_location_x1));
StaticShearStress_BeforeEnterDynamic = element_location_tau0';
StaticNormalStress_BeforeEnterDynamic = element_location_sig0';
RecentFullFlag = 0;
%static slip array, only used when RecentSRProfile is full

%Step since the beginning of dynamic calculation, starting from 1
StepSinceDynaSim = 1;

% Record the total slip happens in the dynamic simulation
DynaSlip_start = nan(1,length(element_location_x1));
DynaSlip_end = nan(1,length(element_location_x1));

%% Calculate rupture process matrix multiplation
% a note here, may need to rethink the time in the future, Mar22-2022
% Here is what I think now,
% the time in the process start from 0
% while when calculate kernel, I didn't calculate the first step, because
% it is 0 everywhere except the collocation point, which is \mu/2beta and
% whould actually go into the solving process.
Slip_start(:)=0;
Slip_end(:)=0;

ElementKernel_static_Kshear = ElementKernel_static_K12;
ElementKernel_static_Knormal = ElementKernel_static_K22;

% prepare stress kernel matrix for checking if needed to get into dynamic 07/03/2022

% Kernel_K12_MaxDyna = ...
%     max(reshape(Kernel_K12_all,[length(element_location_x1) MaxEndIndex_AllElem length(element_location_x1)]),[],2);
% Kernel_K12_MaxDyna = reshape(Kernel_K12_MaxDyna,[length(element_location_x1) length(element_location_x1)]);

% Kernel_K11_all = ;
% Kernel_K12_all =;
% Kernel_K13_all =

% start the clock
tic

SimulationStateFlag = 0;% start from quasi-dyanmic rupture

% Flags for getting into the dynamic phase
GoToDyna_Flag = 0;%could be 0 or 1;
%GoOutFromDyna_Flag = 0;%could be 0 or 1;%only use if in the dynamic phase


% array variable for debug
norm_Diff_OldNewSolveV = nan(1,length(t));
norm_Diff_OldNewFindMidV_ExtShear = nan(1,length(t));
norm_Diff_OldNewFindMidV_SR = nan(1,length(t));
norm_Diff_OldNewFindMidV_Theta = nan(1,length(t));

% save filename
%filename_SS = ['SaveStateData_S',num2str(SaveStateData.StartStep),'_E',num2str(SaveStateData.EndStep),'_FullyDynamic.mat'];
filename_SS = ['SaveStateData_S',num2str(SaveStateData.StartStep),'_E',num2str(SaveStateData.EndStep),...
    '_AssymLoadingDc',num2str(Fault(1).Dc),'.mat'];

%For degug
count_DynaOut = 0;
clear TestShearStress_Diff_max TestNormalStress_Diff_max

% set minimum time step for quasi-dynamic simulation
max_timestep_QS = max(element_location_Dc'./Fault(1).Bulk_backslip_rate);

%% Create kernel for estimate wave intensity, for GoToDyna_flag evaluation (2022.11.07)
% version 1 calculate the dynamic stress perturbation of the next
% MaxEndIndex_AllElem time steps (hypothetical)

%length(element_location_x1)*MaxEndIndex_AllElem*length(element_location_x1);
% for me or other people to better understand what is done here, I do the
% matrix manipulation step by step
% first go to linear, , Space*(Time*Space) ---> (Space*Time*Space)*1
% GoToDyna_Test_K11 = reshape(Kernel_K11_all,[length(element_location_x1)*MaxEndIndex_AllElem*length(element_location_x1) 1]);
% GoToDyna_Test_K12 = reshape(Kernel_K12_all,[length(element_location_x1)*MaxEndIndex_AllElem*length(element_location_x1) 1]);
% GoToDyna_Test_K22 = reshape(Kernel_K22_all,[length(element_location_x1)*MaxEndIndex_AllElem*length(element_location_x1) 1]);
% 
% % then (Space*Time*Space)*1 ---> (Space*Time)*Space
% GoToDyna_Test_K11 = reshape(GoToDyna_Test_K11,[length(element_location_x1)*MaxEndIndex_AllElem length(element_location_x1)]);
% GoToDyna_Test_K12 = reshape(GoToDyna_Test_K12,[length(element_location_x1)*MaxEndIndex_AllElem length(element_location_x1)]);
% GoToDyna_Test_K22 = reshape(GoToDyna_Test_K22,[length(element_location_x1)*MaxEndIndex_AllElem length(element_location_x1)]);
% 
% % for calculation (multipulate with diag(SR_Prev)), needed to be transpose
% % now Space*(Space*Time)
% GoToDyna_Test_K11 = GoToDyna_Test_K11';
% GoToDyna_Test_K12 = GoToDyna_Test_K12';
% GoToDyna_Test_K22 = GoToDyna_Test_K22';

%DynaRatio_Thre = 0.1;

% version 2 calculate constant slip rate stress kernel for MaxEndIndex_AllElem time steps (hypothetical)
% for me or other people to better understand what is done here, I do the
% matrix manipulation step by step
% first go to 3d, , Space*(Time*Space) ---> Space*Time*Space
ConstantSR_K11 = reshape(Kernel_K11_all,[length(element_location_x1) MaxEndIndex_AllElem length(element_location_x1)]);
ConstantSR_K12 = reshape(Kernel_K12_all,[length(element_location_x1) MaxEndIndex_AllElem length(element_location_x1)]);
ConstantSR_K22 = reshape(Kernel_K22_all,[length(element_location_x1) MaxEndIndex_AllElem length(element_location_x1)]);

ConstantSR_K11 = reshape(sum(ConstantSR_K11,2),[length(element_location_x1) length(element_location_x1)]);
ConstantSR_K12 = reshape(sum(ConstantSR_K12,2),[length(element_location_x1) length(element_location_x1)]);
ConstantSR_K22 = reshape(sum(ConstantSR_K22,2),[length(element_location_x1) length(element_location_x1)]);

% Nov-2022, threshold of hypothetical stress perturbation, determine
% whether go into dynamic simulation
Hypo_DynaStaticDiff_Thre = 1e-3;



%%
for k=1:length(t)
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Section for debug %%%%%%%%%%%%%%%%%%%%%%%%%%%
    if k==518
        Max_ExtShearStress_StepStart_now = max(ExtShearStress_StepStart);
    end 
    
    
    % check if this step needs to be "quasi-dynamic" or fully dynamic"
    % modify the flag based on the flag at last step (quasidynamic or
    % dynamic) and some transition condition (whether the flag need to change at this step)
    
    % check if the SimulationStateFlag need to be changed from this step
    
    % calculate the slip rate at the start of this time step using the quasi-static assumption, see if the RD
    % factor is significantly large, if so (not), go to (out of) the
    % dynamic simulation.
        
    if k==1
        [SR_Temp,~] = SolveForV_QuasiStatic(ExtShearStress_StepStart,Theta_StepStart,...
                mu,ct,element_location_a',element_location_b',...
                element_location_mu0',element_location_V0',...
                element_location_Theta0',NormalStress_StepStart,...
                FricOnly_threRatio_forTransition,RDOnly_threRatio);

        Approx_dt = step_control_coef*min(element_location_Dc'./SR_Temp);
        SR_prev = SR_Temp; %assign SR_prev a value
        dt_prev = Approx_dt;
        MinDt = Approx_dt;
        DeltaSlip_LastStep = SR_prev*dt_prev;
    else
        if SimulationStateFlag == 0
            Approx_dt = Dt_arr(k-1);
        else %SimulationStateFlag == 1
            Approx_dt = step_control_coef*min(element_location_Dc'./SR_prev);
        end
    end
    
    Approx_dt_arr(k) = Approx_dt;

%% Modified on Nov 7, 2022. GoToDyna_Flag is only needed if SimulationStateFlag == 0

    if SimulationStateFlag == 0
        %%%%%%%%%%%%%%  try using approx_dt to determine GooDyna_Flag
    %     if Approx_dt < 10*dt_dyna
    %         GoToDyna_Flag = 1;
    %     else
    %         GoToDyna_Flag = 0;
    %     end
        
        %%%%%%%%%%%%%%  try using Minimum dt in the last ode45 (for dyna, use dt_dyna) to determine GooDyna_Flag
%         if MinDt < 100*dt_dyna
%             GoToDyna_Flag = 1;
%         else
%             GoToDyna_Flag = 0;
%         end
        
        %%%%%%%%%%%%%% try using criteria between static and radiation damping to determine
        % whether go to dyna state or not
        %DeltaSlip_LastStep = Slip_end - Slip_start;
    %     
    %     RD_stress_estimate = -mu/(2*ct)*SR_prev;
    %     ShearStatic_stress_estimate = (ElementKernel_static_Kshear*DeltaSlip_LastStep')';
    %     %NormalStatic_stress_estimate = ElementKernel_static_Knormal*DeltaSlip_LastStep;
    %     
    %     RD_Ratio_max = max(abs(RD_stress_estimate./ShearStatic_stress_estimate));
    %     if RD_Ratio_max > 0.5
    %         GoToDyna_Flag = 1;
    %     else
    %         GoToDyna_Flag = 0;
    %     end
 
        %%%%%%%%%%%%%% try using criteria between quasidynamic and fully-dynamic to determine
        % whether go to dyna state or not Version 1
        % projection for one dt_dyna
        % compare these two stress field
        % 1. SR_prev for dt_prev time step, causing dynamic stress pertubation
        % 2. Static stress left from last step
%         DynamicStressField_estimate_K11 = sum(diag(SR_prev)*GoToDyna_Test_K11,1);
%         DynamicStressField_estimate_K11 = reshape(DynamicStressField_estimate_K11,[length(element_location_x1) MaxEndIndex_AllElem]);
% 
%         DynamicStressField_estimate_K12 = sum(diag(SR_prev)*GoToDyna_Test_K12,1);
%         DynamicStressField_estimate_K12 = reshape(DynamicStressField_estimate_K12,[length(element_location_x1) MaxEndIndex_AllElem]);  
% 
%         DynamicStressField_estimate_K22 = sum(diag(SR_prev)*GoToDyna_Test_K22,1);
%         DynamicStressField_estimate_K22 = reshape(DynamicStressField_estimate_K22,[length(element_location_x1) MaxEndIndex_AllElem]);        
% 
%         Max_DynamicStressField_estimate_K11 = max(DynamicStressField_estimate_K11,2);
%         Max_DynamicStressField_estimate_K12 = max(DynamicStressField_estimate_K12,2);
%         Max_DynamicStressField_estimate_K22 = max(DynamicStressField_estimate_K22,2);
% 
%         TestShearStress_LeftStatic = ExtShearStress_StepStart';
%         TestNormalStress_LeftStatic = NormalStress_StepStart';
%         
%         DynaRatio_ShearStress = abs(Max_DynamicStressField_estimate_K12)./abs(TestShearStress_LeftStatic);
%         DynaRatio_NormalStress = abs(Max_DynamicStressField_estimate_K22)./abs(TestNormalStress_LeftStatic);
% 
%         MaxDynaRatio = max(max(DynaRatio_ShearStress),max(DynaRatio_NormalStress));
%         
% 
%          if MaxDynaRatio > DynaRatio_Thre
%              GoToDyna_Flag = 1;
%          else
%              GoToDyna_Flag = 0;
%          end

        % whether go to dyna state or not Version 2
        %  MaxEndIndex_AllElem dt_dyna
        % compare some combinitions of the following stress field
        % 1. constant SR_prev for MaxEndIndex_AllElem dt_dyna time steps before, causing dynamic stress pertubation increment
        % 2. constant SR_prev for MaxEndIndex_AllElem dt_dyna time steps before, causing static stress pertubation increment
        % 3. stress at this moment
         if MinDt < MaxEndIndex_AllElem*dt_dyna
%              HypoSRProfile = ones(MaxEndIndex_AllElem,length(element_location_x1))*diag(SR_prev);
%              HypoSRProfile_linear = reshape(HypoSRProfile,[1 MaxEndIndex_AllElem*length(element_location_x1)]);

             HypoStress_11_dynaIncrement = ConstantSR_K11 * SR_prev';        
             HypoStress_12_dynaIncrement = ConstantSR_K12 * SR_prev';         
             HypoStress_22_dynaIncrement = ConstantSR_K22 * SR_prev'; 

             HypoStress_11_staticIncrement = ...
                ElementKernel_static_K11*(SR_prev'.*MaxEndIndex_AllElem*dt_dyna);
             HypoStress_12_staticIncrement = ...
                ElementKernel_static_K12*(SR_prev'.*MaxEndIndex_AllElem*dt_dyna);
             HypoStress_22_staticIncrement = ...
                ElementKernel_static_K22*(SR_prev'.*MaxEndIndex_AllElem*dt_dyna);   

             DiffIncrement_hypo_11 = HypoStress_11_dynaIncrement - HypoStress_11_staticIncrement;
             DiffIncrement_hypo_12 = HypoStress_12_dynaIncrement - HypoStress_12_staticIncrement;
             DiffIncrement_hypo_22 = HypoStress_22_dynaIncrement - HypoStress_22_staticIncrement;

             RelativeDiff_hypo_11 = abs(DiffIncrement_hypo_11./HypoStress_11_dynaIncrement);
             max_RelativeDiff_hypo_11 = max(RelativeDiff_hypo_11(~isinf(RelativeDiff_hypo_11)));
             RelativeDiff_hypo_12 = abs(DiffIncrement_hypo_12./HypoStress_12_dynaIncrement);
             max_RelativeDiff_hypo_12 = max(RelativeDiff_hypo_12(~isinf(RelativeDiff_hypo_12)));
             RelativeDiff_hypo_22 = abs(DiffIncrement_hypo_12./HypoStress_22_dynaIncrement);
             max_RelativeDiff_hypo_22 = max(RelativeDiff_hypo_22(~isinf(RelativeDiff_hypo_22)));

             max_RelativeDiff_hypo_whole = max([max_RelativeDiff_hypo_11 max_RelativeDiff_hypo_12 max_RelativeDiff_hypo_22]);
             max_RelativeDiff_hypo_whole_arr(k) = max_RelativeDiff_hypo_whole;

             % use stress at this point as demominator
             RelativeDiff_hypoCompNow_Shear = abs(DiffIncrement_hypo_12./ExtShearStress_StepStart');
             RelativeDiff_hypoCompNow_Normal = abs(DiffIncrement_hypo_12./NormalStress_StepStart');
             
             max_RelativeDiff_hypoCompNow_whole = max(max(RelativeDiff_hypoCompNow_Shear),max(RelativeDiff_hypoCompNow_Normal));
             max_RelativeDiff_hypoCompNow_whole_arr(k) = max_RelativeDiff_hypoCompNow_whole;

             %if max_RelativeDiff_hypo_whole > Hypo_DynaStaticDiff_Thre
             if max_RelativeDiff_hypoCompNow_whole > Hypo_DynaStaticDiff_Thre
                 GoToDyna_Flag = 1;
             else
                 GoToDyna_Flag = 0;
             end             
         else
             GoToDyna_Flag = 0;
         end

    end    
    
    
    if ~isreal(Theta_StepStart)
        disp('stop');
    end    
        
        
    if SimulationStateFlag == 0 % check if need to go to dynamic simulation
         if GoToDyna_Flag == 1
             SimulationStateFlag = 1;
             TransitionStepAfterDynamic_count = 0; % although I will also set it to 0 later, but still do it here for safe
             
             % Initiate the dynamic rupture array
             RecentSRProfile = zeros(MaxEndIndex_AllElem,length(element_location_x1));%starting from recent to old

             RecentSlipProfile_1StepEarlier = zeros(MaxEndIndex_AllElem,length(element_location_x1));
             %starting from recent to old, one timestep offset (earlier than)RecentSRProfile
             % The use of is to obtain StaticSlip_BeforeRecent;
             % It does not include the slip before this round of the
             % dynamic simulation
             
             
             StaticSlip_BeforeRecent(:) = 0;
             % StaticSlip_BeforeRecent not including the slip before dynamic simualtion,
             
             % added June-2022
             StaticSlip_BeforeEnterDynamic = Slip_start; 
            
             StaticShearStress_BeforeEnterDynamic = ExtShearStress_StepStart;
             StaticNormalStress_BeforeEnterDynamic = NormalStress_StepStart;
             
             RecentFullFlag = 0;
             
             StepSinceDynaSim = 1;
             
             DynaSlip_start(:) = 0;
             DynaSlip_end(:) = 0;
         end
    elseif SimulationStateFlag == 1
        % new in Jun-2022, using a separate flag to control whether the code need to go out of
        % the "dynamic" simulation, criterion being whether dynamic stress
        % is still different enough from static stress
        TestShearStress_IfStatic = ...
            ElementKernel_static_Kshear*(Slip_start-StaticSlip_BeforeEnterDynamic)' + StaticShearStress_BeforeEnterDynamic';
        TestNormalStress_IfStatic = ...
            -(ElementKernel_static_Knormal*(Slip_start-StaticSlip_BeforeEnterDynamic)'- StaticNormalStress_BeforeEnterDynamic');
        
        TestShearStress_IfDynamic = ExtShearStress_StepStart';
        TestNormalStress_IfDynamic = NormalStress_StepStart';
        
        TestShearStress_Diff = TestShearStress_IfStatic - TestShearStress_IfDynamic;
        TestNormalStress_Diff = TestNormalStress_IfStatic - TestNormalStress_IfDynamic;
        
        count_DynaOut = count_DynaOut + 1;
        TestShearStress_Diff_max(count_DynaOut) = max(abs(TestShearStress_Diff./TestShearStress_IfDynamic));
        TestNormalStress_Diff_max(count_DynaOut) = max(abs(TestNormalStress_Diff./TestNormalStress_IfDynamic));
        
        
        
        GoOutFromDyna_Condition = ...
            max(abs(TestShearStress_Diff./TestShearStress_IfDynamic)) < Relative_DynaStaticDiff_thre & ...
            max(abs(TestNormalStress_Diff./TestNormalStress_IfDynamic)) < Relative_DynaStaticDiff_thre;
        
        % have to stay in the dynamic simulation until all the waves have
        % gone
        %GoOutFromDyna_Condition = GoOutFromDyna_Condition & StepSinceDynaSim > MaxEndIndex_AllElem;
        GoOutFromDyna_Condition = GoOutFromDyna_Condition & StepSinceDynaSim > 100;

         if GoOutFromDyna_Condition
             SimulationStateFlag = 0;
         end
        
%         if GoOutFromDyna_Condition && TransitionStepAfterDynamic_count > TransitionStepAfterDynamic_Thre
%             SimulationStateFlag = 0;
%             TransitionStepAfterDynamic_count = 0;
%         else
%             TransitionStepAfterDynamic_count = TransitionStepAfterDynamic_count + 1;
%         end                 
        
        %Relative_DynaStaticDiff_thre
        
    else
        error('Made mistakes somewhere')
    end
    
    % for debug
    if k==240
        disp('debuging...');
    end
 
    
    if ForceQuasiDynamic_Flag == 1
        SimulationStateFlag = 0;
    end
    % Start simulation of this timestep depends on the quasi flag
    if SimulationStateFlag == 0 % perform quasi-dynami simulation
         % Calculate the slip rate at this time step
         
         if k==1
             SR_StepStart = SolveForV_QuasiDyna_PreStepTrial_Mult(ExtShearStress_StepStart,Theta_StepStart,...
                    mu,ct,element_location_a',element_location_b',...
                    element_location_mu0',element_location_V0',...
                    element_location_Theta0',NormalStress_StepStart,...
                    SR_prev,...
                    FricOnly_threRatio,RDOnly_threRatio,...
                    GradientMethodExit_DiffThreRatio,iter_Threhold_GradientMethod);
         else
             SR_StepStart = SR_prev;
         end

    %%%%%%%%%%%%%%%%%%%% prescribe a sudden theta drop at one point (in quasi-dynami phase)
        if k==-1
%             PerturbationWidth = 20*dx;
%             PerturbationAmplog = -10;
%             PerturbationShape = exp(-0.5*((Fault(ind).Element_x1-Fault(ind).Element_x1(1))./PerturbationWidth).^2); 
%             Theta_StepStart = Theta_StepStart.*10.^(PerturbationShape'.*PerturbationAmplog);
             
            PerturbationAmp = 1e-10;
            PerturbationWidthIndex = 20;
            PerturbationShape = ElementOnesArray;
            PerturbationCondition = Fault(ind).Element_x1-Fault(ind).Element_x1(PerturbationWidthIndex) <= 0;
            Theta_StepStart(PerturbationCondition') = Theta_StepStart(PerturbationCondition').*PerturbationAmp;

            %fix tau, change SR
            %need to recalculate the SR start for ode45
            SR_StepStart = SolveForV_QuasiDyna_PreStepTrial_Mult(ExtShearStress_StepStart,Theta_StepStart,...
                    mu,ct,element_location_a',element_location_b',...
                    element_location_mu0',element_location_V0',...
                    element_location_Theta0',NormalStress_StepStart,...
                    SR_prev,...
                    FricOnly_threRatio,RDOnly_threRatio,...
                    GradientMethodExit_DiffThreRatio,iter_Threhold_GradientMethod);
            
            %fix SR, change tau
            %need to recalculate the tau start for ode45
%             friction_AftPtb=element_location_mu0'...
%                 .*(SR_StepStart ./element_location_V0').^(element_location_a'./element_location_mu0')...
%                 .*(Theta_StepStart./element_location_Theta0').^(element_location_b'./element_location_mu0');
%             ExtShearStress_StepStart = NormalStress_StepStart.*friction_AftPtb + mu/(2*ct)*SR_StepStart;

        end  
                        
         ShearStress_StepStart = ExtShearStress_StepStart - mu/(2*ct)*SR_StepStart; 
         
         
        % Building the y variable for ode45
         
        % find mid step info
    
%         [ExtShearStress_StepEnd,NormalStress_StepEnd,SR_StepEnd,Slip_end,Theta_StepEnd,dt_now] = ...
%             FindVEnd_ode45_Mult(ExtShearStress_StepStart,NormalStress_StepStart,...
%                 SR_StepStart,Slip_start,Theta_StepStart,...
%                 ElementKernel_static_Kshear,ElementKernel_static_Knormal,...
%                 step_control_coef,element_location_Dc',...
%                 mu,ct,element_location_a',element_location_b',...
%                 element_location_tau_rate',element_location_sigma_rate');
        
        dt_now = min(min(abs(step_control_coef*element_location_Dc'./SR_StepStart)),max_timestep_QS);

        [ExtShearStress_StepEnd,NormalStress_StepEnd,SR_StepEnd,Slip_end,Theta_StepEnd,MinDt] = ...
            FindVEnd_ode45_Mult_MinDt(ExtShearStress_StepStart,NormalStress_StepStart,SR_StepStart,Slip_start,Theta_StepStart,...
                ElementKernel_static_Kshear,ElementKernel_static_Knormal,...
                dt_now,element_location_Dc',...
                mu,ct,element_location_a',element_location_b',...
                element_location_mu0',element_location_V0',element_location_Theta0',...
                element_location_tau_rate',element_location_sigma_rate');            
         
            
         ShearStress_StepEnd = ExtShearStress_StepEnd - (mu/(2*ct))*SR_StepEnd;   
%          pos=~isfinite(Theta_StepEnd);
%          Theta_StepEnd(pos)=Theta_StepStart(pos)+dt_now;


        %%% option: forcing symetric results, only usable for planar fault 
        if ForceSysmetric_Flag == 1
            if mod(Fault(1).NumOfElement,2) == 1
                error(['when forcing symmertic results, number of element of the planar fault need ' ...
                    'to be even, now NumOfElement = ',num2str(Fault(1).NumOfElement)]);
            elseif mod(Fault(1).NumOfElement,2) == 0
                
                SR_StepEnd(1:Fault(1).NumOfElement/2) = flip(SR_StepEnd(Fault(1).NumOfElement/2+1:end));
                Slip_end(1:Fault(1).NumOfElement/2) = flip(Slip_end(Fault(1).NumOfElement/2+1:end));
                ExtShearStress_StepEnd(1:Fault(1).NumOfElement/2) = flip(ExtShearStress_StepEnd(Fault(1).NumOfElement/2+1:end));
                ShearStress_StepEnd(1:Fault(1).NumOfElement/2) = flip(ShearStress_StepEnd(Fault(1).NumOfElement/2+1:end));
                NormalStress_StepEnd(1:Fault(1).NumOfElement/2) = flip(NormalStress_StepEnd(Fault(1).NumOfElement/2+1:end));
                Theta_StepEnd(1:Fault(1).NumOfElement/2) = flip(Theta_StepEnd(Fault(1).NumOfElement/2+1:end));

            else
                error('illegal input of ForceSysmetric_Flag')
            end
        end
            
         % saving the step
        if mod(k,NumOfTimeStep_perSave) ~= 0 && k~=length(t)
            SaveStateInd = mod(k,NumOfTimeStep_perSave);
            SaveStateData.SlipRate(SaveStateInd,:) = SR_StepEnd;
            SaveStateData.Slip(SaveStateInd,:) = Slip_end;
            SaveStateData.ExtShearStress(SaveStateInd,:) = ExtShearStress_StepEnd;
            SaveStateData.ShearStress(SaveStateInd,:) = ShearStress_StepEnd;
            SaveStateData.NormalStress(SaveStateInd,:) = NormalStress_StepEnd;
            SaveStateData.StateVariable(SaveStateInd,:) = Theta_StepEnd;
            SaveStateData.dt_now(SaveStateInd,:) = dt_now;
            SaveStateData.SimulationStateFlag(SaveStateInd,:) = SimulationStateFlag;
        elseif mod(k,NumOfTimeStep_perSave) == 0 || k==length(t)
            % First, write this timestep into matrix
            SaveStateInd = NumOfTimeStep_perSave;
            SaveStateData.SlipRate(SaveStateInd,:) = SR_StepEnd;
            SaveStateData.Slip(SaveStateInd,:) = Slip_end;
            SaveStateData.ExtShearStress(SaveStateInd,:) = ExtShearStress_StepEnd;
            SaveStateData.ShearStress(SaveStateInd,:) = ShearStress_StepEnd;
            SaveStateData.NormalStress(SaveStateInd,:) = NormalStress_StepEnd;
            SaveStateData.StateVariable(SaveStateInd,:) = Theta_StepEnd;
            SaveStateData.dt_now(SaveStateInd,:) = dt_now;
            SaveStateData.SimulationStateFlag(SaveStateInd,:) = SimulationStateFlag;

            % Second, save this matrix to file
            %save(filename_SS,'SaveStateData','-v7.3');

            % Lastly, restart start and end timestep
            SaveStateData.StartStep = round(k/NumOfTimeStep_perSave)*NumOfTimeStep_perSave + 1;
            SaveStateData.EndStep = SaveStateData.StartStep+NumOfTimeStep_perSave-1;
        end
        
        if ~isreal(Theta_StepEnd)
            disp('stop');
        end          
        
        
        % Update temporary arrays for the next timestep
        DeltaSlip_LastStep = Slip_end - Slip_start;
        
        SR_prev = SR_StepEnd;
        Slip_start = Slip_end;
        ExtShearStress_StepStart = ExtShearStress_StepEnd;
        Theta_StepStart = Theta_StepEnd;
        NormalStress_StepStart = NormalStress_StepEnd;
        
        dt_prev = dt_now;
        
        %%%%%%%%%% for debug, testing if elastic stress equal friction %%%%%%%%%%%
        elastic_stress_test_now = ShearStress_StepEnd;
        friction_test_now = element_location_mu0'.*(SR_StepEnd./element_location_V0')...
            .^(element_location_a'./element_location_mu0')...
            .*(Theta_StepEnd./element_location_Theta0').^(element_location_b'./element_location_mu0');
        
        Diff_ElasticStress_FaultStrength = ...
            elastic_stress_test_now-friction_test_now.*NormalStress_StepEnd;
        
        maxDiffRatio_forTest = max(abs(Diff_ElasticStress_FaultStrength./elastic_stress_test_now));
        
        
       
        

    elseif SimulationStateFlag == 1 % perform dynami simulation
        
        dt_now = dt_dyna;
        MinDt = dt_now;
        
        % Calculating the external elastic stress profile at this time step
        % assuming no elastic field propagating out an element
        % dynamic component
        RecentSRProfile_linear = reshape(RecentSRProfile,[1 MaxEndIndex_AllElem*length(element_location_x1)]);
        ExtStress_11_dyna_StepEnd = (Kernel_K11_all * RecentSRProfile_linear')';        
        ExtStress_12_dyna_StepEnd = (Kernel_K12_all * RecentSRProfile_linear')';         
        ExtStress_22_dyna_StepEnd = (Kernel_K22_all * RecentSRProfile_linear')';        

        ExtStress_11_StepEnd = ExtStress_11_dyna_StepEnd;
        ExtStress_12_StepEnd = ExtStress_12_dyna_StepEnd;
        ExtStress_22_StepEnd = ExtStress_22_dyna_StepEnd;
        
        %static component
        if RecentFullFlag == 1
        
            StaticSlip_BeforeRecent = RecentSlipProfile_1StepEarlier(end,:);

            Stress_11_static_StepEnd = (ElementKernel_static_K11 * StaticSlip_BeforeRecent')';        
            Stress_12_static_StepEnd = (ElementKernel_static_K12 * StaticSlip_BeforeRecent')';         
            Stress_22_static_StepEnd = (ElementKernel_static_K22 * StaticSlip_BeforeRecent')'; 

            ExtStress_11_StepEnd = ExtStress_11_StepEnd + Stress_11_static_StepEnd;
            ExtStress_12_StepEnd = ExtStress_12_StepEnd + Stress_12_static_StepEnd;
            ExtStress_22_StepEnd = ExtStress_22_StepEnd + Stress_22_static_StepEnd;
        
        end

               
        % for now, since the geometry is simple (normal vector all being [0 1])
        % I just explicitly write out the calculation for shear and normal
        % stress (:)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Note(June-16-2022):It looks a bit confusing but could be right
        % StaticSlip_BeforeRecent and StaticShearStress_BeforeRecent are at
        % different time point, StaticSlip_BeforeRecent is before the
        % Recent Matrix, StaticShearStress_BeforeRecent is at the beginning
        % of the dynamic phase. 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        ExtShearStress_StepEnd = ExtStress_12_StepEnd + StaticShearStress_BeforeEnterDynamic;
        %ExtShearStress_StepEnd = ExtStress_12_StepEnd + element_location_tau0';
        
        
        % notice that in the input, compress is positive, while in the kernel
        % convention, negative is compress;
        % Here, for the convenience of slip rate calculation, NormalStress
        % output use "compress being normal", so there is a minus in front of
        % the expression
        ExtNormalStress_StepEnd = -(ExtStress_22_StepEnd - StaticNormalStress_BeforeEnterDynamic);
        %ExtNormalStress_StepEnd = -(ExtStress_22_StepEnd - element_location_sig0');
        
        % assume that dt is small enough so the Stress_StepEnd is the same
        % through out the time step is constant, and StressStepEnd is the same as
        % StressStepMid or StressStepStart
        % Note!!!! Sylvain suggest to improve this step and not using this
        % assumption by calculating stress kernel at half step.
        
        ExtShearStress_StepStart = ExtShearStress_StepEnd;
        NormalStress_StepStart = NormalStress_StepEnd;       
        ExtShearStress_StepMid = ExtShearStress_StepEnd;
        NormalStress_StepMid = NormalStress_StepEnd;
        
        
        % Solve for slip rate of at the start of this step using theta at the start 
        
%          SR_StepStart = SolveForV_Dyna(ExtShearStress_StepStart,Theta_StepStart,...
%                 mu,ct,element_location_a',element_location_b',...
%                 element_location_mu0',element_location_V0',...
%                 element_location_Theta0',NormalStress_StepStart,...
%                 FricOnly_threRatio,RDOnly_threRatio,...
%                 GradientMethodExit_DiffThreRatio,iter_Threhold_GradientMethod);
            
%          SR_StepStart = SolveForV_Dyna_PreStepTrial(ExtShearStress_StepStart,Theta_StepStart,...
%                 mu,ct,element_location_a',element_location_b',...
%                 element_location_mu0',element_location_V0',...
%                 element_location_Theta0',NormalStress_StepStart,...
%                 SR_prev,...
%                 GradientMethodExit_DiffThreRatio,iter_Threhold_GradientMethod);            
            
         SR_StepStart = SolveForV_Dyna_PreStepTrial_Mult_NonNegV2(ExtShearStress_StepStart,Theta_StepStart,...
                mu,ct,element_location_a',element_location_b',...
                element_location_mu0',element_location_V0',...
                element_location_Theta0',NormalStress_StepStart,...
                SR_prev,...
                GradientMethodExit_DiffThreRatio,iter_Threhold_GradientMethod);             
         ShearStress_StepStart = ExtShearStress_StepStart - mu/(2*ct)*SR_StepStart;
         ShearStress_StepEnd = ShearStress_StepStart;

        % 07/17 smoothing the slip rate
        SR_StepStart = smoothdata(SR_StepStart,'gaussian',3);         
         
         if sum(isinf(SR_StepStart)) ~= 0
             inf_ind = find(isinf(SR_StepStart));
         end
         
         % use predictor-corrector to solve for slip rate of this time
         % step. Need to use that because theta evolution is non linear
         
         %%% 03-29, decide to use explicit method
         
%          [SR_StepMid,Theta_StepMid] = FindVMid_Dyna(ExtShearStress_StepStart,...
%                 NormalStress_StepStart,SR_StepStart,Theta_StepStart,...
%                 dt_dyna,element_location_Dc',mu,ct,element_location_a',element_location_b',...
%                 element_location_mu0',element_location_V0',element_location_Theta0',...
%                 FricOnly_threRatio,RDOnly_threRatio,...
%                 GradientMethodExit_DiffThreRatio,iter_Threhold_GradientMethod,...
%                 PredCorr_ThreRatio,iter_Threhold_PredCorr); 
            
         %SR_now = SR_StepMid;  
         SR_now = SR_StepStart;
         SR_StepEnd = SR_now;
          
         % and use the slip rate to calculate other parameters at the end
         % of the time step.
         Slip_end = Slip_start + SR_now*dt_now;
         
         % keep track of the slip happens in this round of the dynamic simulation
         DynaSlip_end = DynaSlip_start + SR_now*dt_now;
         
         % calculate state variable at the end of the step
%          Theta_StepEnd = ...
%                 element_location_Dc'./SR_now +...
%                 (Theta_StepStart-element_location_Dc'./SR_now)...
%                 .*exp(-SR_now.*dt_now./element_location_Dc'); 
         

% avoid big number minus big number

% Analytical evolution         
         Theta_StepEnd = ...
                (element_location_Dc'./SR_now).*(1 - exp(-SR_now.*dt_now./element_location_Dc')) + ...
                Theta_StepStart.*exp(-SR_now.*dt_now./element_location_Dc'); 

         if ~isreal(Theta_StepEnd)
            disp('stop');
         end 


         if sum(isinf(1./Theta_StepEnd)) ~= 0
             inf_ind = find(isinf(1./Theta_StepEnd));
         end

        %%% option: forcing symetric results, only usable for planar fault 
        if ForceSysmetric_Flag == 1
            if mod(Fault(1).NumOfElement,2) == 1
                error(['when forcing symmertic results, number of element of the planar fault need ' ...
                    'to be even, now NumOfElement = ',num2str(Fault(1).NumOfElement)]);
            elseif mod(Fault(1).NumOfElement,2) == 0
                
                SR_StepEnd(1:Fault(1).NumOfElement/2) = flip(SR_StepEnd(Fault(1).NumOfElement/2+1:end));
                Slip_end(1:Fault(1).NumOfElement/2) = flip(Slip_end(Fault(1).NumOfElement/2+1:end));
                ExtShearStress_StepEnd(1:Fault(1).NumOfElement/2) = flip(ExtShearStress_StepEnd(Fault(1).NumOfElement/2+1:end));
                ShearStress_StepEnd(1:Fault(1).NumOfElement/2) = flip(ShearStress_StepEnd(Fault(1).NumOfElement/2+1:end));
                NormalStress_StepEnd(1:Fault(1).NumOfElement/2) = flip(NormalStress_StepEnd(Fault(1).NumOfElement/2+1:end));
                Theta_StepEnd(1:Fault(1).NumOfElement/2) = flip(Theta_StepEnd(Fault(1).NumOfElement/2+1:end));

            else
                error('illegal input of ForceSysmetric_Flag')
            end
        end         
            
         % saving the step
        if mod(k,NumOfTimeStep_perSave) ~= 0 && k~=length(t)
            SaveStateInd = mod(k,NumOfTimeStep_perSave);
            SaveStateData.SlipRate(SaveStateInd,:) = SR_StepEnd;
            SaveStateData.Slip(SaveStateInd,:) = Slip_end;
            SaveStateData.ExtShearStress(SaveStateInd,:) = ExtShearStress_StepEnd;
            SaveStateData.ShearStress(SaveStateInd,:) = ShearStress_StepEnd;
            SaveStateData.NormalStress(SaveStateInd,:) = NormalStress_StepEnd;
            SaveStateData.StateVariable(SaveStateInd,:) = Theta_StepEnd;
            SaveStateData.dt_now(SaveStateInd,:) = dt_now;
            SaveStateData.SimulationStateFlag(SaveStateInd,:) = SimulationStateFlag;
        elseif mod(k,NumOfTimeStep_perSave) == 0 || k==length(t)
            % First, write this timestep into matrix
            SaveStateInd = NumOfTimeStep_perSave;
            SaveStateData.SlipRate(SaveStateInd,:) = SR_StepEnd;
            SaveStateData.Slip(SaveStateInd,:) = Slip_end;
            SaveStateData.ExtShearStress(SaveStateInd,:) = ExtShearStress_StepEnd;
            SaveStateData.ShearStress(SaveStateInd,:) = ShearStress_StepEnd;
            SaveStateData.NormalStress(SaveStateInd,:) = NormalStress_StepEnd;
            SaveStateData.StateVariable(SaveStateInd,:) = Theta_StepEnd;
            SaveStateData.dt_now(SaveStateInd,:) = dt_now;
            SaveStateData.SimulationStateFlag(SaveStateInd,:) = SimulationStateFlag;

            % Second, save this matrix to file
            %save(filename_SS,'SaveStateData','-v7.3');

            % Lastly, restart start and end timestep
            SaveStateData.StartStep = (round(k/NumOfTimeStep_perSave)-1)*NumOfTimeStep_perSave + 1;
            SaveStateData.EndStep = SaveStateData.StartStep+NumOfTimeStep_perSave-1;
        end            
        
        %for debug
%         if SimulationStateFlag == 0
%             a = SimulationStateFlag;
%         end
        
        
        % Update temporary arrays for the next timestep
        
        % array for dynamic calculation
        RecentSRProfile(2:end,:) = RecentSRProfile(1:end-1,:);
        RecentSRProfile(1,:) = SR_StepEnd;

        RecentSlipProfile_1StepEarlier(2:end,:)=RecentSlipProfile_1StepEarlier(1:end-1,:);
        RecentSlipProfile_1StepEarlier(1,:)= Slip_start-StaticSlip_BeforeEnterDynamic;
    
        if StepSinceDynaSim==MaxEndIndex_AllElem 
            % starting from this time step, Recent SR profile is full, need to consider static stress
            % using Segall solution

            % 
            RecentFullFlag = 1;
        end
        StepSinceDynaSim = StepSinceDynaSim + 1;        
        
        
        % update infomation for next step
        DeltaSlip_LastStep = Slip_end - Slip_start;
        
        SR_prev = SR_StepEnd;
        Slip_start = Slip_end;
        DynaSlip_start = DynaSlip_end;
        ExtShearStress_StepStart = ExtShearStress_StepEnd;
        Theta_StepStart = Theta_StepEnd;
        NormalStress_StepStart = NormalStress_StepEnd;
        
        
        dt_prev = dt_now;
        
    else
        error('Made mistakes somewhere')
    end
    
    % save t info
    if k==1
        t(k) = dt_now;
    else
        t(k) = t(k-1)+dt_now;
    end
    Dt_arr(k) = dt_now;
   
    if mod(k,100)==1
        elapsed_time = toc;
        disp(['Total Elase Time = ',num2str(elapsed_time),' s']);
        disp(['Step = ',num2str(k),', Time = ',num2str(t(k)),' s']);
    end    
    
end

%% because I only use one savestate file, I can also save it outside the loop
% save time array;
SaveStateData.t = t;
hstar
L_b
Ru
Rb
RuDotRb
if OpenEnd_Flag == 1 
    filename_SS_thisExperiment = ['MatchSvet_SaveState_Ru',num2str(Ru,'%.2f'),'_Rb',num2str(Rb,'%.2f'),'_sig',num2str(Fault(1).sig0./1e6,'%.2f'),...
        '_b',num2str(Fault(1).b ,'%.4f'),'_a',num2str(Fault(1).a ,'%.4f'),'_L',num2str(Fault(1).Dc ,'%.1e'),...
        '_QDFlag',num2str(ForceQuasiDynamic_Flag,'%d'),'_SimStep',num2str(Simulation_TimeStep,'%d'),'_aging_OpenEnd'];
else
        filename_SS_thisExperiment = ['MatchSvet_SaveState_Ru',num2str(Ru,'%.2f'),'_Rb',num2str(Rb,'%.2f'),'_sig',num2str(Fault(1).sig0./1e6,'%.2f'),...
        '_b',num2str(Fault(1).b ,'%.4f'),'_a',num2str(Fault(1).a ,'%.4f'),'_L',num2str(Fault(1).Dc ,'%.1e'),...
        '_QDFlag',num2str(ForceQuasiDynamic_Flag,'%d'),'_SimStep',num2str(Simulation_TimeStep,'%d'),'_aging_FixedEnd'];
end

if ForceSysmetric_Flag == 1
    filename_SS_thisExperiment = [filename_SS_thisExperiment,'_FSym','.mat'];
else
    filename_SS_thisExperiment = [filename_SS_thisExperiment,'.mat'];
end

if savedata_flag == 1
    save(filename_SS_thisExperiment,'SaveStateData','-v7.3');
end


%% testing loading time
% tic;
% ind_i = 8;
% Kernel_thisElement = mfile_kn.ElementKernel(1,:);
% elapsed_time_loading = toc

%% visualize rupture process


% figure(2)
% % str_Dc = num2str(Dc);
% % str_sigma = num2str(sigma);
% % str_spacing = num2str(spacing_factor);
% 
% StateFlag = SaveStateData.SimulationStateFlag;
% %log_speed_matrix = log10(SaveStateData.SlipRate);
% %xlabel(['spacing/',str_spacing,'Dc'],'Fontsize',20)
% h=colorbar();
% ylabel(h,'State flag','Fontsize',20)
% ylabel('Time step (integer)')
% %title(['Dc=',str_Dc,', Shear Impedence, Speed'],'Fontsize',20)
% hold on
% %log_speed_matrix = nan(10000,NumElement);
% % for k=1:1:10000
% %     log_speed_matrix(k,:) = log10(TimeStep(k).Speed_start');
% %     set(gca,'Fontsize',15)
% % end
% imagesc(StateFlag)
% %imagesc(log_speed_matrix(1:1888,:))
% hold off
% axis tight
% %set(gca,'clim',[-12 1])


%% plot time series of shear stress average

Import_Svetlizky2019a_fig1b;

FaultAver_ShearStress = mean(SaveStateData.ShearStress,2);
figure(43)

plot(Svetlizky2019aFig1b.Time-Svetlizky2019aFig1b.Time(1),Svetlizky2019aFig1b.Fs_KN,...
    'linewidth',5,'Color',[0.7 0.7 0.7])

hold on

plot(t,FaultAver_ShearStress.*(200*5.5/1e6)./1e3,'linewidth',2,'LineStyle','-')

xlabel('Time (s)')
ylabel('Net Shear Force (kN)')

%xlim([t(20000) t(20000)+(1110-1050)])
ylim([1.8 3.25])
%yticks([2:0.4:3.6])

%legend(p,legend_str)
%grid on
set(gca,'Fontsize',20,'Fontweight','bold')
set(gca, 'FontName', 'Helvetica')
set(gcf, 'Renderer', 'Painters');% make eps clear

hold off

%% plot cycle in timestep, log velocity

figure(51)

imagesc(element_location_x1,1:MAX_TIMESTEP,log10(SaveStateData.SlipRate))
%imagesc(element_location_x1,1:600,log10(SaveStateData.SlipRate(1:600,:)))
h=colorbar();
colormap(cycles);
ylabel(h,'log(slip rate), m/s')
ylabel('Time step (integer)')
xlabel('Along fault distance (m)')
title([num2str(MAX_TIMESTEP),' step simulation, Fully dynamic'])


axis tight
set(gca,'Ydir','normal')
set(gca,'clim',[-10 1])
set(gca,'Fontsize',20,'Fontweight','bold')
set(gca, 'FontName', 'Helvetica')
set(h,'Fontsize',20,'Fontweight','bold')
set(h, 'FontName', 'Helvetica')
set(gcf, 'Renderer', 'Painters');% make eps clear




%% one event, contact area
AScaled = (SaveStateData.StateVariable.*Fault(1).V0./Fault(1).Dc).^(Fault(1).b./Fault(1).mu0);
AScaled = AScaled./max(max(AScaled));

% make a mesh grid for space time

[x1_oridata_mat,t_oridata_mat] = meshgrid(element_location_x1,t);
AScaled_oridata_mat = AScaled;
ShearStress_oridata_mat = SaveStateData.ShearStress; 
Slip_oridata_mat = SaveStateData.Slip;
% [x1_oridata_mat,t_oridata_mat] = meshgrid(element_location_x1,t(1:1800));
% AScaled_oridata_mat = AScaled(1:1800,:);
% ShearStress_oridata_mat = SaveStateData.ShearStress(1:1800,:); 
% Slip_oridata_mat = SaveStateData.Slip(1:1800,:);

% interpolate

interp_tstart = t(39335)-48000*dt_dyna;
interp_tend = t(39335);
interp_tarr = linspace(interp_tstart,interp_tend,20000);
[x1_interp_mat,t_interp_mat] = meshgrid(element_location_x1,interp_tarr);

AScaled_interp_mat...
    = interp2(x1_oridata_mat,t_oridata_mat,AScaled_oridata_mat,x1_interp_mat,t_interp_mat);
ShearStress_interp_mat...
    = interp2(x1_oridata_mat,t_oridata_mat,ShearStress_oridata_mat,x1_interp_mat,t_interp_mat);
Slip_interp_mat...
    = interp2(x1_oridata_mat,t_oridata_mat,Slip_oridata_mat,x1_interp_mat,t_interp_mat);

figure(81)
%imagesc(element_location_x1,interp_tarr-interp_tarr(1),AScaled_interp_mat);
%imagesc(element_location_x1,interp_tarr-interp_tarr(1),AScaled_interp_mat./max(max(AScaled_interp_mat)));
imagesc(element_location_x1*1e3,(interp_tarr-interp_tarr(1))*1e3,AScaled_interp_mat./AScaled_interp_mat(1,:));

plot_timeshift =  0;%ms
CR_sign_Loc = 150;%ms
CR_sign_Time = 7.5;%ms
VrSignLength = 30;
VrtextFontSize = 25;

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
set(gca,'Ydir','normal')
h = colorbar;
colormap('jet')
%set(gca,'clim',[0.5 1.1])
xlabel('Along fault distance (mm)')
ylabel('Time (ms)')
title('First event, contact area, FD')
ylabel(h,'contact area, normalized')

set(gca,'Fontsize',20,'Fontweight','bold')
set(gca, 'FontName', 'Helvetica')
set(h,'Fontsize',20,'Fontweight','bold')
set(h, 'FontName', 'Helvetica')
set(gcf, 'Renderer', 'Painters');% make eps clear









%% time series of real area contact at 150
figure(92)
plot(t,AScaled(:,150),'linewidth',3)
ylabel('contact area, normalized')
xlabel('Time, s')
set(gca,'Fontsize',20,'Fontweight','bold')
set(gca, 'FontName', 'Helvetica')
set(gcf, 'Renderer', 'Painters');% make eps clear

%%

% figure(91);clf
% 
% cla
% % % % % 
% %only plot half of the simulation
% % PlotRangeCond = ...
% %     element_location_x1 > element_location_x1(length(element_location_x1)/2);
% 
% PlotRangeCond =  element_location_x1 > -9999999;
% 
% 
% %location_mat = repmat((element_location_x1(PlotRangeCond)-Fault.Length/2)'.*1e3,MAX_TIMESTEP,1);
% location_mat = repmat(element_location_x1(PlotRangeCond)'.*1e3,MAX_TIMESTEP,1);
% 
% step=1;
% pcolor(location_mat(1:step:end,:),...
%     SaveStateData.Slip(1:step:end,PlotRangeCond)*1e3,...
%     (SaveStateData.SlipRate(1:step:end,PlotRangeCond))), shading flat
% h=colorbar();
% colormap(cycles);
% ylabel(h,'log(slip rate), m/s')
% ylabel('Slip (mm)')
% xlabel('Along fault distance (mm)')
% %title([num2str(Sim1.MAX_TIMESTEP),' step simulation, Slow'])
% title('Simulation 1, Slow')
% 
% axis tight
% xlim([0 200])
% xticks(0:50:200)
% title('Sim1')
% set(gca,'TickDir','Out')
% set(gca,'Ydir','normal')
% %set(gca,'clim',[-10 1])
% set(gca,'Fontsize',20,'Fontweight','bold')
% set(gca, 'FontName', 'Helvetica')
% set(h,'Fontsize',20,'Fontweight','bold')
% set(h, 'FontName', 'Helvetica')
% set(gcf, 'Renderer', 'Painters');% make eps clear