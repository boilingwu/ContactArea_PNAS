% making plots for the note
% plots for the note has number in 100
clear all
clc
%%
data_folder = '/Users/baoning/Dropbox/2DDynamicCycle/GRL_Draft/AprMatchSvet_ForPaper/May20_L400Simulation/W400mmExperiment/OpenEndAging_PStress/';
DClr = load([data_folder,'ColorScheme_MatlabOrder.mat']);

Test_ID = 'PS29';
filename_FuncksK2 = ['ForTestUnivFuncsK2_',Test_ID,'_2ndExplore.mat'];
filename_CalGamma = ['ForTestCalGamma_',Test_ID,'_2ndExplore.mat'];

filename_PredictV = ['PredictingvData_',Test_ID,'_2ndExplore.mat'];

SigxyUseFinish_Flag = 1;
Tau_Pk2ResUseFinish_Flag = 1;


%%
ct = 1361;% s wave speed
%cl = 2680;%Plane Strain
cl = 2345;%Plane Stress, effective dilitational wave speed in order to use plain strain solution
rho = 1170;
mu = rho*ct*ct;
lambda = rho*cl*cl-2*mu;
nu = (cl^2-2*ct^2)/(2*(cl^2-ct^2));
E_young = mu*(3*lambda+2*mu)/(lambda+mu);


% estimating cR
%cR = 0.918*ct;
% actually calculating cR
func_cR_allPara = @(cR,ct,cl) (2.*(1./cR).^2-1./ct.^2).^2 ...
                                -4.*(1./cR).^2.*((1./cR).^2-1./cl.^2).^0.5.*((1./cR).^2-1./ct.^2).^0.5;
func_cR = @(cR) func_cR_allPara(cR,ct,cl);
x0 = 0.92*ct;
cR = fzero(func_cR,x0);
%cR./ct

%% Calculate A2 universal function, equation (1) in Svetlizky et al., 2017, Supp
% Freund 1990 around (5.3.10)

Vrup = [0:0.001:0.999]*cR;
alpha_d = sqrt(1-Vrup.^2./cl.^2);
alpha_s = sqrt(1-Vrup.^2./ct.^2);
FactorD = 4.*alpha_d.*alpha_s - (1+alpha_s.^2).^2;

A2_univfunc = (Vrup.^2.*alpha_s)./((1-nu).*ct.^2.*FactorD);

% figure(1)
% plot(Vrup./cR,A2_univfunc)
% xlabel('Vrup, cR')
% ylabel('A2_univfunc')

%% Load information from prepared data

%ForTestUnivFuncsK2 = load('ForTestUnivFuncsK2_case6_2ndExplore.mat');
ForTestUnivFuncsK2 = load(filename_FuncksK2);

CrackL = ForTestUnivFuncsK2.CrackL;

if SigxyUseFinish_Flag == 0
    % use final stress
    CrackDeltaSigxy = ForTestUnivFuncsK2.CrackDeltaSigxy;
elseif SigxyUseFinish_Flag == 1
    %%%% A MAJOR CHANGE in 2ndExplore, use finishing
    % use finishing stress
    CrackDeltaSigxy = ForTestUnivFuncsK2.CrackDeltaSigxy_Usefinishing;
else
    error('Wrong Flag')
end


delta_l = CrackL(2) - CrackL(1);
CrackL = CrackL + delta_l./2; % use middle point as the x eval location

%ForTestCalGamma = load('ForTestCalGamma_case6_2ndExplore.mat');
ForTestCalGamma = load(filename_CalGamma);
Gamma_StressSlip_RupLoc = ForTestCalGamma.Gamma_StressSlip_RupLoc;
Gamma_StressSlip_RupLoc_MinusAvgStress = ForTestCalGamma.Gamma_StressSlip_RupLoc_MinusAvgStress;
Gamma_StressSlip_RupLoc_ShorterCurve = ForTestCalGamma.Gamma_StressSlip_RupLoc_ShorterCurve;
Gamma_StressSlip_RupLoc_ShorterCurve_MinusAvgStress = ForTestCalGamma.Gamma_StressSlip_RupLoc_ShorterCurve_MinusAvgStress;
%Gamma_StressSlip_RupLoc = ForTestCalGamma.Gamma_StressSlip_RupLoc_ShorterCurve;
%Gamma_StressSlip_RupLoc = ForTestCalGamma.Gamma_StressSlip_RupLoc_ShorterCurve_MinusAvgStress;


XLoc = ForTestCalGamma.XLoc;

if Tau_Pk2ResUseFinish_Flag == 0
    % use finishing stress
    Tau_Pk2Res_RupLoc = ForTestCalGamma.Tau_Pk2Res_RupLoc;
elseif Tau_Pk2ResUseFinish_Flag == 1
    % use final stress
    Tau_Pk2Res_RupLoc = ForTestCalGamma.Tau_Pk2Res_UseFinishing_RupLoc;
else
    error('Wrong Flag')    
end

ADropWidth_SmoothWidth = 3;
ADropWidth = ForTestCalGamma.ADropWidth;
ADropWidth_forCalGamma = smooth(ADropWidth,ADropWidth_SmoothWidth);
ADropWidth_forCalGamma = ADropWidth_forCalGamma';
%ADropWidth_forCalGamma(:) = 0.02; % fancy version

RupSpeed = ForTestCalGamma.RupSpeed;
RupSpeed_SmoothforCalGamma = ForTestCalGamma.RupSpeed_smooth_forNote';


%% calculating K2S and G2S, equation (4) in Svetlizky et al., 2017, Supp



K2S = nan(size(CrackL));
for ind=1:length(CrackL)
   
    % Judging from the expression, I think when l = 0, K2S =0, because the
    % integral is from 0 to 0, even though 2/sqrt(pi*l) and the kernel
    % explode
    if ind==1
       K2S(ind) = 0;
       continue;
    end
     % there is singularity, for now, assume the integration is regularized
     % so I only do numerical integration to ind-1
    CurrentEvalL = CrackL(ind);
    CurrentCrackFromBegin2L_sArray = CrackL(1:ind-1);
    CurrentCrackFromBegin2L_SigArray = CrackDeltaSigxy(1:ind-1);
    F_CurrentEvalL = 1+0.3*(1-(CurrentCrackFromBegin2L_sArray./CurrentEvalL).^(5/4));    
    K2S_intKernel = ...
        CurrentCrackFromBegin2L_SigArray.*F_CurrentEvalL./sqrt(1-(CurrentCrackFromBegin2L_sArray./CurrentEvalL).^2);
    if length(CurrentCrackFromBegin2L_sArray) == 1
        Integral_CumtrapzTemp = K2S_intKernel.*delta_l;
    else
        Integral_CumtrapzTemp = cumtrapz(CurrentCrackFromBegin2L_sArray,K2S_intKernel);
    end    
    K2S_int = Integral_CumtrapzTemp(end);
    K2S(ind) = 2*K2S_int/sqrt(pi*CurrentEvalL);
end

figure(2)
plot(CrackL.*1e3,K2S,'-*')
xlabel('CrackL, mm')
ylabel('K2S')

xlim([0 200])

%% calculate G2s in equation (3) in Svetlizky et al., 2017, Supp

%vq = interp1(x,v,xq)
% K2S
% CrackL

% interpolate, only calculate location that has  for 
K2S_forXLoc = interp1(CrackL,K2S,XLoc);

% 10/24/2023: there are two possibilities, since we are doing plain strain,
% havem't figure out which one is correct. Depends on the equation (4) in
% Svetlizky et al., 2017, Supp is thinking about plain stress or plain
% strain.

% Need to check the use of wave speed & elastic modulus is consistent in
% this fracture mechanics analysis later
%G2S_forXLoc = 1./E_young.*K2S_forXLoc.^2;
G2S_forXLoc = (1-nu.^2)./E_young.*K2S_forXLoc.^2; % This one fit the data better, and if the equation (4) in
                                                   % Svetlizky et al.,
                                                   % 2017, Supp is thinking
                                                   % about plain strain,
                                                   % then the analysis here
                                                   % seems to be consistent
                                                   % throughout



%% Calculate fracture energy using Svet17 method.
%\Gamma in the "Measurements of \Gamma" section in Svetlizky et al., 2017, Supp

Intgral_in_K2GammaCal = sqrt(pi);
K2_inGammaCal = Tau_Pk2Res_RupLoc.*sqrt(ADropWidth_forCalGamma).*sqrt(2/pi).*Intgral_in_K2GammaCal;

%Calculate f2_cf for Gamma calculation
alpha_d_Gamma = sqrt(1-RupSpeed_SmoothforCalGamma.^2./cl.^2);
alpha_s_Gamma = sqrt(1-RupSpeed_SmoothforCalGamma.^2./ct.^2);
FactorD_Gamma = 4.*alpha_d_Gamma.*alpha_s_Gamma - (1+alpha_s_Gamma.^2).^2;

f2_cf_Gamma = (RupSpeed_SmoothforCalGamma.^2.*alpha_s_Gamma)./((1-nu).*ct.^2.*FactorD_Gamma);



%Gamma = 1./E_young.*f2_cf_Gamma.*K2_inGammaCal.^2;
Gamma = (1-nu.^2)./E_young.*f2_cf_Gamma.*K2_inGammaCal.^2;

% figure(3)
% plot(XLoc.*1e3,Gamma,'-*')
% xlabel('XLoc, mm')
% ylabel('Gamma')
% 
% xlim([0 200])

%% plot G2S and Gamma separately
clear p legend_str ind_p
figure(104)

ind_p = 0;

ind_p = ind_p + 1;
p(ind_p) = plot(XLoc.*1e3,G2S_forXLoc,'-d','linewidth',2,'Color',DClr.c(ind_p,:));
legend_str{ind_p} = 'G_2^S with Svet17 Method';

hold on

ind_p = ind_p + 1;
p(ind_p) = plot(XLoc.*1e3,Gamma,'-d','linewidth',2,'Color',DClr.c(ind_p,:));
legend_str{ind_p} = '\Gamma with Svet17 Method';

ind_p = ind_p + 1;
p(ind_p) = plot(XLoc.*1e3,Gamma_StressSlip_RupLoc,'-d','linewidth',2,'Color',DClr.c(ind_p,:));
legend_str{ind_p} = '\Gamma with Stress-slip, final slip local';

ind_p = ind_p + 1;
p(ind_p) = plot(XLoc.*1e3,Gamma_StressSlip_RupLoc_MinusAvgStress,'-d','linewidth',2,'Color',DClr.c(ind_p,:));
legend_str{ind_p} = '\Gamma with Stress-slip, final slip Avg';
% 
% ind_p = ind_p + 1;
% p(ind_p) = plot(XLoc.*1e3,Gamma_StressSlip_RupLoc_ShorterCurve,'-d','linewidth',2,'Color',DClr.c(ind_p,:));
% legend_str{ind_p} = '\Gamma with Stress-slip, finishing slip Local';
% 
ind_p = ind_p + 1;
p(ind_p) = plot(XLoc.*1e3,Gamma_StressSlip_RupLoc_ShorterCurve_MinusAvgStress,'-d','linewidth',2,'Color',DClr.c(ind_p,:));
legend_str{ind_p} = '\Gamma with Stress-slip, finishing slip Avg';


hold off

% for exploration
% plot(XLoc.*1e3,Gamma_StressSlip_RupLoc_MinusAvgStress,'-d','linewidth',2,'Color',DClr.c(3,:))
% hold off

xlabel('Along fault distance, mm')
ylabel('J/m^2')

legend(p,legend_str,'location','best')
xlim([0 200])

grid on

set(gca,'Fontsize',20,'Fontweight','bold')
set(gca, 'FontName', 'Helvetica')
set(gcf, 'Renderer', 'Painters');% make eps clear

% plot Cf/cR against G2S/Gamma
%cR = 0.92*ct;
% figure(5)
% plot(G2S_forXLoc./Gamma,RupSpeed./cR,'*')
% xlabel('G2S/\Gamma')
% ylabel('Cf/CR')
% 
% xlim([0 10])
% ylim([0 1])
% 
% set(gca,'xscale','linear')






%% about the k2(Cf) in equation (2) in Svetlizky et al., 2017, Supp
% Freund 1990 around (6.4.26), (6.4.32), (6.4.35) !!!!!! This is wrong
% !!!! Should use (6.4.41), (6.4.42) for k2(Cf); (6.4.32), (6.4.35) are for
% k1(Cf) , the only difference is need to change k_h and k_h_modified expression
% also see in ipad note

% Defined earlier
% cR = 0.92*ct;
% Vrup = [0:0.01:1]*0.92*ct;
h = 1./Vrup;
c = 1./cR;
a = 1./cl;
b = 1./ct;

a_plus = a./(1+a./h);
a_minus = a./(1-a./h);
b_plus = b./(1+b./h);
b_minus = b./(1-b./h);

fun_alpha = @(zeta,a,h) (a.^2-zeta.^2+a.^2.*zeta.^2./h.^2-2.*a.^2.*zeta./h).^0.5;
fun_beta = @(zeta,b,h) (b.^2-zeta.^2+b.^2.*zeta.^2./h.^2-2.*b.^2.*zeta./h).^0.5;

fun_V_minus = @(eta,a,b,h) (4.*eta^2.*fun_beta(eta,b,h).*abs(fun_alpha(eta,a,h)))./...
    (2.*eta.^2-b.^2-b.^2.*eta.^2./h.^2+2.*b.^2.*eta./h).^2;
fun_V_plus = @(eta,a,b,h) (4.*eta.^2.*fun_beta(eta,b,h).*abs(fun_alpha(eta,a,h)))./...
    (2.*eta.^2-b.^2-b.^2.*eta.^2./h.^2-2.*b.^2.*eta./h).^2;
% put an abs of fun_beta, so the output is real. Different from the given
% solution of (6.4.18), but can generate plot of Fig 6.10
fun_V_plus_modified = @(eta,a,b,h) (4.*eta.^2.*abs(fun_beta(eta,b,h)).*abs(fun_alpha(eta,a,h)))./...
    (2.*eta.^2-b.^2-b.^2.*eta.^2./h.^2-2.*b.^2.*eta./h).^2;

fun_SKernel_plus = @(zeta,eta,a,b,h) atan(fun_V_plus(eta,a,b,h))./(eta+zeta);
fun_SKernel_minus = @(zeta,eta,a,b,h) atan(fun_V_minus(eta,a,b,h))./(eta-zeta);

fun_SKernel_plus_modified = @(zeta,eta,a,b,h) atan(fun_V_plus_modified(eta,a,b,h))./(eta+zeta);

S_plus = nan(size(Vrup));
S_plus_modified = nan(size(Vrup));
for ind = 1:length(S_plus)
    Vrup_now = Vrup(ind);
    h_now = 1./Vrup_now;
    a_minus_now = a_minus(ind);
    b_minus_now = b_minus(ind);
    S_plus(ind) = exp(-1/pi* integral(@(eta) fun_SKernel_plus(h_now,eta,a,b,h_now),a_minus_now,b_minus_now));
    S_plus_modified(ind) = exp(-1/pi* integral(@(eta) fun_SKernel_plus_modified(h_now,eta,a,b,h_now),a_minus_now,b_minus_now));
end

% This is for k1 !!!!
k_h = (1-c./h)./(S_plus.*sqrt(1-a./h));
k_h_modified = (1-c./h)./(S_plus_modified.*sqrt(1-a./h));
k2_h = (1-c./h)./(S_plus.*sqrt(1-b./h));
k2_h_modified = (1-c./h)./(S_plus_modified.*sqrt(1-b./h));

k2_h_approx = (1-c./h)./(sqrt(1-b./h));

% figure(6)
% plot(Vrup./cR,S_plus)


% figure(7)
% plot(Vrup./cR,k_h)
% xlabel('Vrup, cR')
% ylabel('kh')
% title('Original form with fun_beta')
% 
% figure(8)
% plot(Vrup./cR,k_h_modified)
% xlabel('Vrup, cR')
% ylabel('kh')
% title('modified form with abs(fun_beta)')

% figure(9)
% plot(Vrup./cR,k2_h)
% xlabel('Vrup, cR')
% ylabel('kh')
% title('Original form with fun_beta')
%
figure(10)
plot(Vrup./cR,k2_h_modified)
hold on
plot(Vrup./cR,k2_h_approx)
xlabel('Vrup, cR')
ylabel('kh')
title('modified form with abs(\beta), vs approximate form')
legend('full','approx')


%% Predicting rupture speed with theoretical relation and test with numerical results
% combine A2 (f2) universal function with k(h) to get the theoretical relation between Vrup/cR versus G2S/Gamma
% Freund 1990 around (6.4.32) to (6.4.35). Svetlizky et al. 2017 SM
% equation (3)

% G2S/Gamma = 1/(A2(v)*k(v)), solving v


G2V = (A2_univfunc.*k2_h_modified.^2);
G2S_Gamma_ratio_theory = 1./G2V;

% figure(21)
% plot(Vrup./cR,G2V);
% xlabel('Vrup./cR')
% ylabel('GV')

% theory vs data, Gamma calculated by Svet17
clear p
figure(22)
p(1) = plot(G2S_Gamma_ratio_theory,Vrup./cR);
hold on
p(2) = plot(G2S_forXLoc./Gamma,RupSpeed_SmoothforCalGamma./cR,'*');
%%%% These are old plot for exploration
% Gamma_fancy = Gamma;
% Gamma_fancy(:) = 3;
%p(2) = plot(G2S_forXLoc./Gamma,RupSpeed./cR,'*');
%p(2) = plot(G2S_forXLoc./Gamma_fancy,RupSpeed./cR,'*');
%p(2) = plot(G2S_forXLoc./Gamma_StressSlip_RupLoc,RupSpeed./cR,'*');
%p(2) = plot(G2S_forXLoc./Gamma_StressSlip_RupLoc,RupSpeed_smooth./cR,'*');
%p(2) = plot(G2S_forXLoc./Gamma_fancy,RupSpeed_smooth./cR,'*');
%Gamma_case2 = load("Gamma_case2.mat");
%Gamma_case2_resmple = interp1(Gamma_case2.XLoc,Gamma_case2.Gamma,XLoc);
%p(2) = plot(G2S_forXLoc./Gamma_case2_resmple,RupSpeed_smooth./cR,'*');
%Multiple G2S/Gamma with a constant
%p(2) = plot(0.8*G2S_forXLoc./Gamma,RupSpeed_SmoothforCalGamma./cR,'*');

hold off


legend_str = {'theory (Svet17)','simulation 1'}
xlabel('G2S/\Gamma')
ylabel('Vrup./cR')

legend(p,legend_str,'location','best')
xlim([0 10])
ylim([0 1])

set(gca,'xscale','linear')


%% theory vs data, Gamma calculated by Stress slip
clear p
figure(23)
p(1) = plot(G2S_Gamma_ratio_theory,Vrup./cR);
hold on
p(2) = plot(G2S_forXLoc./Gamma_StressSlip_RupLoc,RupSpeed_SmoothforCalGamma./cR,'*');
% show okay fit at initial part
p(3) = plot(G2S_forXLoc./Gamma_StressSlip_RupLoc_MinusAvgStress,RupSpeed_SmoothforCalGamma./cR,'*');

%%%% These are old plot for exploration
% Gamma_fancy = Gamma;
% Gamma_fancy(:) = 3;
%p(2) = plot(G2S_forXLoc./Gamma,RupSpeed./cR,'*');
%p(2) = plot(G2S_forXLoc./Gamma_fancy,RupSpeed./cR,'*');
%p(2) = plot(G2S_forXLoc./Gamma_StressSlip_RupLoc,RupSpeed./cR,'*');
%p(2) = plot(G2S_forXLoc./Gamma_StressSlip_RupLoc,RupSpeed_smooth./cR,'*');
%p(2) = plot(G2S_forXLoc./Gamma_fancy,RupSpeed_smooth./cR,'*');
%Gamma_case2 = load("Gamma_case2.mat");
%Gamma_case2_resmple = interp1(Gamma_case2.XLoc,Gamma_case2.Gamma,XLoc);
%p(2) = plot(G2S_forXLoc./Gamma_case2_resmple,RupSpeed_smooth./cR,'*');
%p(2) = plot(G2S_forXLoc./Gamma_StressSlip_RupLoc_MinusAvgStress,RupSpeed_SmoothforCalGamma./cR,'*');
%p(3) = plot(G2S_forXLoc./Gamma_StressSlip_RupLoc_ShorterCurve_MinusAvgStress,RupSpeed_SmoothforCalGamma./cR,'*');
%Multiple G2S/Gamma with a constant
%p(2) = plot(0.7*G2S_forXLoc./Gamma_StressSlip_RupLoc,RupSpeed_SmoothforCalGamma./cR,'*');
%Multiple G2S/Gamma with a original speed
%p(2) = plot(G2S_forXLoc./Gamma_StressSlip_RupLoc,RupSpeed./cR,'*');
%p(3) = plot(G2S_forXLoc./Gamma_StressSlip_RupLoc_MinusAvgStress,RupSpeed_SmoothforCalGamma./cR,'*');

hold off


legend_str = {'theory (Stress-slip)','Stress slip Local','Minus Average Stress'}
xlabel('G2S/\Gamma')
ylabel('Vrup./cR')

legend(p,legend_str,'location','best')
xlim([0 10])
ylim([0 1])

set(gca,'xscale','linear')

%% save pedicting v plot data for the merge plot

save(filename_PredictV,'XLoc','A2_univfunc','k2_h_modified','Vrup','cR','G2S_forXLoc',...
    'Gamma_StressSlip_RupLoc_MinusAvgStress','Gamma','Gamma_StressSlip_RupLoc','RupSpeed_SmoothforCalGamma');
