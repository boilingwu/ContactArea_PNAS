% Plot MaxVr against Ru*Rb

clear all 
clc

% Directory to load data,
%SimsData_Dir = '/Users/baoning/Dropbox/2DDynamicCycle/GRL_Draft/AprMatchSvet_ForPaper/May20_L400Simulation/W400mmExperiment/OpenEndAging_PStress/';
FracAnalysis_RootDir = './SecondAttempt/';
DClr = load(['ColorScheme_MatlabOrder.mat']);

%% Load Ru and Rb number

% PS
PSX_RuAndRb = ...
    [9.37	0.25
    6.25	0.25
    10.93	0.25
    13.12	0.25
    5.46	0.62
    10.93	0.33
    6.56	0.25
    5.46	0.25
    8.2	0.25
    10.93	0.25
    5.46	0.6
    4.37	0.6
    3.28	0.6
    5.46	0.25
    7.29	0.25
    10.93	0.25
    6.56	0.38
    6.56	0.43
    6.56	0.33
    6.56	0.3
    8.2	0.3
    5.46	0.27
    6.56	0.27
    8.2	0.27
    10.93	0.27
    6.56	0.23
    8.2	0.23
    5.83	0.23
    8.74	0.23];

% PFX
PFX_RuAndRb = ...
    [6.31	0.62
    8.2	0.5
    7.71	0.5
    6.56	0.6
    8.2	0.5];

Merge_RuAndRb = cat(1,PSX_RuAndRb,PFX_RuAndRb);
Merge_RuDotRb = nan(length(Merge_RuAndRb(:,1)),1);
Merge_RuDotRb(:) = Merge_RuAndRb(:,1).*Merge_RuAndRb(:,2);


%% load max Vr

Merge_MaxVr = nan(size(Merge_RuDotRb));
Model_count = 0;

%PSX model
ModelID_SArray = [1:29];
ModelType_now = 'S';

for ModelNum_now = ModelID_SArray
    
    Model_count = Model_count + 1;

    % set up names based on the chosen model
    ModelID_now = ['P',ModelType_now,num2str(ModelNum_now,'%d')];
    ModelFolderName_now = ['Model',ModelID_now];
    
    filename_FuncksK2 = ['ForTestUnivFuncsK2_',ModelID_now,'_2ndExplore.mat'];
    filename_CalGamma = ['ForTestCalGamma_',ModelID_now,'_2ndExplore.mat'];
    filename_PredictV = ['PredictingvData_',ModelID_now,'_2ndExplore.mat'];
    
    % load_data
    PredictVrupData_now = load([FracAnalysis_RootDir,ModelFolderName_now ,'/',filename_PredictV]);   
    
    
    Merge_MaxVr(Model_count) = max(PredictVrupData_now.RupSpeed_SmoothforCalGamma);

end


%PFX model
ModelID_FArray = [1:5];
ModelType_now = 'F';

for ModelNum_now = ModelID_FArray
    
    Model_count = Model_count + 1;

    % set up names based on the chosen model
    ModelID_now = ['P',ModelType_now,num2str(ModelNum_now,'%d')];
    ModelFolderName_now = ['Model',ModelID_now];
    
    filename_FuncksK2 = ['ForTestUnivFuncsK2_',ModelID_now,'_2ndExplore.mat'];
    filename_CalGamma = ['ForTestCalGamma_',ModelID_now,'_2ndExplore.mat'];
    filename_PredictV = ['PredictingvData_',ModelID_now,'_2ndExplore.mat'];
    
    % load_data
    PredictVrupData_now = load([FracAnalysis_RootDir,ModelFolderName_now ,'/',filename_PredictV]);   
    
    
    Merge_MaxVr(Model_count) = max(PredictVrupData_now.RupSpeed_SmoothforCalGamma);

end

%% plot MaxVr against Ru*Rb

% loading theoretical curve from one experiment (all data have, just load from one)
Vrup = PredictVrupData_now.Vrup;
cR = PredictVrupData_now.cR;
A2_univfunc = PredictVrupData_now.A2_univfunc;
k2_h_modified = PredictVrupData_now.k2_h_modified;
G2V = (A2_univfunc.*k2_h_modified.^2);
G2S_Gamma_ratio_theory = 1./G2V;


%%%%% remove the simulations that aren't used in Fracture Mechanics Analysis
%RemoveModelInd = [4 5 6 11 12 13 31 32];
%RemoveModelInd = [4 5 6 11 12 13 31 32];
RemoveModelInd = [1 2 4 5 6 11 12 13 31 32];
Merge_RuDotRb_RemoveUnused = Merge_RuDotRb;
Merge_MaxVr_RemoveUnused = Merge_MaxVr;
Merge_RuDotRb_RemoveUnused(RemoveModelInd) = nan;
Merge_MaxVr_RemoveUnused(RemoveModelInd) = nan;


clear p legend_str
figure(1)

% plot parameters
dot_size = 200;
dotEdge_color = [0.1 0.1 0.1];
%dotFace_color = DClr.c(1,:);
dotFace_color = [0.5 0.5 0.5];
%dot_color = DClr.c(1,:);

% scatter(Merge_RuDotRb,Merge_MaxVr./cR,60,'o',...
%         'MarkerFaceColor',dot_color,'MarkerEdgeColor',dot_color,...
%         'MarkerEdgeAlpha',1, 'MarkerFaceAlpha',1);
p(1) = scatter(Merge_RuDotRb_RemoveUnused,Merge_MaxVr_RemoveUnused./cR,dot_size,'o',...
        'MarkerFaceColor',dotFace_color,'MarkerEdgeColor',dotEdge_color,...
        'MarkerEdgeAlpha',1, 'MarkerFaceAlpha',0.8);
legend_str{1} = 'All simulations';

hold on

p(2) = scatter(Merge_RuDotRb(1),Merge_MaxVr(1)./cR,dot_size*1.5,'s',...
        'MarkerFaceColor',DClr.c(2,:),'MarkerEdgeColor',dotEdge_color,...
        'MarkerEdgeAlpha',1, 'MarkerFaceAlpha',1);
legend_str{2} = 'Simulation A';

p(3) = scatter(Merge_RuDotRb(2),Merge_MaxVr(2)./cR,dot_size*1.5,'d',...
        'MarkerFaceColor',DClr.c(3,:),'MarkerEdgeColor',dotEdge_color,...
        'MarkerEdgeAlpha',1, 'MarkerFaceAlpha',1);
legend_str{3} = 'Simulation B';

hold off

% plot set up
xlabel('Ru*Rb')
ylabel('Max Vr/cR')
%title(['Model ',ModelID_now])

%legend(p,legend_str,'location','best')
%xlim([0 10])
ylim([0 1])

legend(p,legend_str,'location','best')

grid on
box on
ax = gca;
ax.LineWidth = 2;


set(gca,'xscale','linear')

set(gca,'Fontsize',20,'Fontweight','bold')
set(gca, 'FontName', 'Helvetica')
set(gcf, 'Renderer', 'Painters');% make eps clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot MaxVr against [W/(mu*L/((b-a)*sig)] * (b-a)/b = [W/(mu*L/(b*sig)] * [(b-a)/b]^2
% the Current Ru*Rb number is calculated use the following equation in
% the simulation code, 
    %%%%%%%%%%%%%%%%%%%%
    % hstar = 2/pi*(mu/(1-nu)*Fault(ind).Dc)/((Fault(ind).b-Fault(ind).a)*Fault(ind).sig0)
    % 
    % L_b = 2/pi*(mu/(1-nu)*Fault(ind).Dc)/((Fault(ind).b)*Fault(ind).sig0)
    % 
    % hstar_min = 2/pi*(mu/(1-nu)*Fault(ind).Dc)/((Fault(ind).b-Fault(ind).a)*Fault(ind).sig0*max(ElementSigShape))
    % 
    % L_b_min = 2/pi*(mu/(1-nu)*Fault(ind).Dc)/((Fault(ind).b)*Fault(ind).sig0*max(ElementSigShape))
    % 
    % hstar_max = 2/pi*(mu/(1-nu)*Fault(ind).Dc)/((Fault(ind).b-Fault(ind).a)*Fault(ind).sig0*min(ElementSigShape))
    % 
    % L_b_max = 2/pi*(mu/(1-nu)*Fault(ind).Dc)/((Fault(ind).b)*Fault(ind).sig0*min(ElementSigShape))
    % 
    % Fault(ind).hstar=hstar;
    % Fault(ind).L_b=L_b;
    % 
    % Ru = Fault(1).Length/Fault(1).hstar;
    % Rb = (Fault(1).b - Fault(1).a)/Fault(1).b;
    % RuDotRb = Ru*Rb
    %%%%%%%%%%%%%%%%%%%%%%%
% Therefore, we need to multiple a coefficient to convert Ru*Rb to
% [W/(mu*L/((b-a)*sig)] * (b-a)/b
% Note that Fault(1).Length is the full simulation length in full space (400mm),
% for half-space, the length is (200mm)
%%%%%
%%% the coefficent should be
% [W/(mu*L/((b-a)*sig)] * (b-a)/b = Ru*Rb * pi/2*(1-nu)/2
%
%
%%% the simulation parameters



%%% material setup
ct = 1361;
%cl = 2680;%Plane Strain
cl = 2345;%Plane Stress
rho = 1170;
mu = rho*ct*ct;
lambda = rho*cl*cl-2*mu;
nu = (cl^2-2*ct^2)/(2*(cl^2-ct^2));


%%% scale

Merge_NewRatio = Merge_RuDotRb * (pi/2*(1-nu)/2);


%%%%% remove the simulations that aren't used in Fracture Mechanics Analysis
%RemoveModelInd = [4 5 6 11 12 13 31 32];
%RemoveModelInd = [4 5 6 11 12 13 31 32];
RemoveModelInd = [1 2 4 5 6 11 12 13 31 32];
Merge_NewRatio_RemoveUnused = Merge_NewRatio;
Merge_MaxVr_RemoveUnused = Merge_MaxVr;
Merge_NewRatio_RemoveUnused(RemoveModelInd) = nan;
Merge_MaxVr_RemoveUnused(RemoveModelInd) = nan;


clear p legend_str
figure(2)

% plot parameters
dot_size = 200;
dotEdge_color = [0.1 0.1 0.1];
%dotFace_color = DClr.c(1,:);
dotFace_color = [0.5 0.5 0.5];
%dot_color = DClr.c(1,:);

% scatter(Merge_RuDotRb,Merge_MaxVr./cR,60,'o',...
%         'MarkerFaceColor',dot_color,'MarkerEdgeColor',dot_color,...
%         'MarkerEdgeAlpha',1, 'MarkerFaceAlpha',1);
p(1) = scatter(Merge_NewRatio_RemoveUnused,Merge_MaxVr_RemoveUnused./cR,dot_size,'o',...
        'MarkerFaceColor',dotFace_color,'MarkerEdgeColor',dotEdge_color,...
        'MarkerEdgeAlpha',1, 'MarkerFaceAlpha',0.8);
legend_str{1} = 'All simulations';

hold on

p(2) = scatter(Merge_NewRatio(1),Merge_MaxVr(1)./cR,dot_size*1.5,'s',...
        'MarkerFaceColor',DClr.c(2,:),'MarkerEdgeColor',dotEdge_color,...
        'MarkerEdgeAlpha',1, 'MarkerFaceAlpha',1);
legend_str{2} = 'Simulation A';

p(3) = scatter(Merge_NewRatio(2),Merge_MaxVr(2)./cR,dot_size*1.5,'d',...
        'MarkerFaceColor',DClr.c(3,:),'MarkerEdgeColor',dotEdge_color,...
        'MarkerEdgeAlpha',1, 'MarkerFaceAlpha',1);
legend_str{3} = 'Simulation B';

hold off

% plot set up
xlabel('Ru*Rb*pi/2*(1-nu)/2')
ylabel('Max Vr/cR')
%title(['Model ',ModelID_now])

%legend(p,legend_str,'location','best')
%xlim([0 10])
ylim([0 1])

legend(p,legend_str,'location','best')

grid on
box on
ax = gca;
ax.LineWidth = 2;


set(gca,'xscale','linear')

set(gca,'Fontsize',20,'Fontweight','bold')
set(gca, 'FontName', 'Helvetica')
set(gcf, 'Renderer', 'Painters');% make eps clear




