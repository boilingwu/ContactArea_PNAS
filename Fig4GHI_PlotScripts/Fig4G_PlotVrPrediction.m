% Collapse curves of Fracture Mechanics predictions for models
clear all
clc



% Directory to load data,
%SimsData_Dir = '/Users/baoning/Dropbox/2DDynamicCycle/GRL_Draft/AprMatchSvet_ForPaper/May20_L400Simulation/W400mmExperiment/OpenEndAging_PStress/';
FracAnalysis_RootDir = './SecondAttempt/';
DClr = load(['ColorScheme_MatlabOrder.mat']);



%% plot one simulation, prepare data
% Simulation ID
ModelType_now = 'S';
%ModelType_now = 'F';

ModelNum_now = 5;


% set up names based on the chosen model
ModelID_now = ['P',ModelType_now,num2str(ModelNum_now,'%d')];
ModelFolderName_now = ['Model',ModelID_now];

filename_FuncksK2 = ['ForTestUnivFuncsK2_',ModelID_now,'_2ndExplore.mat'];
filename_CalGamma = ['ForTestCalGamma_',ModelID_now,'_2ndExplore.mat'];
filename_PredictV = ['PredictingvData_',ModelID_now,'_2ndExplore.mat'];

% load_data
PredictVrupData_now = load([FracAnalysis_RootDir,ModelFolderName_now ,'/',filename_PredictV]);

% loading theoretical curve from one experiment (all data have, just load from one)
Vrup = PredictVrupData_now.Vrup;
cR = PredictVrupData_now.cR;
A2_univfunc = PredictVrupData_now.A2_univfunc;
k2_h_modified = PredictVrupData_now.k2_h_modified;
G2V = (A2_univfunc.*k2_h_modified.^2);
G2S_Gamma_ratio_theory = 1./G2V;


% plot one simulation, plot
% theory vs data, Gamma calculated by Svet17
clear p legend_str p_ind
p_ind = 0;

figure(22)

p_ind = p_ind + 1;
p(p_ind) = plot(G2S_Gamma_ratio_theory,Vrup./cR,'-','linewidth',3,'Color','k'); % theory
legend_str{p_ind} = 'theory';

hold on

p_ind = p_ind + 1;
p(p_ind) = plot(PredictVrupData_now.G2S_forXLoc./PredictVrupData_now.Gamma,...
    PredictVrupData_now.RupSpeed_SmoothforCalGamma./cR,...
    'o','MarkerFaceColor',DClr.c(1,:),'MarkerEdgeColor',DClr.c(1,:)); %Svet Gamma
legend_str{p_ind} = ['Simulation',ModelID_now,'Svet Method'];

p_ind = p_ind + 1;
p(p_ind) = plot(PredictVrupData_now.G2S_forXLoc./PredictVrupData_now.Gamma_StressSlip_RupLoc,...
    PredictVrupData_now.RupSpeed_SmoothforCalGamma./cR,...
    'o','MarkerFaceColor',DClr.c(2,:),'MarkerEdgeColor',DClr.c(2,:)); %Stress-Slip Local
legend_str{p_ind} = ['Simulation',ModelID_now,'StressSlipLocal Method'];

p_ind = p_ind + 1;
p(p_ind) = plot(PredictVrupData_now.G2S_forXLoc./PredictVrupData_now.Gamma_StressSlip_RupLoc_MinusAvgStress,...
    PredictVrupData_now.RupSpeed_SmoothforCalGamma./cR,...
    'o','MarkerFaceColor',DClr.c(3,:),'MarkerEdgeColor',DClr.c(3,:)); %Stress Slip MinAvg
legend_str{p_ind} = ['Simulation',ModelID_now,'StressSlipMinAvg Method'];


% 
% p(4) = plot(PredictVrupData_case5.G2S_forXLoc./PredictVrupData_case5.Gamma,...
%     PredictVrupData_case5.RupSpeed_SmoothforCalGamma./cR,...
%     'o','MarkerFaceColor',DClr.c(3,:),'MarkerEdgeColor',DClr.c(3,:)); %case 5
% 
% p(5) = plot(PredictVrupData_case6.G2S_forXLoc./PredictVrupData_case6.Gamma,...
%     PredictVrupData_case6.RupSpeed_SmoothforCalGamma./cR,...
%     'o','MarkerFaceColor',DClr.c(4,:),'MarkerEdgeColor',DClr.c(4,:)); %case 6

hold off

xlabel('G2S/\Gamma')
ylabel('v/c_R')
title(['Model ',ModelID_now])

legend(p,legend_str,'location','best')
xlim([0 10])
ylim([0 1])

grid on
set(gca,'xscale','linear')

set(gca,'Fontsize',20,'Fontweight','bold')
set(gca, 'FontName', 'Helvetica')
set(gcf, 'Renderer', 'Painters');% make eps clear



%% Plot multiple models, v1


% some plot parameters
SvetDotEdgeAlpha = 1;
SvetDotFaceAlpha = 0.5;


clear p legend_str p_ind
p_ind = 0;
figure(101)

% plot theory
p_ind = p_ind + 1;
p(p_ind) = plot(G2S_Gamma_ratio_theory,Vrup./cR,'-','linewidth',5,'Color','k'); % theory
legend_str{p_ind} = 'theory';

hold on

%%%% Part 1: Plot PSX Model, Svet Method
ModelID_SArray_Svet = [1:3,7:10,14:29];
ModelType_now = 'S';
plot_count = 0;
for ModelNum_now = ModelID_SArray_Svet
    
    plot_count = plot_count + 1;

    % set up names based on the chosen model
    ModelID_now = ['P',ModelType_now,num2str(ModelNum_now,'%d')];
    ModelFolderName_now = ['Model',ModelID_now];
    
    filename_FuncksK2 = ['ForTestUnivFuncsK2_',ModelID_now,'_2ndExplore.mat'];
    filename_CalGamma = ['ForTestCalGamma_',ModelID_now,'_2ndExplore.mat'];
    filename_PredictV = ['PredictingvData_',ModelID_now,'_2ndExplore.mat'];
    
    % load_data
    PredictVrupData_now = load([FracAnalysis_RootDir,ModelFolderName_now ,'/',filename_PredictV]);   
    
    
    p_ind = p_ind + 1;
    % p(p_ind) = plot(PredictVrupData_now.G2S_forXLoc./PredictVrupData_now.Gamma,...
    %     PredictVrupData_now.RupSpeed_SmoothforCalGamma./cR,...
    %     'o','Color',DClr.c(mod(p_ind,7)+1,:),...
    %     'MarkerFaceColor',DClr.c(mod(p_ind,7)+1,:),'MarkerEdgeColor',DClr.c(mod(p_ind,7)+1,:),...
    %     'MarkerSize',1);
    p(p_ind) = scatter(PredictVrupData_now.G2S_forXLoc./PredictVrupData_now.Gamma,...
        PredictVrupData_now.RupSpeed_SmoothforCalGamma./cR,30,...
        'o',...
        'MarkerFaceColor',DClr.c(mod(p_ind,7)+1,:),'MarkerEdgeColor',DClr.c(mod(p_ind,7)+1,:),...
        'MarkerEdgeAlpha',SvetDotEdgeAlpha, 'MarkerFaceAlpha',SvetDotFaceAlpha);

    legend_str{p_ind} = ['Simulation',ModelID_now,'Svet Method'];    

end


%%%% Part 2: Plot PFX Model, Svet Method
ModelID_FArray_Svet = [5];
ModelType_now = 'F';
plot_count = 0;
for ModelNum_now = ModelID_FArray_Svet
    
    plot_count = plot_count + 1;

    % set up names based on the chosen model
    ModelID_now = ['P',ModelType_now,num2str(ModelNum_now,'%d')];
    ModelFolderName_now = ['Model',ModelID_now];
    
    filename_FuncksK2 = ['ForTestUnivFuncsK2_',ModelID_now,'_2ndExplore.mat'];
    filename_CalGamma = ['ForTestCalGamma_',ModelID_now,'_2ndExplore.mat'];
    filename_PredictV = ['PredictingvData_',ModelID_now,'_2ndExplore.mat'];
    
    % load_data
    PredictVrupData_now = load([FracAnalysis_RootDir,ModelFolderName_now ,'/',filename_PredictV]);   
    
    
    p_ind = p_ind + 1;
    % p(p_ind) = plot(PredictVrupData_now.G2S_forXLoc./PredictVrupData_now.Gamma,...
    %     PredictVrupData_now.RupSpeed_SmoothforCalGamma./cR,...
    %     'o','Color',DClr.c(mod(p_ind,7)+1,:),...
    %     'MarkerFaceColor',DClr.c(mod(p_ind,7)+1,:),'MarkerEdgeColor',DClr.c(mod(p_ind,7)+1,:),...
    %     'MarkerSize',1);
    p(p_ind) = scatter(PredictVrupData_now.G2S_forXLoc./PredictVrupData_now.Gamma,...
        PredictVrupData_now.RupSpeed_SmoothforCalGamma./cR,30,...
        'o',...
        'MarkerFaceColor',DClr.c(mod(p_ind,7)+1,:),'MarkerEdgeColor',DClr.c(mod(p_ind,7)+1,:),...
        'MarkerEdgeAlpha',SvetDotEdgeAlpha, 'MarkerFaceAlpha',SvetDotFaceAlpha);

    legend_str{p_ind} = ['Simulation',ModelID_now,'Svet Method'];    

end


%%%% Part 3: Plot PSX Model, Stress-Slip Method
ModelID_SArray_StressSlip = [];
ModelType_now = 'S';
plot_count = 0;
for ModelNum_now = ModelID_SArray_StressSlip
    
    plot_count = plot_count + 1;

    % set up names based on the chosen model
    ModelID_now = ['P',ModelType_now,num2str(ModelNum_now,'%d')];
    ModelFolderName_now = ['Model',ModelID_now];
    
    filename_FuncksK2 = ['ForTestUnivFuncsK2_',ModelID_now,'_2ndExplore.mat'];
    filename_CalGamma = ['ForTestCalGamma_',ModelID_now,'_2ndExplore.mat'];
    filename_PredictV = ['PredictingvData_',ModelID_now,'_2ndExplore.mat'];
    
    % load_data
    PredictVrupData_now = load([FracAnalysis_RootDir,ModelFolderName_now ,'/',filename_PredictV]);   
    
    
    p_ind = p_ind + 1;
    % p(p_ind) = plot(PredictVrupData_now.G2S_forXLoc./PredictVrupData_now.Gamma,...
    %     PredictVrupData_now.RupSpeed_SmoothforCalGamma./cR,...
    %     'o','Color',DClr.c(mod(p_ind,7)+1,:),...
    %     'MarkerFaceColor',DClr.c(mod(p_ind,7)+1,:),'MarkerEdgeColor',DClr.c(mod(p_ind,7)+1,:),...
    %     'MarkerSize',1);
    p(p_ind) = scatter(PredictVrupData_now.G2S_forXLoc./PredictVrupData_now.Gamma_StressSlip_RupLoc,...
        PredictVrupData_now.RupSpeed_SmoothforCalGamma./cR,30,...
        's',...
        'MarkerFaceColor',DClr.c(mod(p_ind,7)+1,:),'MarkerEdgeColor',DClr.c(mod(p_ind,7)+1,:),...
        'MarkerEdgeAlpha',SvetDotEdgeAlpha, 'MarkerFaceAlpha',SvetDotFaceAlpha);

    legend_str{p_ind} = ['Simulation',ModelID_now,'Stress-Slip Method'];    

end

%%%% Part 4: Plot PFX Model, Stress-Slip Method
ModelID_FArray_StressSlip = [1 4];
ModelType_now = 'F';
plot_count = 0;
for ModelNum_now = ModelID_FArray_StressSlip
    
    plot_count = plot_count + 1;

    % set up names based on the chosen model
    ModelID_now = ['P',ModelType_now,num2str(ModelNum_now,'%d')];
    ModelFolderName_now = ['Model',ModelID_now];
    
    filename_FuncksK2 = ['ForTestUnivFuncsK2_',ModelID_now,'_2ndExplore.mat'];
    filename_CalGamma = ['ForTestCalGamma_',ModelID_now,'_2ndExplore.mat'];
    filename_PredictV = ['PredictingvData_',ModelID_now,'_2ndExplore.mat'];
    
    % load_data
    PredictVrupData_now = load([FracAnalysis_RootDir,ModelFolderName_now ,'/',filename_PredictV]);   
    
    
    p_ind = p_ind + 1;
    % p(p_ind) = plot(PredictVrupData_now.G2S_forXLoc./PredictVrupData_now.Gamma,...
    %     PredictVrupData_now.RupSpeed_SmoothforCalGamma./cR,...
    %     'o','Color',DClr.c(mod(p_ind,7)+1,:),...
    %     'MarkerFaceColor',DClr.c(mod(p_ind,7)+1,:),'MarkerEdgeColor',DClr.c(mod(p_ind,7)+1,:),...
    %     'MarkerSize',1);
    p(p_ind) = scatter(PredictVrupData_now.G2S_forXLoc./PredictVrupData_now.Gamma_StressSlip_RupLoc,...
        PredictVrupData_now.RupSpeed_SmoothforCalGamma./cR,60,...
        's',...
        'MarkerFaceColor',DClr.c(mod(p_ind,7)+1,:),'MarkerEdgeColor',DClr.c(mod(p_ind,7)+1,:),...
        'MarkerEdgeAlpha',SvetDotEdgeAlpha, 'MarkerFaceAlpha',SvetDotFaceAlpha);

    legend_str{p_ind} = ['Simulation',ModelID_now,'Stress-Slip Method'];    

end


hold off

% plot set up
xlabel('G2S/\Gamma')
ylabel('v/c_R')
%title(['Model ',ModelID_now])

%legend(p,legend_str,'location','best')
xlim([0 10])
ylim([0 1])

grid on
set(gca,'xscale','linear')

set(gca,'Fontsize',20,'Fontweight','bold')
set(gca, 'FontName', 'Helvetica')
set(gcf, 'Renderer', 'Painters');% make eps clear


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Plot multiple models, v2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% some plot parameters
SvetDotEdgeAlpha = 0;
SvetDotFaceAlpha = 0.3;

StressSlipDotEdgeAlpha = 0;
StressSlipDotFaceAlpha = 0.5;

SvetDotColor = [0.5 0.5 0.5];
StressSlipDotColor = [0.4 0.4 0.4];

SvetDotSize = 30;
StressSlipDotSize = 60;

SvetDotType = 'o';
StressSlipDotType = '^';


clear p legend_str p_ind
figure(111)

% plot theory
p(1) = plot(G2S_Gamma_ratio_theory,Vrup./cR,'-','linewidth',5,'Color','k'); % theory
legend_str{1} = 'LEFM theory';

hold on

%%%% Part 1: Plot PSX Model, Svet Method
ModelID_SArray_Svet = [1:3,7:10,14:29];
ModelType_now = 'S';
plot_count = 0;
%p_ind = p_ind + 1;
for ModelNum_now = ModelID_SArray_Svet
    
    plot_count = plot_count + 1;

    % set up names based on the chosen model
    ModelID_now = ['P',ModelType_now,num2str(ModelNum_now,'%d')];
    ModelFolderName_now = ['Model',ModelID_now];
    
    filename_FuncksK2 = ['ForTestUnivFuncsK2_',ModelID_now,'_2ndExplore.mat'];
    filename_CalGamma = ['ForTestCalGamma_',ModelID_now,'_2ndExplore.mat'];
    filename_PredictV = ['PredictingvData_',ModelID_now,'_2ndExplore.mat'];
    
    % load_data
    PredictVrupData_now = load([FracAnalysis_RootDir,ModelFolderName_now ,'/',filename_PredictV]);   
    
    
    
    % p(p_ind) = plot(PredictVrupData_now.G2S_forXLoc./PredictVrupData_now.Gamma,...
    %     PredictVrupData_now.RupSpeed_SmoothforCalGamma./cR,...
    %     'o','Color',DClr.c(mod(p_ind,7)+1,:),...
    %     'MarkerFaceColor',DClr.c(mod(p_ind,7)+1,:),'MarkerEdgeColor',DClr.c(mod(p_ind,7)+1,:),...
    %     'MarkerSize',1);
    p_temp = scatter(PredictVrupData_now.G2S_forXLoc./PredictVrupData_now.Gamma,...
        PredictVrupData_now.RupSpeed_SmoothforCalGamma./cR,SvetDotSize,...
        SvetDotType,...
        'MarkerFaceColor',SvetDotColor,'MarkerEdgeColor',SvetDotColor,...
        'MarkerEdgeAlpha',SvetDotEdgeAlpha, 'MarkerFaceAlpha',SvetDotFaceAlpha);
    
    if plot_count == 1
        p(2) = p_temp;
        legend_str{2} = ['Simulations, Svet19 method'];    
    end

end


%%%% Part 2: Plot PFX Model, Svet Method
ModelID_FArray_Svet = [5];
ModelType_now = 'F';
plot_count = 0;
for ModelNum_now = ModelID_FArray_Svet
    
    plot_count = plot_count + 1;

    % set up names based on the chosen model
    ModelID_now = ['P',ModelType_now,num2str(ModelNum_now,'%d')];
    ModelFolderName_now = ['Model',ModelID_now];
    
    filename_FuncksK2 = ['ForTestUnivFuncsK2_',ModelID_now,'_2ndExplore.mat'];
    filename_CalGamma = ['ForTestCalGamma_',ModelID_now,'_2ndExplore.mat'];
    filename_PredictV = ['PredictingvData_',ModelID_now,'_2ndExplore.mat'];
    
    % load_data
    PredictVrupData_now = load([FracAnalysis_RootDir,ModelFolderName_now ,'/',filename_PredictV]);   
    
    
    % p(p_ind) = plot(PredictVrupData_now.G2S_forXLoc./PredictVrupData_now.Gamma,...
    %     PredictVrupData_now.RupSpeed_SmoothforCalGamma./cR,...
    %     'o','Color',DClr.c(mod(p_ind,7)+1,:),...
    %     'MarkerFaceColor',DClr.c(mod(p_ind,7)+1,:),'MarkerEdgeColor',DClr.c(mod(p_ind,7)+1,:),...
    %     'MarkerSize',1);
    p_temp = scatter(PredictVrupData_now.G2S_forXLoc./PredictVrupData_now.Gamma,...
        PredictVrupData_now.RupSpeed_SmoothforCalGamma./cR,SvetDotSize,...
        SvetDotType,...
        'MarkerFaceColor',SvetDotColor,'MarkerEdgeColor',SvetDotColor,...
        'MarkerEdgeAlpha',SvetDotEdgeAlpha, 'MarkerFaceAlpha',SvetDotFaceAlpha);
  

end


%%%% Part 3: Plot PSX Model, Stress-Slip Method
ModelID_SArray_StressSlip = [];
ModelType_now = 'S';
plot_count = 0;
for ModelNum_now = ModelID_SArray_StressSlip
    
    plot_count = plot_count + 1;

    % set up names based on the chosen model
    ModelID_now = ['P',ModelType_now,num2str(ModelNum_now,'%d')];
    ModelFolderName_now = ['Model',ModelID_now];
    
    filename_FuncksK2 = ['ForTestUnivFuncsK2_',ModelID_now,'_2ndExplore.mat'];
    filename_CalGamma = ['ForTestCalGamma_',ModelID_now,'_2ndExplore.mat'];
    filename_PredictV = ['PredictingvData_',ModelID_now,'_2ndExplore.mat'];
    
    % load_data
    PredictVrupData_now = load([FracAnalysis_RootDir,ModelFolderName_now ,'/',filename_PredictV]);   
    
    
    % p(p_ind) = plot(PredictVrupData_now.G2S_forXLoc./PredictVrupData_now.Gamma,...
    %     PredictVrupData_now.RupSpeed_SmoothforCalGamma./cR,...
    %     'o','Color',DClr.c(mod(p_ind,7)+1,:),...
    %     'MarkerFaceColor',DClr.c(mod(p_ind,7)+1,:),'MarkerEdgeColor',DClr.c(mod(p_ind,7)+1,:),...
    %     'MarkerSize',1);
    p_temp = scatter(PredictVrupData_now.G2S_forXLoc./PredictVrupData_now.Gamma_StressSlip_RupLoc,...
        PredictVrupData_now.RupSpeed_SmoothforCalGamma./cR,StressSlipDotSize,...
        StressSlipDotType,...
        'MarkerFaceColor',StressSlipDotColor,'MarkerEdgeColor',StressSlipDotColor,...
        'MarkerEdgeAlpha',StressSlipDotEdgeAlpha, 'MarkerFaceAlpha',StressSlipDotFaceAlpha);



end

%%%% Part 4: Plot PFX Model, Stress-Slip Method
ModelID_FArray_StressSlip = [1 4];
ModelType_now = 'F';
plot_count = 0;
for ModelNum_now = ModelID_FArray_StressSlip
    
    plot_count = plot_count + 1;

    % set up names based on the chosen model
    ModelID_now = ['P',ModelType_now,num2str(ModelNum_now,'%d')];
    ModelFolderName_now = ['Model',ModelID_now];
    
    filename_FuncksK2 = ['ForTestUnivFuncsK2_',ModelID_now,'_2ndExplore.mat'];
    filename_CalGamma = ['ForTestCalGamma_',ModelID_now,'_2ndExplore.mat'];
    filename_PredictV = ['PredictingvData_',ModelID_now,'_2ndExplore.mat'];
    
    % load_data
    PredictVrupData_now = load([FracAnalysis_RootDir,ModelFolderName_now ,'/',filename_PredictV]);   
    
    
    % p(p_ind) = plot(PredictVrupData_now.G2S_forXLoc./PredictVrupData_now.Gamma,...
    %     PredictVrupData_now.RupSpeed_SmoothforCalGamma./cR,...
    %     'o','Color',DClr.c(mod(p_ind,7)+1,:),...
    %     'MarkerFaceColor',DClr.c(mod(p_ind,7)+1,:),'MarkerEdgeColor',DClr.c(mod(p_ind,7)+1,:),...
    %     'MarkerSize',1);
    p_temp = scatter(PredictVrupData_now.G2S_forXLoc./PredictVrupData_now.Gamma_StressSlip_RupLoc,...
        PredictVrupData_now.RupSpeed_SmoothforCalGamma./cR,StressSlipDotSize,...
        StressSlipDotType,...
        'MarkerFaceColor',StressSlipDotColor,'MarkerEdgeColor',StressSlipDotColor,...
        'MarkerEdgeAlpha',StressSlipDotEdgeAlpha, 'MarkerFaceAlpha',StressSlipDotFaceAlpha);
    
    if plot_count == 1
        p(3) = p_temp;
        legend_str{3} = ['Simulations, Stress-Slip method'];    
    end     

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot the two models that are shown in details
%%% Fig3
ModelType_now = 'S';
ModelNum_now = 1;


% set up names based on the chosen model
ModelID_now = ['P',ModelType_now,num2str(ModelNum_now,'%d')];
ModelFolderName_now = ['Model',ModelID_now];
filename_FuncksK2 = ['ForTestUnivFuncsK2_',ModelID_now,'_2ndExplore.mat'];
filename_CalGamma = ['ForTestCalGamma_',ModelID_now,'_2ndExplore.mat'];
filename_PredictV = ['PredictingvData_',ModelID_now,'_2ndExplore.mat'];

% load_data
PredictVrupData_now = load([FracAnalysis_RootDir,ModelFolderName_now ,'/',filename_PredictV]);

p(4) = plot(PredictVrupData_now.G2S_forXLoc./PredictVrupData_now.Gamma,...
    PredictVrupData_now.RupSpeed_SmoothforCalGamma./cR,...
    'o','MarkerFaceColor',DClr.c(2,:),'MarkerEdgeColor',DClr.c(2,:)); %Svet Gamma
legend_str{4} = ['Simulation A, Svet19 method'];

%%%%%%%%%%%%%%%%%%%%%%%%
%%% Fig4

ModelType_now = 'S';
ModelNum_now = 2;

% set up names based on the chosen model
ModelID_now = ['P',ModelType_now,num2str(ModelNum_now,'%d')];
ModelFolderName_now = ['Model',ModelID_now];
filename_FuncksK2 = ['ForTestUnivFuncsK2_',ModelID_now,'_2ndExplore.mat'];
filename_CalGamma = ['ForTestCalGamma_',ModelID_now,'_2ndExplore.mat'];
filename_PredictV = ['PredictingvData_',ModelID_now,'_2ndExplore.mat'];

% load_data
PredictVrupData_now = load([FracAnalysis_RootDir,ModelFolderName_now ,'/',filename_PredictV]);

p(5) = plot(PredictVrupData_now.G2S_forXLoc./PredictVrupData_now.Gamma,...
    PredictVrupData_now.RupSpeed_SmoothforCalGamma./cR,...
    'o','MarkerFaceColor',DClr.c(3,:),'MarkerEdgeColor',DClr.c(3,:)); %Svet Gamma
legend_str{5} = ['Simulation B, Svet19 method'];

hold off

% plot set up
xlabel('G2S/\Gamma')
ylabel('Vr/cR')
%title(['Model ',ModelID_now])

legend(p,legend_str,'location','best')
xlim([0 8])
ylim([0 1])

grid on
box on
ax = gca;
ax.LineWidth = 2;

set(gca,'xscale','linear')

set(gca,'Fontsize',20,'Fontweight','bold')
set(gca, 'FontName', 'Helvetica')
set(gcf, 'Renderer', 'Painters');% make eps clear
