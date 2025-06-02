% plot relative Ar drops vs tau drop

clear all
close all
clc

%data_folder = '/Users/baoning/Dropbox/2DDynamicCycle/GRL_Draft/AprMatchSvet_ForPaper/May20_L400Simulation/W400mmExperiment/OpenEndAging_PStress/';
DClr = load(['ColorScheme_MatlabOrder.mat']);


%% set up data

SvetFast_ArDrop = 0.4027;
SvetFast_TauDrop = 0.2326;

SvetSlow_ArDrop = 0.3502;
SvetSlow_TauDrop = 0.1247;

ModelDrops_mat = ...
    [0.3186	0.1683
0.3346	0.1662
0.3085	0.1674
0.2963	0.165
0.3832	0.349
0.234	0.161
0.1767	0.0837
0.1794	0.0822
0.1744	0.0837
0.1643	0.0826
0.1681	0.1517
0.1681	0.1506
0.1627	0.1433
0.2359	0.1096
0.2318	0.1124
0.2154	0.1114
0.2453	0.1607
0.2157	0.1586
0.267	0.1636
0.2907	0.165
0.2815	0.1654
0.3134	0.1629
0.3143	0.1657
0.3049	0.167
0.2849	0.1653
0.3586	0.1658
0.3481	0.1687
0.4386	0.2101
0.4376	0.2214
0.364	0.3552
0.3285	0.2988
0.3357	0.2978
0.1614	0.1524
0.3505	0.2989];

%%%%% remove the simulations that aren't used in Fracture Mechanics Analysis
%RemoveModelInd = [4 5 6 11 12 13 31 32];
%RemoveModelInd = [4 5 6 11 12 13 31 32];
RemoveModelInd = [1 2 4 5 6 11 12 13 31 32];
%RemoveModelInd = [1 2];
ModelDrops_mat_RemoveUnused = ModelDrops_mat;
ModelDrops_mat_RemoveUnused(RemoveModelInd,:) = nan;



%% make plot

clear p legend_str p_ind

p_ind = 0;
figure(1)

ModelEdgeClr = [0.7 0.7 0.7];
ModelFaceClr = ModelEdgeClr;

SvetFastEdgeClr = DClr.c(1,:);
SvetFastFaceClr = SvetFastEdgeClr;

SvetSlowEdgeClr = DClr.c(4,:);
SvetSlowFaceClr = SvetSlowEdgeClr;

ModelSlowEdgeClr = DClr.c(3,:);
ModelSlowFaceClr = ModelSlowEdgeClr;

ModelFastEdgeClr = DClr.c(2,:);
ModelFastFaceClr = ModelFastEdgeClr;

p_ind = p_ind + 1;
p(p_ind) = scatter(ModelDrops_mat_RemoveUnused(:,1),ModelDrops_mat_RemoveUnused(:,2),100,'o','MarkerEdgeColor',ModelEdgeClr,...
              'MarkerFaceColor',ModelFaceClr,...
              'LineWidth',1);
legend_str{p_ind} = 'All simulations';

hold on

p_ind = p_ind + 1;
p(p_ind) = scatter(ModelDrops_mat(1,1),ModelDrops_mat(1,2),300,'s','MarkerEdgeColor',ModelFastEdgeClr,...
              'MarkerFaceColor',ModelFastFaceClr,...
              'LineWidth',1);
legend_str{p_ind} = 'Simulation A';

p_ind = p_ind + 1;
p(p_ind) = scatter(ModelDrops_mat(2,1),ModelDrops_mat(2,2),300,'d','MarkerEdgeColor',ModelSlowEdgeClr,...
              'MarkerFaceColor',ModelSlowFaceClr,...
              'LineWidth',1);
legend_str{p_ind} = 'Simulation B';

p_ind = p_ind + 1;
p(p_ind) = scatter(SvetFast_ArDrop,SvetFast_TauDrop,300,'^','MarkerEdgeColor',SvetFastEdgeClr,...
              'MarkerFaceColor',SvetFastFaceClr,...
              'LineWidth',1);
legend_str{p_ind} = 'Lab Event 1';

p_ind = p_ind + 1;
p(p_ind) = scatter(SvetSlow_ArDrop,SvetSlow_TauDrop,300,'v','MarkerEdgeColor',SvetSlowEdgeClr,...
              'MarkerFaceColor',SvetSlowFaceClr,...
              'LineWidth',1);
legend_str{p_ind} = 'Lab Event 2';


p_ind = p_ind + 1;
p(p_ind) = plot([0 0.5],[0 0.5],'-','Color','k','LineWidth',2);
legend_str{p_ind} = '1:1 line';

hold off

xlabel('\DeltaAr/Ar');
ylabel('\Delta\tau/\tau');

xlim([0 0.5])
ylim([0 0.5])

legend(p,legend_str,'location','best')




%axis equal

grid on
box on
ax = gca;
ax.LineWidth = 2;

set(gca,'Fontsize',20,'Fontweight','bold')
set(gca, 'FontName', 'Helvetica')
set(gcf, 'Renderer', 'Painters');% make eps clear