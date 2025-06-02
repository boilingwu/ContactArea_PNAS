
% load simulation data

c0_workspace = load('vstep_workspace_c0.00_AgingLaw_Review2.mat');
D_series_aging = c0_workspace.D_series;
f_series_aging = c0_workspace.f_series;
A_series_aging = c0_workspace.A_series;
LoadPointD_results_aging = c0_workspace.LoadPointD_results;

c1_workspace = load('vstep_workspace_c0.00_SlipLaw_Review2.mat');
D_series_slip = c1_workspace.D_series;
f_series_slip = c1_workspace.f_series;
A_series_slip = c1_workspace.A_series;
LoadPointD_results_slip = c1_workspace.LoadPointD_results;
clearvars -except ...
    D_series_aging f_series_aging A_series_aging LoadPointD_results_aging ...
    D_series_slip f_series_slip A_series_slip LoadPointD_results_slip


% %% plot model against DK94 fig8 data
% load ColorScheme_MatlabOrder.mat
% ImportDK94_friction;
% ImportDK94_PA
% 
% data_linewidth = 6;
% aging_linewidth = 4;
% agingNagata_linewidth = 3;
% 
% data_color = [0.7 0.7 0.7];
% 
% 
% D_DK94Friction = DK94Fig8DataDigitizedFriction.D_DK94Friction;
% f_DK94Friction = DK94Fig8DataDigitizedFriction.f_DK94Friction;
% D_DK94PA = DK94Fig8DataDigitizedPercentArea.D_DK94PA;
% PA_DK94PA = DK94Fig8DataDigitizedPercentArea.PA_DK94PA;
% 
% figure(3)
% 
% yyaxis left
% plot(D_DK94Friction,f_DK94Friction,...
%     'linewidth',data_linewidth,'Linestyle','-','Color',data_color)
% hold on
% plot(D_series_aging.*1e6 + 400,f_series_aging,'linewidth',aging_linewidth,'Linestyle','-','Color',c(1,:))
% plot(D_series_slip.*1e6 + 400,f_series_slip,'linewidth',agingNagata_linewidth,'Linestyle','-.','Color',c(3,:))
% hold off
% 
% set(gca,'Ycolor','k');
% xlabel('Block Displacement, \mu m')
% ylabel('COF')
% ylim([0.5 0.9])
% 
% yyaxis right
% plot(D_DK94PA,PA_DK94PA,...
%     'linewidth',3,'Linestyle','-','Color',data_color)
% hold on
% plot(D_series_aging.*1e6 + 400,A_series_aging*100,'linewidth',aging_linewidth,'Linestyle','-','Color',c(2,:))
% plot(D_series_slip.*1e6 + 400,A_series_slip*100,'linewidth',agingNagata_linewidth,'Linestyle','-.','Color',c(6,:))
% hold off
% 
% set(gca,'Ycolor','k');
% %xlabel('Displacement, \mu m')
% ylabel('area percentage')
% ylim([0.3 0.8])
% 
% xlim([400 800])
% set(gca,'Fontsize',20,'Fontweight','bold')
% set(gca, 'FontName', 'Helvetica')
% set(gcf, 'Renderer', 'Painters');% make eps clear
% 
% 
% %% plot model against DK94 fig8 data, using loadpoint velocity
% clear p legend_str
% figure(4)
% 
% yyaxis left
% p(1) = plot(D_DK94Friction,f_DK94Friction,...
%     'linewidth',data_linewidth,'Linestyle','-','Color',data_color);
% hold on
% p(2) = plot(LoadPointD_results_aging.*1e6 + 400,f_series_aging,...
%     'linewidth',aging_linewidth,'Linestyle','-','Color',c(1,:));
% p(3) = plot(LoadPointD_results_slip.*1e6 + 400,f_series_slip,...
%     'linewidth',agingNagata_linewidth,'Linestyle','--','Color',c(2,:));
% hold off
% 
% ytickLoc = [0.65:0.05:0.9];
% ytickLabel = num2cell(ytickLoc);
% ytickLabel = string(ytickLabel);
% ytickLabel(1:2:end) = nan;
% yticks(ytickLoc)
% yticklabels(ytickLabel)
% 
% set(gca,'Ycolor','k');
% xlabel('Loadpoint Displacement, \mu m')
% ylabel('Friction')
% ylim([0.4 0.9])
% 
% yyaxis right
% plot(D_DK94PA,PA_DK94PA,...
%     'linewidth',data_linewidth,'Linestyle','-','Color',data_color)
% hold on
% plot(LoadPointD_results_aging.*1e6 + 400,A_series_aging*100,'linewidth',aging_linewidth,'Linestyle','-','Color',c(1,:))
% plot(LoadPointD_results_slip.*1e6 + 400,A_series_slip*100,'linewidth',agingNagata_linewidth,'Linestyle','--','Color',c(2,:))
% hold off
% 
% ytickLoc = [0.3:0.05:0.5];
% ytickLabel = num2cell(ytickLoc);
% ytickLabel = string(ytickLabel);
% ytickLabel(2:2:end) = nan;
% yticks(ytickLoc)
% yticklabels(ytickLabel)
% 
% set(gca,'Ycolor','k');
% %xlabel('Displacement, \mu m')
% ylabel('area percentage')
% ylim([0.3 0.8])
% grid on
% 
% % legend_str = {'lab data','Aging model','Nagata-aging model'};
% % legend(p,legend_str,'location','best')
% 
% xlim([400 800])
% set(gca,'Fontsize',22,'Fontweight','bold')
% set(gca, 'FontName', 'Helvetica')
% set(gcf, 'Renderer', 'Painters');% make eps clear

%% plot model against DK94 fig8 data, using loadpoint velocity, just aging law
clear p legend_str

load ColorScheme_MatlabOrder.mat
ImportDK94_friction;
ImportDK94_PA

data_linewidth = 6;
aging_linewidth = 4;
slip_linewidth = 3;

data_color = [0.7 0.7 0.7];


D_DK94Friction = DK94Fig8DataDigitizedFriction.D_DK94Friction;
f_DK94Friction = DK94Fig8DataDigitizedFriction.f_DK94Friction;
D_DK94PA = DK94Fig8DataDigitizedPercentArea.D_DK94PA;
PA_DK94PA = DK94Fig8DataDigitizedPercentArea.PA_DK94PA;

figure(5)

yyaxis left
p(1) = plot(D_DK94Friction,f_DK94Friction,...
    'linewidth',data_linewidth,'Linestyle','-','Color',data_color);
hold on
p(2) = plot(LoadPointD_results_aging.*1e6 + 400,f_series_aging,...
    'linewidth',aging_linewidth,'Linestyle','-','Color',c(1,:));
hold off

ytickLoc = [0.65:0.05:0.9];
ytickLabel = num2cell(ytickLoc);
ytickLabel = string(ytickLabel);
ytickLabel(1:2:end) = nan;
yticks(ytickLoc)
yticklabels(ytickLabel)

set(gca,'Ycolor','k');
xlabel('Loadpoint Displacement, \mu m')
ylabel('Friction')
ylim([0.4 0.9])

yyaxis right
plot(D_DK94PA,PA_DK94PA,...
    'linewidth',data_linewidth,'Linestyle','-','Color',data_color)
hold on
plot(LoadPointD_results_aging.*1e6 + 400,A_series_aging*100,'linewidth',aging_linewidth,'Linestyle','-','Color',c(1,:))
%plot(LoadPointD_results_slip.*1e6 + 400,A_series_slip*100,'linewidth',slip_linewidth,'Linestyle','--','Color',c(2,:))
hold off

ytickLoc = [0.3:0.05:0.5];
ytickLabel = num2cell(ytickLoc);
ytickLabel = string(ytickLabel);
ytickLabel(2:2:end) = nan;
yticks(ytickLoc)
yticklabels(ytickLabel)

set(gca,'Ycolor','k');
%xlabel('Displacement, \mu m')
ylabel('area percentage')
ylim([0.3 0.8])
grid on

legend_str = {'lab data','Aging model','Nagata-aging model'};
legend(p,legend_str,'location','best')

xlim([400 800])
set(gca,'Fontsize',22,'Fontweight','bold')
set(gca, 'FontName', 'Helvetica')
set(gcf, 'Renderer', 'Painters');% make eps clear

%% plot model against DK94 fig8 data, using loadpoint velocity, just slip law
clear p legend_str

load ColorScheme_MatlabOrder.mat
ImportDK94_friction;
ImportDK94_PA

data_linewidth = 6;
aging_linewidth = 4;
slip_linewidth = 4;

data_color = [0.7 0.7 0.7];


D_DK94Friction = DK94Fig8DataDigitizedFriction.D_DK94Friction;
f_DK94Friction = DK94Fig8DataDigitizedFriction.f_DK94Friction;
D_DK94PA = DK94Fig8DataDigitizedPercentArea.D_DK94PA;
PA_DK94PA = DK94Fig8DataDigitizedPercentArea.PA_DK94PA;

figure(6)

yyaxis left
p(1) = plot(D_DK94Friction,f_DK94Friction,...
    'linewidth',data_linewidth,'Linestyle','-','Color',data_color);
hold on
p(2) = plot(LoadPointD_results_slip.*1e6 + 400,f_series_slip,...
    'linewidth',slip_linewidth,'Linestyle','-','Color',c(2,:));
hold off

ytickLoc = [0.65:0.05:0.9];
ytickLabel = num2cell(ytickLoc);
ytickLabel = string(ytickLabel);
ytickLabel(1:2:end) = nan;
yticks(ytickLoc)
yticklabels(ytickLabel)

set(gca,'Ycolor','k');
xlabel('Loadpoint Displacement, \mu m')
ylabel('Friction')
ylim([0.4 0.9])

yyaxis right
plot(D_DK94PA,PA_DK94PA,...
    'linewidth',data_linewidth,'Linestyle','-','Color',data_color)
hold on
%plot(LoadPointD_results_aging.*1e6 + 400,A_series_aging*100,'linewidth',aging_linewidth,'Linestyle','-','Color',c(1,:))
plot(LoadPointD_results_slip.*1e6 + 400,A_series_slip*100,'linewidth',slip_linewidth,'Linestyle','-','Color',c(2,:))
hold off

ytickLoc = [0.3:0.05:0.5];
ytickLabel = num2cell(ytickLoc);
ytickLabel = string(ytickLabel);
ytickLabel(2:2:end) = nan;
yticks(ytickLoc)
yticklabels(ytickLabel)

set(gca,'Ycolor','k');
%xlabel('Displacement, \mu m')
ylabel('area percentage')
ylim([0.3 0.8])
grid on

% legend_str = {'lab data','Aging model','Nagata-aging model'};
% legend(p,legend_str,'location','best')

xlim([400 800])
set(gca,'Fontsize',22,'Fontweight','bold')
set(gca, 'FontName', 'Helvetica')
set(gcf, 'Renderer', 'Painters');% make eps clear


%% plot Ar fit vs residue, using loadpoint velocity, just aging law
clear p legend_str

load ColorScheme_MatlabOrder.mat
ImportDK94_friction;
ImportDK94_PA

data_linewidth = 6;
aging_linewidth = 4;
slip_linewidth = 3;

data_color = [0.7 0.7 0.7];


D_DK94Friction = DK94Fig8DataDigitizedFriction.D_DK94Friction;
f_DK94Friction = DK94Fig8DataDigitizedFriction.f_DK94Friction;
D_DK94PA = DK94Fig8DataDigitizedPercentArea.D_DK94PA;
PA_DK94PA = DK94Fig8DataDigitizedPercentArea.PA_DK94PA;

figure(7)


p(1) = plot(D_DK94PA,PA_DK94PA,...
    'linewidth',data_linewidth,'Linestyle','-','Color',data_color);
hold on
p(2) = plot(LoadPointD_results_aging.*1e6 + 400,A_series_aging*100,'linewidth',aging_linewidth,'Linestyle','-','Color',c(1,:));

%CalculateResidue
[C_temp,ia,ic] = unique(LoadPointD_results_aging);%remove redundent x points
ModelAr_AtDataD_aging = interp1(LoadPointD_results_aging(ia).*1e6 + 400,A_series_aging(ia)*100,D_DK94PA);

p(3) = plot(D_DK94PA,ModelAr_AtDataD_aging-PA_DK94PA,'linewidth',aging_linewidth,'Linestyle','-','Color','k');
hold off

ytickLoc = [-0.1:0.05:0.5];
ytickLabel = num2cell(ytickLoc);
ytickLabel = string(ytickLabel);
ytickLabel(2:2:end) = nan;
yticks(ytickLoc)
yticklabels(ytickLabel)

set(gca,'Ycolor','k');
%xlabel('Displacement, \mu m')
ylabel('area percentage')
ylim([-0.1 0.5])
grid on

legend_str = {'lab data','aging law model','residue'};
legend(p,legend_str,'location','best');

xlim([400 800])
set(gca,'Fontsize',22,'Fontweight','bold')
set(gca, 'FontName', 'Helvetica')
set(gcf, 'Renderer', 'Painters');% make eps clear


%% plot Ar fit vs residue, using loadpoint velocity, just slip law
clear p legend_str

load ColorScheme_MatlabOrder.mat
ImportDK94_friction;
ImportDK94_PA

data_linewidth = 6;
aging_linewidth = 4;
slip_linewidth = 4;

data_color = [0.7 0.7 0.7];


D_DK94Friction = DK94Fig8DataDigitizedFriction.D_DK94Friction;
f_DK94Friction = DK94Fig8DataDigitizedFriction.f_DK94Friction;
D_DK94PA = DK94Fig8DataDigitizedPercentArea.D_DK94PA;
PA_DK94PA = DK94Fig8DataDigitizedPercentArea.PA_DK94PA;

figure(8)


p(1) = plot(D_DK94PA,PA_DK94PA,...
    'linewidth',data_linewidth,'Linestyle','-','Color',data_color);
hold on
p(2) = plot(LoadPointD_results_slip.*1e6 + 400,A_series_slip*100,'linewidth',slip_linewidth,'Linestyle','-','Color',c(2,:));

%CalculateResidue
[C_temp,ia,ic] = unique(LoadPointD_results_slip);%remove redundent x points
ModelAr_AtDataD_slip = interp1(LoadPointD_results_slip(ia).*1e6 + 400,A_series_slip(ia)*100,D_DK94PA);

p(3) = plot(D_DK94PA,ModelAr_AtDataD_slip-PA_DK94PA,'linewidth',slip_linewidth,'Linestyle','-','Color','k');
hold off

ytickLoc = [-0.1:0.05:0.5];
ytickLabel = num2cell(ytickLoc);
ytickLabel = string(ytickLabel);
ytickLabel(2:2:end) = nan;
yticks(ytickLoc)
yticklabels(ytickLabel)

set(gca,'Ycolor','k');
%xlabel('Displacement, \mu m')
ylabel('area percentage')
ylim([-0.1 0.5])
grid on

legend_str = {'lab data','slip law model','residue'};
legend(p,legend_str,'location','best')

xlim([400 800])
set(gca,'Fontsize',22,'Fontweight','bold')
set(gca, 'FontName', 'Helvetica')
set(gcf, 'Renderer', 'Painters');% make eps clear


%% plot model against DK94 fig8 data, using loadpoint velocity, friction, both aging and slip law



clear p legend_str

load ColorScheme_MatlabOrder.mat
ImportDK94_friction;
ImportDK94_PA

data_linewidth = 3;
aging_linewidth = 2;
slip_linewidth = 1.5;

data_color = [0.7 0.7 0.7];


D_DK94Friction = DK94Fig8DataDigitizedFriction.D_DK94Friction;
f_DK94Friction = DK94Fig8DataDigitizedFriction.f_DK94Friction;
D_DK94PA = DK94Fig8DataDigitizedPercentArea.D_DK94PA;
PA_DK94PA = DK94Fig8DataDigitizedPercentArea.PA_DK94PA;

figure(11)

%yyaxis left
p(1) = plot(D_DK94Friction,f_DK94Friction,...
    'linewidth',data_linewidth,'Linestyle','-','Color',data_color);
hold on

p(2) = plot(LoadPointD_results_aging.*1e6 + 400,f_series_aging,...
    'linewidth',aging_linewidth,'Linestyle',':','Color',c(1,:));

p(3) = plot(LoadPointD_results_slip.*1e6 + 400,f_series_slip,...
    'linewidth',slip_linewidth,'Linestyle','-','Color',c(2,:));

hold off

ytickLoc = [0.55:0.05:0.9];
ytickLabel = num2cell(ytickLoc);
ytickLabel = string(ytickLabel);
ytickLabel(1:2:end) = nan;
yticks(ytickLoc)
yticklabels(ytickLabel)

set(gca,'Ycolor','k');
xlabel('Loadpoint Displacement, \mu m')
ylabel('Friction')
ylim([0.65 0.9])

pbaspect([4.2 1 1])

grid on

legend_str = {'lab data','Aging model','Slip model'};
legend(p,legend_str,'location','best','orientation','horizontal')

xlim([400 800])
set(gca,'Fontsize',15,'Fontweight','bold')
set(gca,'Linewidth',1,'Fontweight','bold')
set(gca, 'FontName', 'Helvetica')
set(gcf, 'Renderer', 'Painters');% make eps clear

%% plot model against DK94 fig8 data, using loadpoint velocity, contact area, both aging and slip law

clear p legend_str

load ColorScheme_MatlabOrder.mat
ImportDK94_friction;
ImportDK94_PA

data_linewidth = 3;
aging_linewidth = 2;
slip_linewidth = 1.5;

data_color = [0.7 0.7 0.7];


D_DK94Friction = DK94Fig8DataDigitizedFriction.D_DK94Friction;
f_DK94Friction = DK94Fig8DataDigitizedFriction.f_DK94Friction;
D_DK94PA = DK94Fig8DataDigitizedPercentArea.D_DK94PA;
PA_DK94PA = DK94Fig8DataDigitizedPercentArea.PA_DK94PA;

figure(12)

%yyaxis left
p(1) = plot(D_DK94PA,PA_DK94PA,...
    'linewidth',data_linewidth,'Linestyle','-','Color',data_color);
hold on

p(2) = plot(LoadPointD_results_aging.*1e6 + 400,A_series_aging*100,'linewidth',aging_linewidth,'Linestyle',':','Color',c(1,:));

p(3) = plot(LoadPointD_results_slip.*1e6 + 400,A_series_slip*100,'linewidth',slip_linewidth,'Linestyle','-','Color',c(2,:));

hold off

ytickLoc = [0.25:0.05:0.50];
ytickLabel = num2cell(ytickLoc);
ytickLabel = string(ytickLabel);
ytickLabel(1:2:end) = nan;
yticks(ytickLoc)
yticklabels(ytickLabel)

set(gca,'Ycolor','k');
xlabel('Loadpoint Displacement, \mu m')
ylabel('Contact Area (%)')
ylim([0.3 0.5])

pbaspect([4.2 1 1])

grid on

legend_str = {'lab data','Aging model','Slip model'};
legend(p,legend_str,'location','best','orientation','horizontal')

%set(gca,'TickDir','out');

xlim([400 800])
set(gca,'Fontsize',15,'Fontweight','bold')
set(gca,'Linewidth',1,'Fontweight','bold')
set(gca, 'FontName', 'Helvetica')
set(gcf, 'Renderer', 'Painters');% make eps clear



