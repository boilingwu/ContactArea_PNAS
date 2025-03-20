% find the slip rate as a vector given the external stress vector

%%%%%%%%%%%%% Inputs %%%%%%%%%%%%%%%%%%
% tau_ext, external shear stress at each element at this moment, vector
% Theta_now, state variable at each element at this moment, vector
% G, shear modulus, scalar
% Cs, shear wave speed, scalar
% a, rate-state a, vector
% b, rate-state b, vector
% mu0, rate-state mu0, vector
% V0, rate-state V0, vector
% Theta0, rate-state Theta0, vector
% sigma, normal stress at each element at this moment, vector

% FricOnly_threRatio, the ratio threshold to determing if
% frictional-velocity dependent is dominant, scalar

% RDOnly_threRatio, the ratio threshold to determing if
% radiation damping-velocity dependent is dominant, scalar

% GradientMethodExit_DiffThreRatio, the ratio threshold for
% Newton approach, scalar

% iter_Threhold_GradientMethod, the ratio threshold for
% Newton approach, scalar

%%%%%%%%%%%%% Output %%%%%%%%%%%%%%%%%%
% speed_now, speed under this condition, vector 



function [speed_now,GoToDyna_Flag] = SolveForV_QuasiStatic(tau_ext,Theta_now,...
    G,Cs,a,b,mu0,V0,Theta0,sigma,...
    FricOnly_threRatio,RDOnly_threRatio)
    
%     size_tau_ext = size(tau_ext);
%     size_Theta_now = size(Theta_now);
%     size_sigma = size(sigma);
    if ~isequaln(size(tau_ext),size(Theta_now),size(sigma),size(a),...
            size(b),size(mu0),size(V0),size(Theta0))
        error('Input arrays dimenstions do not match Theta_now each other!')
    end

    % temporal variable
%     Speed_IfFricOnly = nan;
%     RDamping_IfFricOnly = nan;
%     Speed_IfRDOnly = nan;
%     RSFricStress_IfRDOnly = nan;
%     Speed_BothRDFric_trial = nan;
%     Stress_BothRDFric_trial = nan;
%     DiffStress_BothRDFric_trial = nan;
%     slope_trialnow = nan;
%     iter_count_GradientMethod = nan;
%     max_speed = nan;

    % speed that if not considering radiation dampping
    Speed_IfFricOnly = V0.*exp(tau_ext./(a.*sigma)-mu0./a).*(Theta_now./Theta0).^(-b./a); 
    RDamping_IfFricOnly = G/(2*Cs).*Speed_IfFricOnly;
    speed_now = Speed_IfFricOnly;
    %if external stress now is significantly larger than the RDamping_FricOnly,
    %than we can ignore RDamping effect  
    if max(abs(RDamping_IfFricOnly./tau_ext)) < FricOnly_threRatio        
        GoToDyna_Flag = 0;
    else
        GoToDyna_Flag = 1;
    end     
end