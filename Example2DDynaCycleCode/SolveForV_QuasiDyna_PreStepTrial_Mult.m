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



function speed_now = SolveForV_QuasiDyna_PreStepTrial_Mult(tau_ext,Theta_now,...
    G,Cs,a,b,mu0,V0,Theta0,sigma,...
    SR_prev,...
    FricOnly_threRatio,RDOnly_threRatio,...
    GradientMethodExit_DiffThreRatio,iter_Threhold_GradientMethod)
    
%     size_tau_ext = size(tau_ext);
%     size_Theta_now = size(Theta_now);
%     size_sigma = size(sigma);
    if ~isequaln(size(tau_ext),size(Theta_now),size(sigma),size(a),...
            size(b),size(mu0),size(V0),size(Theta0))
        error('Input arrays dimenstions do not match Theta_now each other!')
    end

    Speed_BothRDFric_trial = SR_prev;
    Speed_BothRDFric_trial(Speed_BothRDFric_trial<0)=10^(-12);
    
    % friction coefficient
    friction=mu0.*(Speed_BothRDFric_trial./V0).^(a./mu0).*(Theta_now./Theta0).^(b./mu0);
    
    % multiplicative, power-law form
    Stress_BothRDFric_trial = G/(2*Cs).*Speed_BothRDFric_trial + ...
        friction.*sigma;
    DiffStress_BothRDFric_trial = tau_ext - Stress_BothRDFric_trial;
    
    
    Error_DiffStress_ratio = max(abs(DiffStress_BothRDFric_trial./tau_ext));
    iter_count_GradientMethod = 0;

    Error_DiffStress_ratio_pre = 9999999999;
    while(1==1)
        if Error_DiffStress_ratio < GradientMethodExit_DiffThreRatio
           speed_now = Speed_BothRDFric_trial;
           break
        elseif abs((Error_DiffStress_ratio_pre-Error_DiffStress_ratio)./Error_DiffStress_ratio)...
                < 1e-7
           speed_now = Speed_BothRDFric_trial;
           break                    
        end
        iter_count_GradientMethod = iter_count_GradientMethod + 1;
        if iter_count_GradientMethod > iter_Threhold_GradientMethod
            error('Too many iteration, abort!')
        end                 
        slope_trialnow = G/(2*Cs) + a.*sigma.*(friction./mu0)./Speed_BothRDFric_trial;
        Speed_BothRDFric_trial = Speed_BothRDFric_trial + ...
        DiffStress_BothRDFric_trial./slope_trialnow;    
    
        Speed_BothRDFric_trial(Speed_BothRDFric_trial<0)=10^(-12);
    
        
    
        %update friction and stress
        friction=mu0.*(Speed_BothRDFric_trial./V0).^(a./mu0).*(Theta_now./Theta0).^(b./mu0);
        Stress_BothRDFric_trial = G/(2*Cs).*Speed_BothRDFric_trial + ...
            friction.*sigma;
        DiffStress_BothRDFric_trial = tau_ext - Stress_BothRDFric_trial;    
    

        Error_DiffStress_ratio = max(abs(DiffStress_BothRDFric_trial./tau_ext));
        Error_DiffStress_ratio_pre = Error_DiffStress_ratio;
    end

   
end