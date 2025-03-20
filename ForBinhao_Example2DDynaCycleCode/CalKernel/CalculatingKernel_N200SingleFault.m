% 2D dynamic rupture code for slip weakening law, a first very simple model
clear all
clc

%% material setup
ct = 1361;
cl = 2680;
rho = 1170;
mu = rho*ct*ct;
lambda = rho*cl*cl-2*mu;
nu = (cl^2-2*ct^2)/(2*(cl^2-ct^2));

%% Discretisaztion setup
dx = 1e-3;%mm
discretization_size = dx;
delta_s = dx;
CFL = 0.45;
dt = CFL*dx/cl;
delta_t = dt;

%% Fault geometry setup a simple step over
% making geometry, all fault along x3
%fix parameter

%% %%%%%%%%%%%%%%% Fault No1, geometry %%%%%%%%%%%%%%%%%%%%
ind = 1;
Fault(ind).NumOfElement = 200;
Fault(ind).Length = Fault(ind).NumOfElement*discretization_size;
Fault(ind).x1_center_start = 0;
Fault(ind).x1_center_end = ...
    Fault(ind).x1_center_start+(Fault(ind).NumOfElement-1)*discretization_size;
Fault(ind).x2 = 0;
Fault(ind).n = [0 1]';
% Fault(ind).Vl = 1e-9;
% Fault(ind).dip = 90;
% Fault(ind).rake = 0;% the last role read in, Sylvain call it the rake

% allocate array
Fault(ind).Element_No = [1:Fault(ind).NumOfElement]';%
ElementOnesArray = Fault(ind).Element_No./Fault(ind).Element_No;

Fault(ind).Element_x1 = ...
    [Fault(ind).x1_center_start:discretization_size:Fault(ind).x1_center_start+(Fault(ind).NumOfElement-1)*discretization_size]';%
Fault(ind).Element_x2 = Fault(ind).x2 * ElementOnesArray;%
Fault(ind).Element_width = discretization_size * ElementOnesArray;%
Fault(ind).Element_n = cat(1,ElementOnesArray'*Fault(ind).n(1),ElementOnesArray'*Fault(ind).n(2));%
Fault(ind).Element_n = Fault(ind).Element_n';
Fault(ind).Element_FaultNo = ind * ElementOnesArray;%

% frictional setting
% Fault(ind).tau0 = 3e6; %MPa
% Fault(ind).mu0 = 0.6;
% Fault(ind).mu_d = 0.1;
% Fault(ind).sig0 = 10;
% Fault(ind).D0 = 0.06;
% 
% 
% Fault(ind).Element_tau0 = Fault(ind).tau0 * ElementOnesArray;%
% Fault(ind).Element_mu0 = Fault(ind).mu0 * ElementOnesArray;%
% Fault(ind).Element_sig0 = Fault(ind).sig0 * ElementOnesArray;%
% Fault(ind).Element_mu_d = Fault(ind).mu_d * ElementOnesArray;%
% Fault(ind).Element_D0 = Fault(ind).D0 * ElementOnesArray;%


%% Making array for kernel calculation

% element_array
element_location_x1 = [];
element_location_x2 = [];
element_location_n = [];
for ind = 1:length(Fault)
    element_location_x1 = cat(1,element_location_x1,Fault(ind).Element_x1);
    element_location_x2 = cat(1,element_location_x2,Fault(ind).Element_x2);
    element_location_n = cat(1,element_location_n,Fault(ind).Element_n);
end

%% Calculate Time range that needed to be calculated in the kernel
clear ElementCut MaxKernelLength_AllElem_arr
Buffer_Timestep = 500;
%The number of timesteps I assumed that S wave passed and stress decay to thestatic level
for elem_ind = 1:length(element_location_x1)
    
    receiver_x1 = element_location_x1 - element_location_x1(elem_ind);
    receiver_x2 = element_location_x2 - element_location_x2(elem_ind);
    
    ElementCut(elem_ind).r = sqrt(receiver_x1.^2 + receiver_x2.^2);
    ElementCut(elem_ind).shortest_arrivalTime = ElementCut(elem_ind).r./cl - delta_s/cl;
    ElementCut(elem_ind).longest_arrivalTime = ElementCut(elem_ind).r./ct + delta_s/ct;
    
    ElementCut(elem_ind).StartIndex = floor(ElementCut(elem_ind).shortest_arrivalTime./dt) - 1; 
    % minus 1 to play safe.
    ElementCut(elem_ind).StartIndex(ElementCut(elem_ind).StartIndex < 1) = 1;
    ElementCut(elem_ind).EndIndex = ceil(ElementCut(elem_ind).longest_arrivalTime./dt) + 1 + Buffer_Timestep; 
    % plus 1 to play safe
    ElementCut(elem_ind).KernelLength = ElementCut(elem_ind).EndIndex - ElementCut(elem_ind).StartIndex;
    ElementCut(elem_ind).MaxKernelLength = max(ElementCut(elem_ind).KernelLength);
    ElementCut(elem_ind).MaxEndIndex = max(ElementCut(elem_ind).EndIndex);
    MaxKernelLength_AllElem_arr(elem_ind) = ElementCut(elem_ind).MaxKernelLength;
    MaxEndIndex_AllElem_arr(elem_ind) = ElementCut(elem_ind).MaxEndIndex;
    %ElementCut(elem_ind).EndIndex(ElementCut(elem_ind).EndIndex > MAX_TIMESTEP) = MAX_TIMESTEP;
end
MaxKernelLength_AllElem = max(MaxKernelLength_AllElem_arr);
MaxEndIndex_AllElem = max(MaxEndIndex_AllElem_arr(elem_ind));

% use different StartIndex
% t_MemmoryArray = (dt:dt:MaxKernelLength_AllElem*dt)+dt;

% All StartIndex from 1
t_MemmoryArray = (dt:dt:MaxEndIndex_AllElem*dt)+dt;

KernelInfo.MaxEndIndex_AllElem = MaxEndIndex_AllElem;
KernelInfo.t_MemmoryArray = t_MemmoryArray;

% note that t(1) =2dt,represent the time range from dt -> 2dt.
% the time range from 0 to dt is not considered here, (which is radiation damping factor, or instantaneous response)
% Baoning's note on March 20, 2022: I may not think carefully enough about
% the time shift stuff. But I think I have played safe here and below

%% Try to compute the kernel_matrix for each elements
clear ElementKernel

%allocate the matrix-version kernel matrix
KernelRowLength = length(element_location_x1)*length(t_MemmoryArray);
Kernel_K11_all = nan(length(element_location_x1),KernelRowLength);
Kernel_K12_all = nan(length(element_location_x1),KernelRowLength);
Kernel_K22_all = nan(length(element_location_x1),KernelRowLength);

%parpool(4)
parfor elem_ind = 1:length(element_location_x1)
    elem_ind
    receiver_x1 = element_location_x1 - element_location_x1(elem_ind);
    receiver_x2 = element_location_x2 - element_location_x2(elem_ind);
    
    [receiver_x1_mat,t_mat1] = meshgrid(receiver_x1,t_MemmoryArray);
    [receiver_x2_mat,t_mat2] = meshgrid(receiver_x2,t_MemmoryArray);
    if isequal(t_mat1,t_mat2)
        t_mat = t_mat1;
    else
        error('something about the meshing is wrong!')
    end
    
    StartIndex_thisElement = ElementCut(elem_ind).StartIndex';%the outcome should be 1*N
    % shift t matrix
    %t_mat = t_mat + StartIndex_thisElement*dt;
    
    K12 = I12_matrixBothXandT(receiver_x1_mat+0.5*delta_s,receiver_x2_mat,t_mat,ct,cl)...
        -I12_matrixBothXandT(receiver_x1_mat-0.5*delta_s,receiver_x2_mat,t_mat,ct,cl)...
        -I12_matrixBothXandT(receiver_x1_mat+0.5*delta_s,receiver_x2_mat,t_mat-delta_t,ct,cl)...
        +I12_matrixBothXandT(receiver_x1_mat-0.5*delta_s,receiver_x2_mat,t_mat-delta_t,ct,cl);
    %ElementKernel(elem_ind).K12 = K12*(-mu/(2*ct));
    Kernel_K12_all(elem_ind,:) = reshape(K12,[1 KernelRowLength])*(-mu/(2*ct));
    
    K11 = I11_matrixBothXandT(receiver_x1_mat+0.5*delta_s,receiver_x2_mat,t_mat,ct,cl)...
        -I11_matrixBothXandT(receiver_x1_mat-0.5*delta_s,receiver_x2_mat,t_mat,ct,cl)...
        -I11_matrixBothXandT(receiver_x1_mat+0.5*delta_s,receiver_x2_mat,t_mat-delta_t,ct,cl)...
        +I11_matrixBothXandT(receiver_x1_mat-0.5*delta_s,receiver_x2_mat,t_mat-delta_t,ct,cl);
    %ElementKernel(elem_ind).K11 = K11*(-mu/(2*ct));
    Kernel_K11_all(elem_ind,:) = reshape(K11,[1 KernelRowLength])*(-mu/(2*ct));
    
    K22 = I22_matrixBothXandT(receiver_x1_mat+0.5*delta_s,receiver_x2_mat,t_mat,ct,cl)...
        -I22_matrixBothXandT(receiver_x1_mat-0.5*delta_s,receiver_x2_mat,t_mat,ct,cl)...
        -I22_matrixBothXandT(receiver_x1_mat+0.5*delta_s,receiver_x2_mat,t_mat-delta_t,ct,cl)...
        +I22_matrixBothXandT(receiver_x1_mat-0.5*delta_s,receiver_x2_mat,t_mat-delta_t,ct,cl);
    %ElementKernel(elem_ind).K22 = K22*(-mu/(2*ct));
    Kernel_K22_all(elem_ind,:) = reshape(K22,[1 KernelRowLength])*(-mu/(2*ct));
end

%% save variable
%save(['ElementKernel_N',num2str(length(element_location_x1)),'.mat'],'ElementKernel','-v7.3');
save(['ElementCut_N',num2str(length(element_location_x1)),'_SingleFault.mat'],'ElementCut','-v7.3')
save(['KernelInfo_N',num2str(length(element_location_x1)),'_SingleFault.mat'],'KernelInfo','-v7.3')
save(['Fault_N',num2str(length(element_location_x1)),'_SingleFault.mat'],'Fault','-v7.3')

save(['MatrixKernel_N',num2str(length(element_location_x1)),'_SingleFault.mat'],'Kernel_K11_all','Kernel_K12_all','Kernel_K22_all','-v7.3')



%%
%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function answer = I12_matrixBothXandT(x1,x2,t,ct,cl)

% x1, x2, and t need to be three 2D matrix of the same size
% the 1st dimension of the matrix represent location, Let's say there are N
% element
% the 2nd dimension of the matrix represent time

    %%%% Preprocessing, checking if the dimension is correct
    
    if ~isequal(size(x1),size(x2),size(t))
        error(' one of the input array has dimension that is not correct')
    end

    r = sqrt(x1.^2+x2.^2);
    p = ct/cl;
    AlphaT = (cl*t./r).^2-1;
    BetaT = (ct*t./r).^2-1;
    
    answer = t.*nan;
%     if t-r/cl <=0
%         answer = heaviside_BW(x1)*heaviside_BW(t-abs(x2)/ct);
%     
%     elseif t-r/cl > 0 && t-r/ct <= 0
%         answer = heaviside_BW(x1)*heaviside_BW(t-abs(x2)/ct) + ...
%            -1/pi*2*x1/r*p*(2*(3*x2^2-x1^2)/(3*r^2)*p^2*AlphaT^(3/2)+2*x2^2/r^2*p^2*AlphaT^(1/2));
%     
%     elseif t-r/ct > 0
%            answer = heaviside_BW(x1)*heaviside_BW(t-abs(x2)/ct) + ...
%            -1/pi*2*x1/r*p*(2*(3*x2^2-x1^2)/(3*r^2)*p^2*AlphaT^(3/2)+2*x2^2/r^2*p^2*AlphaT^(1/2))...
%            +1/pi*sgn_BW(x1)*(2*abs(x1)/r*(2*(3*x2^2-x1^2)/(3*r^2)*BetaT^(3/2)+2*x2^2/r^2*BetaT^(1/2))-acos(abs(x1)/sqrt(ct^2*t^2-x2^2)));
%     else
%         error('something wrong')
%     end


      cdt_now = t-r/cl <=0;
      answer(cdt_now) = heaviside_BW(x1(cdt_now)).*heaviside_BW(t(cdt_now)-abs(x2(cdt_now))/ct);
      
      cdt_now = t-r/cl > 0 & t-r/ct <= 0;
      answer(cdt_now) = heaviside_BW(x1(cdt_now))'.*heaviside_BW(t(cdt_now)-abs(x2(cdt_now))/ct)' + ...
            -1/pi*2*x1(cdt_now)./r(cdt_now)*p.*(2*(3*x2(cdt_now).^2-x1(cdt_now).^2)./(3*r(cdt_now).^2)*p^2.*AlphaT(cdt_now).^(3/2)+2*x2(cdt_now).^2./r(cdt_now).^2*p^2.*AlphaT(cdt_now).^(1/2));
       
      cdt_now = t-r/ct > 0;  
      answer(cdt_now) = heaviside_BW(x1(cdt_now))'.*heaviside_BW(t(cdt_now)-abs(x2(cdt_now))/ct)' + ...
       -1/pi*2*x1(cdt_now)./r(cdt_now)*p.*(2*(3*x2(cdt_now).^2-x1(cdt_now).^2)./(3*r(cdt_now).^2)*p^2.*AlphaT(cdt_now).^(3/2)+2*x2(cdt_now).^2./r(cdt_now).^2*p^2.*AlphaT(cdt_now).^(1/2))...
       +1/pi*sgn_BW(x1(cdt_now))'.*(2.*abs(x1(cdt_now))./r(cdt_now).*(2*(3*x2(cdt_now).^2-x1(cdt_now).^2)./(3*r(cdt_now).^2).*BetaT(cdt_now).^(3/2)+2*x2(cdt_now).^2./r(cdt_now).^2.*BetaT(cdt_now).^(1/2))-acos(abs(x1(cdt_now))./sqrt(ct^2.*t(cdt_now).^2-x2(cdt_now).^2)));
      

end

function answer = I11_matrixBothXandT(x1,x2,t,ct,cl)

    if ~isequal(size(x1),size(x2),size(t))
        error(' one of the input array has dimension that is not correct')
    end
    
    r = sqrt(x1.^2+x2.^2);
    p = ct/cl;
    AlphaT = (cl*t./r).^2-1;
    BetaT = (ct*t./r).^2-1;
    
    answer = t.*nan;
    
    cdt_now = t-r/cl <=0;
    answer(cdt_now) = t(cdt_now).*0;
    
    cdt_now = t-r/cl > 0 & t-r/ct <= 0;
    answer(cdt_now) = (-1/pi)*2*(x2(cdt_now)./r(cdt_now))*p.*...
        (2*(3*(x1(cdt_now).^2)-(x2(cdt_now).^2))./(3*(r(cdt_now).^2))*(p.^2).*(AlphaT(cdt_now).^(3.0/2))...
        + (1-2*(x2(cdt_now).^2./r(cdt_now).^2)*p.^2).*(AlphaT(cdt_now).^(1.0/2)));
    
    cdt_now = t-r/ct > 0;
    answer(cdt_now) = (-1/pi)*2*(x2(cdt_now)./r(cdt_now))*p.*...
        (2*(3*(x1(cdt_now).^2)-(x2(cdt_now).^2))./(3*(r(cdt_now).^2))*(p.^2).*(AlphaT(cdt_now).^(3.0/2))...
        + (1-2*(x2(cdt_now).^2./r(cdt_now).^2)*p.^2).*(AlphaT(cdt_now).^(1.0/2)))...
        + (1/pi)*2*(x2(cdt_now)./r(cdt_now)).*...
        (2*(3*(x1(cdt_now).^2)-(x2(cdt_now).^2))./(3*(r(cdt_now).^2)).*(BetaT(cdt_now).^(3.0/2))...
        +(1-2*(x2(cdt_now).^2./r(cdt_now).^2)).*(BetaT(cdt_now).^(1.0/2)));

end

function answer = I22_matrixBothXandT(x1,x2,t,ct,cl)
    
    if ~isequal(size(x1),size(x2),size(t))
        error(' one of the input array has dimension that is not correct')
    end

    r = sqrt(x1.^2+x2.^2);
    p = ct/cl;
    AlphaT = (cl*t./r).^2-1;
    BetaT = (ct*t./r).^2-1;
    
    answer = t.*nan;
    
    cdt_now = t-r/cl <=0;
    answer(cdt_now) = t(cdt_now).*0;
    
    cdt_now = t-r/cl > 0 & t-r/ct <= 0;
    answer(cdt_now) = (1/pi)*2*(x2(cdt_now)./r(cdt_now))*p.*...
        (2*(3*(x1(cdt_now).^2)-(x2(cdt_now).^2))./(3*(r(cdt_now).^2))*(p.^2).*(AlphaT(cdt_now).^(3.0/2))...
        + (2*(x1(cdt_now).^2./r(cdt_now).^2)*p.^2-1).*(AlphaT(cdt_now).^(1.0/2)));
   
    cdt_now = t-r/ct > 0;
    answer(cdt_now) = (1/pi)*2*(x2(cdt_now)./r(cdt_now))*p.*...
        (2*(3*(x1(cdt_now).^2)-(x2(cdt_now).^2))./(3*(r(cdt_now).^2))*(p.^2).*(AlphaT(cdt_now).^(3.0/2))...
        + (2*(x1(cdt_now).^2./r(cdt_now).^2)*p.^2-1).*(AlphaT(cdt_now).^(1.0/2)))...
        + (-1/pi)*2*(x2(cdt_now)./r(cdt_now)).*...
        (2*(3*(x1(cdt_now).^2)-(x2(cdt_now).^2))./(3*(r(cdt_now).^2)).*(BetaT(cdt_now).^(3.0/2))...
        + (2*(x1(cdt_now).^2./r(cdt_now).^2)-1).*(BetaT(cdt_now).^(1.0/2)));

end

% function ans_h = heaviside_BW(x)
% %     if x>0
% %         ans = 1;
% %     else 
% %         ans = 0;
% %     end
%     ans_h(x>0) = 1;
%     ans_h(x<=0) = 0;
% end
% 
% function ans_s = sgn_BW(x)
%     ans_s(x>=0) = 1;
%     ans_s(x<0) = -1;
% end