% function for ode, dy/dt = f(t,y)

function dy_dt = RSmarching_RD_mult_1D_SlipNagata(t,y,...
                Kshear_matrix,Knormal_matrix,...
                Dc,G,Cs,a,b,... 
                mu0,V0,Theta0,...% All input vectors in this line is 1*N (need transpose)
                v_load,...
                c)
    % if there is N element, the total number of row(1 colomn) of y is 5N
    % y(1) - y(N): external shear stress 
    % y(N+1) - y(2N): normal stress
    % y(2N+1) - y(3N): slip rate
    % y(3N+1) - y(4N): slip
    % y(4N+1) - y(5N): theta
    

    
    if mod(length(y),5) ~= 0
        error('matrix dimention not right')
    end
    
    N = length(y)/5;
    
    y1 = y(1:N,1);
    y2 = y(N+1:2*N,1);
    y3 = y(2*N+1:3*N,1);
    y4 = y(3*N+1:4*N,1);
    y5 = y(4*N+1:5*N,1);
    
    tau_ext = y1;
    sig = y2;
    v = V0'.*exp(y3);
    D = y4;
    theta = Theta0'.*exp(y5);
    
    dtau_ext_dt = Kshear_matrix*(v-v_load);
    dsig_dt = Knormal_matrix*(v-v_load);
    dD_dt = v;
    
    % typical aging law
    %dy5_dt = 1./theta-v./Dc';
    
    % Nagata Law + aging law
    %dy5_dt = 1./theta-v./Dc'-(c./b.*dtau_ext_dt./sig); 

     % Nagata Law + slip law
    dy5_dt = -v./Dc'.*log(theta.*v./Dc')-(c./b.*dtau_ext_dt./sig);    
    
    %%% below is wrong %%%%
%     dtau_ext_dt = (Kshear_matrix*v+tau_rate'+...
%         G/(2*Cs).*(v./a').*((tau_ext.*dsig_dt)./(sig.*sig)+b'.*dtheta_dt./theta))...
%         ./(1+G/(2*Cs).*(v./a')./sig);   
    %consider radiation dampping, have to solve to equation simultaneously
    %dv_dt = (v./a').*(dtau_dt./sig-(tau.*dsig_dt)./(sig.*sig)-b'.*dtheta_dt./theta);
    %dy3_dt = (b'./a'.*dy5_dt+dsig_dt./(a'.*sig.^2).*(tau_ext-G/(2*Cs).*v)-1./(a'.*sig).*dtau_ext_dt)./(1-1./(a'.*sig).*G/(2*Cs).*v);
    %%% above is wrong %%%%
    
    
    %%% below is logarithm form %%%%
    %dy3_dt =
    %(dtau_ext_dt-(tau_ext-G/(2*Cs).*v).*dsig_dt./sig.^2-b'.*sig.*dy5_dt)./(a'.*sig+G/(2*Cs).*v);
    % wrong, sig.^2 ---> sig
    %dy3_dt = (dtau_ext_dt-(tau_ext-G/(2*Cs).*v).*dsig_dt./sig-b'.*sig.*dy5_dt)./(a'.*sig+G/(2*Cs).*v); 
    
    
    %RD_Factor = G/(2*Cs);
    RD_Factor = 0;
    
    %%% below is multiplicative form %%%%
    friction = (tau_ext-RD_Factor.*v)./sig;
    dy3_dt = (dtau_ext_dt-friction.*dsig_dt-b'.*sig.*dy5_dt.*friction./mu0')./(a'.*sig.*friction./mu0'+RD_Factor.*v); 
    
    dy_dt(1:N,1) = dtau_ext_dt;
    dy_dt(N+1:2*N,1) = dsig_dt;
    dy_dt(2*N+1:3*N,1) = dy3_dt;
    dy_dt(3*N+1:4*N,1) = dD_dt;
    dy_dt(4*N+1:5*N,1) = dy5_dt;
    
    
end