%Function to calculate the custom external force
%Last modified by Anup Teejo Mathew 30/06/2021

function Fext=CustomExtForce(S,q,g,J,t,qd,eta,Jdot)
%%%%NOTE%%%%


%S: Linkage element, 
%q and qd: joint coordinates and their time derivatives, 
%g, J, Jd, and eta: transformation matrix, Jacobian, time derivative of jacobian, and screw velocity at every significant point of the linkage
%t:  time

%Fext should be 6*n column vector where n is the total number of gaussian points of all soft links (nGauss) + number of rigid links.
%(Example: linkage with 2 soft links and 1 rigid link (n=nGauss1+nGauss2+1)
%Fext should be arranged according to the order of precedence
%Fext should be distributed force for a soft link and point force for a rigid link

% Significant points: 1 for every joint, 1 at the center of the rigid link, 1 at the start and end of every soft link and 1 for each Gaussian points

% J   = S.Jacobian(q);         %geometric jacobian of the linkage calculated at every significant points
% g   = S.FwdKinematics(q);    %transformation matrix of the linkage calculated at every significant points
% eta = S.ScrewVelocity(q,qd); %Screwvelocity of the linkage calculated at every significant points
% J   = S.Jacobiandot(q,qd);   %time derivative of geometric jacobian of the linkage calculated at every significant points

%%%END%%%

DL     = S.CP1;
G      = S.G;

i_sig    = 1;
i_sig_nj = 1;


%body
i_sig    = i_sig+1;%joint

rho_water  = 1000;

rho_body   = S.VLinks(1).Rho;
M_here     = S.VLinks(1).Ms;
g_here     = g((i_sig-1)*4+1:i_sig*4,:);
eta_here   = eta((i_sig-1)*6+1:i_sig*6);
DL_here    = DL((i_sig_nj-1)*6+1:i_sig_nj*6,:);
Fextb      = -rho_water/rho_body*M_here*dinamico_Adjoint(ginv(g_here))*G-DL_here*norm(eta_here(4:6))*eta_here;%-rho_water/rho_body

i_sig    = i_sig+1;
i_sig_nj = i_sig_nj+1;

%shaft
i_sig    = i_sig+1;%joint

rho_shaft  = S.VLinks(2).Rho;
M_here     = S.VLinks(2).Ms;
g_here     = g((i_sig-1)*4+1:i_sig*4,:);
eta_here   = eta((i_sig-1)*6+1:i_sig*6);
DL_here    = DL((i_sig_nj-1)*6+1:i_sig_nj*6,:);
Fexts      = -rho_water/rho_shaft*M_here*dinamico_Adjoint(ginv(g_here))*G-DL_here*norm(eta_here(4:6))*eta_here;%-rho_water/rho_shaft

i_sig    = i_sig+1;
i_sig_nj = i_sig_nj+1;

%hook
i_sig    = i_sig+1;%joint

rho_hook  = S.VLinks(3).Rho;
M_hook    = S.VLinks(3).Ms{1};
nGauss    = S.VLinks(3).nGauss{1};
Fexth     = zeros(6*nGauss,1);

for i=1:nGauss
    M_here                    = M_hook((i-1)*6+1:i*6,:);
    g_here                    = g((i_sig-1)*4+1:i_sig*4,:);
    eta_here                  = eta((i_sig-1)*6+1:i_sig*6);
    DL_here                   = DL((i_sig_nj-1)*6+1:i_sig_nj*6,:);
    Fexth((i-1)*6+1:i*6)      = -rho_water/rho_hook*M_here*dinamico_Adjoint(ginv(g_here))*G-DL_here*norm(eta_here(4:6))*eta_here;%-rho_water/rho_hook
    
    i_sig    = i_sig+1;
    i_sig_nj = i_sig_nj+1;
end

%filament
i_sig    = i_sig+1;%joint

rho_filament  = S.VLinks(4).Rho;
M_filament    = S.VLinks(4).Ms{1};
nGauss        = S.VLinks(4).nGauss{1};
Fextf         = zeros(6*nGauss,1);

for i=1:nGauss
    M_here                    = M_filament((i-1)*6+1:i*6,:);
    g_here                    = g((i_sig-1)*4+1:i_sig*4,:);
    eta_here                  = eta((i_sig-1)*6+1:i_sig*6);
    DL_here                   = DL((i_sig_nj-1)*6+1:i_sig_nj*6,:);
    Fextf((i-1)*6+1:i*6)      = -rho_water/rho_filament*M_here*dinamico_Adjoint(ginv(g_here))*G-DL_here*norm(eta_here(4:6))*eta_here;%-rho_water/rho_filament
    
    i_sig    = i_sig+1;
    i_sig_nj = i_sig_nj+1;
end

Fext=[Fextb;Fexts;Fexth;Fextf];
end

