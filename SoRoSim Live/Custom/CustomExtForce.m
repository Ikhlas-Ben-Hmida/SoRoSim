%Function to calculate the custom external force Make sure to edit
%CustomExtPointForce.m file too
%Last modified by Anup Teejo Mathew 30/06/2021

function Fext=CustomExtForce(Tr,q,g,J,t,qd,eta,Jdot)

%%%%NOTE%%%%
%Tr: Linkage element,
%q and qd: joint coordinates and their time derivatives,
%g, J, Jd, and eta: transformation matrix, Jacobian, time derivative of jacobian, and screw velocity at every significant point of the linkage
%t: time
%Currently applicable only for independent bases

%Fext should be 6*n column vector where n is the total number of integration points of all soft links (nip) + number of rigid links.
%(Example: linkage with 2 soft links and 1 rigid link (n=nip1+nip2+1)
%Fext should be arranged according to the order of precedence
%Fext should be distributed force for a soft link and point force for a
%rigid link.

% Significant points: 1 for every joint, 1 at the center of the rigid link, for soft links at every integration points

% J   = S.Jacobian(q);         %geometric jacobian of the linkage calculated at every significant points
% g   = S.FwdKinematics(q);    %transformation matrix of the linkage calculated at every significant points
% eta = S.ScrewVelocity(q,qd); %Screwvelocity of the linkage calculated at every significant points
% J   = S.Jacobiandot(q,qd);   %time derivative of geometric jacobian of the linkage calculated at every significant points

%%%END%%%

n    = Tr.nsig-Tr.N; %everything except rigid joints
Fext = zeros(6*n,1);

%make changes from here

%% Fluid Interaction

DL     = Tr.CP1;
G      = Tr.G;

i_sig    = 1;
i_sig_nj = 1;

Rho_water  = 1000;

for i=1:Tr.N 
    i_sig    = i_sig+1;%joint
    if Tr.VLinks(Tr.LinkIndex(i)).linktype=='r'
        Rho      = Tr.VLinks(Tr.LinkIndex(i)).Rho;
        M_here   = Tr.VLinks(Tr.LinkIndex(i)).M;
        g_here   = g((i_sig-1)*4+1:i_sig*4,:);
        eta_here = eta((i_sig-1)*6+1:i_sig*6);
        DL_here  = DL((i_sig_nj-1)*6+1:i_sig_nj*6,:);
        
        M_here(1:3,1:3)=zeros(3,3);
        
        if Tr.Gravity
            Fext((i_sig_nj-1)*6+1:i_sig_nj*6) = -Rho_water/Rho*M_here*dinamico_Adjoint(ginv(g_here))*G-DL_here*norm(eta_here(4:6))*eta_here;
        else
            Fext((i_sig_nj-1)*6+1:i_sig_nj*6) = -DL_here*norm(eta_here(4:6))*eta_here;
        end

        
        i_sig    = i_sig+1;
        i_sig_nj = i_sig_nj+1;
    end
    for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie-1
        Rho = Tr.VLinks(Tr.LinkIndex(i)).Rho;
        Ms  = Tr.CVTwists{i}(j+1).Ms;
        nip = Tr.CVTwists{i}(j+1).nip;
        
        for ii=1:nip
    
            M_here   = Ms((ii-1)*6+1:ii*6,:);
            g_here   = g((i_sig-1)*4+1:i_sig*4,:);
            eta_here = eta((i_sig-1)*6+1:i_sig*6);
            DL_here  = DL((i_sig_nj-1)*6+1:i_sig_nj*6,:);
            
            if i==8||i==9||i==12||i==13
                DL_here(5,6)=-DL_here(5,6);
                DL_here(6,5)=-DL_here(6,5);
            end
            
            if Tr.Gravity
                Fext((i_sig_nj-1)*6+1:i_sig_nj*6)  = -Rho_water/Rho*M_here*dinamico_Adjoint(ginv(g_here))*G-DL_here*norm(eta_here(4:6))*eta_here;%
            else
                Fext((i_sig_nj-1)*6+1:i_sig_nj*6)  = -DL_here*norm(eta_here(4:6))*eta_here;%
            end

            i_sig    = i_sig+1;
            i_sig_nj = i_sig_nj+1;
        end

    end
end

end