%Function to calculate the custom external point force. Make sure to edit
%CustomExtForce.m file too
%Last modified by Anup Teejo Mathew 30/06/2021

function Fext=CustomExtPointForce(Tr,q,g,J,t,qd,eta,Jdot)

%%%%NOTE%%%%
%Tr: Linkage element,
%q and qd: joint coordinates and their time derivatives,
%g, J, Jd, and eta: transformation matrix, Jacobian, time derivative of jacobian, and screw velocity at every significant point of the linkage
%t: time
%Currently applicable only for independent bases

%Fext should be 6*n column vector where n is the total number of integration points of all soft links (nip) + number of rigid links.
%(Example: linkage with 2 soft links and 1 rigid link (n=nip1+nip2+1)
%Fext should be arranged according to the order of precedence
%Fext should be point force

% Significant points: 1 for every joint, 1 at the center of the rigid link, for soft links at every integration points

% J   = S.Jacobian(q);         %geometric jacobian of the linkage calculated at every significant points
% g   = S.FwdKinematics(q);    %transformation matrix of the linkage calculated at every significant points
% eta = S.ScrewVelocity(q,qd); %Screwvelocity of the linkage calculated at every significant points
% J   = S.Jacobiandot(q,qd);   %time derivative of geometric jacobian of the linkage calculated at every significant points

%%%END%%%

n    = Tr.nsig-Tr.N; %everything except rigid joints
Fext = zeros(6*n,1);

%% Contact Force
if t>14
rr = Tr.VLinks(Tr.LinkIndex(Tr.N)).r(0);

k=50;
b=5;

%rod position
Rp  = Tr.g_ini(end-3:end-1,4);

i_sig_nj = 1;
i_sig = 1;

for i=1:12 %link 12 is the last top flagella
    i_sig    = i_sig+1;%joint
    if Tr.VLinks(Tr.LinkIndex(i)).linktype=='r'
        if i==1
            g_R = g((i_sig-1)*4+1:i_sig*4,:);
            eta_R = eta((i_sig-1)*6+1:i_sig*6);
            PosHS = g_R*Tr.CP3;
            nbp   = size(PosHS,2);
            
            for ii=1:nbp
                
                d_rp    = PosHS(1:3,ii)-Rp; %rod to point
                d_rp(2) = 0; %rod is along y axis
                d_rp_n =norm(d_rp);
                
                d21 = rr-d_rp_n;

                if d21>0 %contact
                    u_rp   = d_rp/d_rp_n; % unit vector from rod to point
                    fc = k*d21*u_rp;
               
                    gi  = [eye(3),PosHS(1:3,ii);0 0 0 1];
                    gRi = ginv(g_R)*gi;
                    eta_i = dinamico_Adjoint(ginv(gRi))*eta_R;
                    d21dot = u_rp'*eta_i(4:6);
                    
                    if d21dot<0 %velocity opposite direction of unit vector (into the rod)
                        fc = fc-b*d21dot*u_rp;
                    end
                   
                    Fext((i_sig_nj-1)*6+1:i_sig_nj*6) = Fext((i_sig_nj-1)*6+1:i_sig_nj*6)+dinamico_coAdjoint(gRi)*[0;0;0;fc];
                end
                
            end
        end
        
        i_sig    = i_sig+1;
        i_sig_nj = i_sig_nj+1;
    end
    for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie-1
        
        nip = Tr.CVTwists{i}(j+1).nip;
        if i==11||i==12
            Xs  = Tr.CVTwists{i}(j+1).Xs;
            r   = Tr.VLinks(Tr.LinkIndex(i)).r{j};

            for ii=1:nip
                g_here = g((i_sig-1)*4+1:i_sig*4,:);
                d_rp    = g_here(1:3,4)-Rp; %rod to point
                d_rp(2) = 0; %y coordinate is irrelevant here
                d_rp_n  = norm(d_rp);
                r_here = r(Xs(ii));
                d21    = rr+r_here-d_rp_n;

                if d21>0 %if contact
                    u_rp = d_rp/d_rp_n;%vector in direction from rod to center of cross section
                    fc = k*d21*u_rp;
                    
                    eta_here = eta((i_sig-1)*6+1:i_sig*6);
                    g_here(1:3,4) = zeros(3,1);
                    eta_d  = dinamico_Adjoint(g_here)*eta_here; % rotated into global frame
                    d21dot = eta_d(4:6)'*u_rp; %normal component add force only if negative
                    
                    if d21dot<0 %velocity opposite direction of unit vector (into the rod)
                        fc = fc-b*d21dot*u_rp;
                    end

                    Ad_g_here_inv = dinamico_Adjoint(ginv(g_here));
                    Fext((i_sig_nj-1)*6+1:i_sig_nj*6) = Fext((i_sig_nj-1)*6+1:i_sig_nj*6)+Ad_g_here_inv*[0;0;0;fc]; %global to local, just rotation
                end
                i_sig    = i_sig+1;
                i_sig_nj = i_sig_nj+1;
            end
        else
            i_sig    = i_sig+nip;
            i_sig_nj = i_sig_nj+nip;
        end
    end
end
end
end