%Function that calculates the generalized actuation matrix (Bq) at every
%Last modified by Anup Teejo Mathew 02.03.2022

function Bq = ActuationMatrix(Tr,q)

if isrow(q)
    q=q';
end

if Tr.Actuated
    
    J              = Tr.Jacobian(q);
    nact           = Tr.nact;
    ndof           = Tr.ndof;
    Bq             = zeros(ndof,nact);
    
    %revolute, prismatic, helical joints
    Bqj1           = Tr.Bqj1;
    n_1dof         = size(Bqj1,2);
    Bq(:,1:n_1dof) = Tr.Bqj1;
    
    %for other joints
    i_jact  = Tr.i_jact;
    i_u     = n_1dof+1;
    i_jactq = Tr.i_jactq;
    
    n_ljact = length(i_jact);
    for iii=n_1dof+1:n_ljact
        i = i_jact(iii);
        i_sig = 1;
        dof_start = 1;
        for ii=1:i-1
            i_sig = i_sig+1;
            dof_start = dof_start+Tr.CVTwists{ii}(1).dof; %joint
            for j=1:Tr.VLinks(Tr.LinkIndex(ii)).npie-1
                i_sig = i_sig+Tr.CVTwists{ii}(j+1).nip;
                dof_start = dof_start+Tr.CVTwists{ii}(j+1).dof;
            end
            if Tr.VLinks(Tr.LinkIndex(ii)).linktype=='r'
                i_sig = i_sig+1;
            end
        end

        if Tr.VLinks(Tr.LinkIndex(i)).jointtype=='C'
            dof_here = dof_start+Tr.CVTwists{i}(1).dof;
            Bq(i_jactq(i_u:i_u+dof_here-1),i_u:i_u+dof_here-1) = [1 0;0 1];
            i_u = i_u+2;
        else 
            dof_here = Tr.CVTwists{i}(1).dof;
            J_here = J((i_sig-1)*6+1:i_sig*6,:);
            S_here = J_here(:,i_jactq(i_u:i_u+dof_here-1));
            B_here = Tr.CVTwists{i}(1).B;

            Bq(i_jactq(i_u:i_u+dof_here-1),i_u:i_u+dof_here-1) = S_here'*B_here;
            i_u = i_u+dof_here;
        end
    end
    
    %cable actuation
    n_sact = Tr.n_sact;
    N      = Tr.N;
    for iso=1:n_sact
        
        dcii = cell(1,N); dcpii = cell(1,N); Sdivii = zeros(N,1); Edivii = zeros(N,1);  
        for i=1:N
            dcii{i}   = Tr.dc{iso,i};
            dcpii{i}  = Tr.dcp{iso,i};
            Sdivii(i) = Tr.Sdiv(iso,i);
            Edivii(i) = Tr.Ediv(iso,i);
        end
        Insideii      = Tr.Inside(iso);

        if Insideii
            
            dof_start = 1;

            for i=1:Tr.N

                dof_start = dof_start+Tr.CVTwists{i}(1).dof;

                for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie-1

                    dof_here  = Tr.CVTwists{i}(j+1).dof;
                    VTwists  = Tr.CVTwists{i};

                    if j>=Sdivii(i)&&j<=Edivii(i)

                        dcd  = dcii{i}{j};
                        dcpd = dcpii{i}{j};

                        q_here  = q(dof_start:dof_start+dof_here-1);
                        xi_star = VTwists(j+1).xi_star;
                        ld      = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
                        Ws      = VTwists(j+1).Ws;
                        nip     = VTwists(j+1).nip;

                        %scaling of quantities
                        Lscale       = ld;

                        for k=1:nip
                            if Ws(k)>0
                                dc_here     = dcd(:,k);
                                dcp_here    = dcpd(:,k);
                                xi_here = xi_star(6*(k-1)+1:6*k,1);
                                xi_here(1:3) = xi_here(1:3)*Lscale; %scaling using the formula: Lscale m = 1 unit

                                dBqdq_here = VTwists(j+1).B((k-1)*6+1:k*6,:);
                                if dof_here>0
                                    xi_here    = dBqdq_here*q_here+xi_here;
                                end
                                xi_here(1:3) = xi_here(1:3)/Lscale; %unscaling
                                xihat_here  = [0 -xi_here(3) xi_here(2) xi_here(4);xi_here(3) 0 -xi_here(1) xi_here(5);-xi_here(2) xi_here(1) 0 xi_here(6)];%4th row is avoided to speedup calculation
                                tang        = xihat_here*[dc_here;1]+dcp_here;
                                tang        = tang/norm(tang);
                                Btau        = [[0 -dc_here(3) dc_here(2);dc_here(3) 0 -dc_here(1);-dc_here(2) dc_here(1) 0]*tang;tang];
                                Btau(1:3)   = Btau(1:3)/Lscale;

                                Bq(dof_start:dof_start+dof_here-1,Tr.n_jact+iso)  = Bq(dof_start:dof_start+dof_here-1,Tr.n_jact+iso)+Ws(k)*dBqdq_here'*Btau*Lscale; %scaled back for addition
                            end
                        end

                    end
                    dof_start = dof_start+dof_here;
                end 
            end
        else
            g = Tr.FwdKinematics(q);
            i_sig     = 1;
            dof_start = 1;
            for i=1:Tr.N
    
                dof_start = dof_start+Tr.CVTwists{i}(1).dof;
                i_sig    = i_sig+1;
                if Tr.VLinks(Tr.LinkIndex(i)).linktype=='r'
                    i_sig = i_sig+1;
                end


                for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie-1

                    dof_here = Tr.CVTwists{i}(j+1).dof; 
                    nip = Tr.CVTwists{i}(j+1).nip;

                    if j==Sdiv(i)
                        dof_act_start = dof_start;
                        dof_act = 0;
                        i_sig_l = i_sig;
                        dc_l = dcii{i}{j}(:,1);
                    end
                    if j==Ediv(i)
                        i_sig_r = i_sig+nip-1;
                        dc_r = dcii{i}{j}(:,end);
                    end
                    if j>=Sdiv(i)&&j<=Ediv(i)
                        dof_act = dof_act+dof_here;
                    end

                    i_sig = i_sig+nip;
                    dof_start = dof_start+dof_here;
                end

                ld       = Tr.VLinks(Tr.LinkIndex(i)).lp{Ediv(i)};
                Lscale   = ld;

                %scaling of quantities
                dc_l       = dc_l/Lscale;
                dc_r       = dc_r/Lscale;

                g_leftinv = [eye(3) -dc_l;0 0 0 1];
                g_right   = [eye(3) dc_r;0 0 0 1];

                g_iL      = ginv(g((i_sig_l-1)*4+1:i_sig_l*4,:))*g((i_sig_r-1)*4+1:(i_sig_r)*4,:);
                g_u       = g_leftinv*g_iL*g_right;
                u         = g_u(1:3,4);
                u         = u/sqrt(u'*u);

                J_iL = J((i_sig_r-1)*6+1:(i_sig_r)*6,:);
                Bq(dof_act_start:dof_act_start+dof_act-1,Tr.n_jact+iso) = J_iL(:,dof_act_start:dof_act_start+dof_act-1)'*dinamico_coAdjoint(g_right)*[0 0 0 u']'*Lscale; %scaling back 

            end
        end
        
    end
    %add more type of actuations if needed
    
else
    Bq = 0;
end
end

