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
    n_jact         = Tr.n_jact;
    
    %revolute, prismatic, helical joints
    Bqj1           = Tr.Bqj1;
    n_1dof         = size(Bqj1,2);
    Bq(:,1:n_1dof) = Tr.Bqj1;
    
    %for other joints
    i_jact         = Tr.i_jact;
    i_tau          = n_1dof+1;
    i_jactq        = Tr.i_jactq;
    
    for i = i_jact(i_tau:end)
        if Tr.VLinks(i).jointtype == 'U'
            Bq(i_jactq(i_tau:i_tau+1),i_tau:i_tau+1) = [1 0;0 1];
            i_tau                                    = i_tau+2;
            
        elseif Tr.VLinks(i).jointtype=='C'
            Bq(i_jactq(i_tau:i_tau+1),i_tau:i_tau+1) = [1 0;0 1];
            i_tau                                    = i_tau+2;
            
        elseif Tr.VLinks(i).jointtype == 'A'
            i_sig                    = 1;
            f                        = 1;
            
            for iso = 1:i-1
                for jj = 1:Tr.VLinks(iso).npie-1
                    i_sig            = i_sig+1+Tr.VLinks(iso).nGauss{jj};
                end
                
                if Tr.VLinks(iso).linktype == 'r'
                    i_sig            = i_sig+1;
                end
                f                    = f+Tr.VLinks(iso).npie;
            end
            
            J_here                                   = J((i_sig-1)*6+1:i_sig*6,:);
            S_here                                   = J_here(:,i_jactq(i_tau:i_tau+2));
            B_here                                   = Tr.Vtwists(f).B;
            
            Bq(i_jactq(i_tau:i_tau+2),i_tau:i_tau+2) = S_here'*B_here;
            i_tau                                    = i_tau+3;
            
        elseif Tr.VLinks(i).jointtype == 'S'
            i_sig = 1;
            
            for iso = 1:i-1
                for jj = 1:Tr.VLinks(iso).npie-1
                    i_sig = i_sig+1+Tr.VLinks(iso).nGauss{jj};
                end
                
                if Tr.VLinks(iso).linktype == 'r'
                    i_sig = i_sig+1;
                end
            end
            
            J_here                                   = J((i_sig-1)*6+1:i_sig*6,:);
            S_here                                   = J_here(:,i_jactq(i_tau:i_tau+2));
            B_here                                   = [eye(3);zeros(3,3)];
            
            Bq(i_jactq(i_tau:i_tau+2),i_tau:i_tau+2) = S_here'*B_here;
            i_tau                                    = i_tau+3;
            
        else %free joint
            i_sig = 1;
            for iso = 1:i-1
                for jj = 1:Tr.VLinks(iso).npie-1
                    i_sig = i_sig+1+Tr.VLinks(iso).nGauss{jj};
                end
                
                if Tr.VLinks(iso).linktype == 'r'
                    i_sig = i_sig+1;
                end
            end
            J_here                                   = J((i_sig-1)*6+1:i_sig*6,:);
            S_here                                   = J_here(:,i_jactq(i_tau:i_tau+5));
            B_here                                   = eye(6);
            
            Bq(i_jactq(i_tau:i_tau+5),i_tau:i_tau+5) = S_here'*B_here;
            i_tau                                    = i_tau+6;
        end
    end
    
    %cable actuation
    n_sact = Tr.n_sact;
    N      = Tr.N;
    for iso=1:n_sact
        
        dcii = cell(1,N); dcpii = cell(1,N); Sdivii = cell(1,N); Edivii = cell(1,N);    
        for i=1:N
            dcii{i}   = Tr.dc{iso,i};
            dcpii{i}  = Tr.dcp{iso,i};
            Sdivii{i} = Tr.Sdiv{iso,i};
            Edivii{i} = Tr.Ediv{iso,i};
        end
        Insideii     = Tr.Inside{iso};

        if Insideii
            
            dof_start = 1;

            for i=1:Tr.N

                dof_start = dof_start+Tr.CVTwists{i}(1).dof;

                for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie-1

                    dof_here  = Tr.CVTwists{i}(j+1).dof;
                    if isempty(dcii{i}{j})
                        dof_start = dof_start+dof_here;
                        continue;
                    end

                    if j>=Sdivii{i}&&j<=Edivii{i}

                        Bdof       = Tr.CVTwists{i}(j+1).Bdof;
                        Bodr       = Tr.CVTwists{i}(j+1).Bodr;

                        dcd  = dcii{i}{j};
                        dcpd = dcpii{i}{j};

                        q_here    = q(dof_start:dof_start+dof_here-1);
                        xi_star   = Tr.CVTwists{i}(j+1).xi_star;
                        ld        = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
                        Ws        = Tr.VLinks(Tr.LinkIndex(i)).Ws{j};
                        nGauss    = Tr.VLinks(Tr.LinkIndex(i)).nGauss{j};
                        
                        q_scale_here  = Tr.q_scale(dof_start:dof_start+dof_here-1);
                        doftheta_here = Bdof(1:3)'*(Bodr(1:3)+1);
                        q_scale_here(1:doftheta_here) = q_scale_here(1:doftheta_here)*lpf;
                        B_scale = repmat(q_scale_here',6*nGauss,1);
                        B = Tr.CVTwists{i}(j+1).B./B_scale; %actual B


                        for ii=2:nGauss-1

                            dc_here     = dcd(:,ii);
                            dcp_here    = dcpd(:,ii);
                            B_here      = B(6*(ii-1)+1:6*ii,:);
                            xi_starhere = xi_star(6*(ii-1)+1:6*ii,1);
                            xi_here     = B_here*q_here+xi_starhere;
                            xihat_here  = dinamico_hat(xi_here);
                            tang        = xihat_here*[dc_here;1]+[dcp_here;0];
                            tang        = tang(1:3)/norm(tang(1:3)); %check with Federico
                            Btau        = [dinamico_tilde(dc_here)*tang;tang];
                            Bq(dof_start:dof_start+dof_here-1,Tr.n_jact+iso)  = Bq(dof_start:dof_start+dof_here-1,Tr.n_jact++iso)+ld*Ws(ii)*B_here'*Btau; 
                            
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

                    if isempty(dcii{i}{j})
                        dof_start = dof_start+dof_here;
                        continue;
                    end


                    nGauss   = Tr.VLinks(Tr.LinkIndex(i)).nGauss{j};

                    if j>=Sdivii{i}&&j<=Edivii{i}
                        dcd      = dcii{i}{j};
                        h_fn     = Tr.VLinks(Tr.LinkIndex(i)).h{j};
                        w_fn     = Tr.VLinks(Tr.LinkIndex(i)).w{j};
                        gi       = Tr.VLinks(Tr.LinkIndex(i)).gi{j};
                        hmax     = gi(2,4)+h_fn(0.5)/2; %check at mid point
                        wmax     = gi(3,4)+w_fn(0.5)/2;

                        if hmax<dcd(2,1)||wmax<dcd(3,1)

                            %scaling of quantities
                            g_leftinv    = [eye(3) -dcd(:,1);0 0 0 1];
                            dc_right     = dcd(:,nGauss);
                            g_right      = [eye(3) dc_right;0 0 0 1];
                            g_iL         = ginv(g((i_sig-1)*4+1:i_sig*4,:))*g((i_sig+nGauss-1-1)*4+1:(i_sig+nGauss-1)*4,:);
                            g_u          = g_leftinv*g_iL*g_right;
                            u            = g_u(1:3,4);
                            u            = u/sqrt(u'*u);

                            J_iL                               = J((i_sig+nGauss-1-1)*6+1:(i_sig+nGauss-1)*6,:);
                            Bq(dof_start:dof_start+dof_here-1) = J_iL(:,dof_start:dof_start+dof_here-1)'*dinamico_coAdjoint(g_right)*[0 0 0 u']'; %scaling back 
                        end
                    end
                    i_sig     = i_sig+nGauss;
                    dof_start = dof_start+dof_here;
                end
            end
        end
        
    end
    %add more type of actuations if needed
    
else
    Bq = 0;
end
end

