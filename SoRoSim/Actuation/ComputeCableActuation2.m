%Function to compute cable actuation when cable is not fully inside (24.05.2021)
function Bq = ComputeCableActuation2(S,dc,Sdiv,Ediv,J,g)

Bq        = zeros(S.ndof,1);
f         = 1;
i_sig     = 1;
dof_start = 1;

for i=1:S.N
    
    dof_here = S.Vtwists(f).dof;
    f        = f+1;
    i_sig    = i_sig+1;
    
    if S.VLinks(i).linktype=='r'
        i_sig = i_sig+1;
    end
    dof_start = dof_start+dof_here;
    
    dci       = dc{i};
    Sdivi     = Sdiv{i};
    Edivi     = Ediv{i};
    
    for j=1:S.VLinks(i).npie-1

            
        dof_here = S.Vtwists(f).dof;
        nGauss   = S.VLinks(i).nGauss{j};
        
        if j>=Sdivi&&j<=Edivi
            dcd      = dci{j};
            h_fn     = S.VLinks(i).h{j};
            w_fn     = S.VLinks(i).w{j};
            g0       = S.g0{f};
            lp       = S.VLinks(i).lp{j};
            hmax     = g0(2,4)+h_fn(0.5)/2; %check at mid point
            wmax     = g0(3,4)+w_fn(0.5)/2;
            Lscale   = lp;
            
            if hmax<dcd(2,1)||wmax<dcd(3,1)

                %scaling of quantities
                dcd       = dcd/Lscale;
                g_leftinv = [eye(3) -dcd(:,1);0 0 0 1];
                

                dc_right     = dcd(:,nGauss);

                g_right      = [eye(3) dc_right;0 0 0 1];
                g_iL         = ginv(g((i_sig-1)*4+1:i_sig*4,:))*g((i_sig+nGauss-1-1)*4+1:(i_sig+nGauss-1)*4,:);
                g_u          = g_leftinv*g_iL*g_right;
                u            = g_u(1:3,4);
                u            = u/sqrt(u'*u);
                
                J_iL                               = J((i_sig+nGauss-1-1)*6+1:(i_sig+nGauss-1)*6,:);
                Bq(dof_start:dof_start+dof_here-1) = J_iL(:,dof_start:dof_start+dof_here-1)'*dinamico_coAdjoint(g_right)*[0 0 0 u']'*Lscale; %scaling back 
            end
        end
        i_sig     = i_sig+nGauss;
        f         = f+1;
        dof_start = dof_start+dof_here;
    end
end

end

