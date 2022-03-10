%Function to compute cable actuation when cable is not fully inside (24.05.2021)
function Bq = ComputeCableActuation2(Tr,dc,Sdiv,Ediv,J,g)

Bq        = zeros(Tr.ndof,1);
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
        
        if isempty(dc{i}{j})
            dof_start = dof_start+dof_here;
            continue;
        end
            
        
        nGauss   = Tr.VLinks(Tr.LinkIndex(i)).nGauss{j};
        
        if j>=Sdiv{i}&&j<=Ediv{i}
            dcd      = dc{i}{j};
            h_fn     = Tr.VLinks(Tr.LinkIndex(i)).h{j};
            w_fn     = Tr.VLinks(Tr.LinkIndex(i)).w{j};
            gi       = Tr.VLinks(Tr.LinkIndex(i)).gi{j};
            ld       = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
            hmax     = gi(2,4)+h_fn(0.5)/2; %check at mid point
            wmax     = gi(3,4)+w_fn(0.5)/2;
            Lscale   = ld;
            
            if hmax<dcd(2,1)||wmax<dcd(3,1)

                %scaling of quantities
                dcd       = dcd/Lscale;
                g_leftinv = [eye(3) -dcd(:,1);0 0 0 1];
                dc_right  = dcd(:,nGauss);

                g_right   = [eye(3) dc_right;0 0 0 1];
                g_iL      = ginv(g((i_sig-1)*4+1:i_sig*4,:))*g((i_sig+nGauss-1-1)*4+1:(i_sig+nGauss-1)*4,:);
                g_u       = g_leftinv*g_iL*g_right;
                u         = g_u(1:3,4);
                u         = u/sqrt(u'*u);
                
                J_iL                               = J((i_sig+nGauss-1-1)*6+1:(i_sig+nGauss-1)*6,:);
                Bq(dof_start:dof_start+dof_here-1) = J_iL(:,dof_start:dof_start+dof_here-1)'*dinamico_coAdjoint(g_right)*[0 0 0 u']'*Lscale; %scaling back 
            end
        end
        i_sig     = i_sig+nGauss;
        dof_start = dof_start+dof_here;
    end
end

end

