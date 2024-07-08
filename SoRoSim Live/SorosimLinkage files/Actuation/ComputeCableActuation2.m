%Function to compute cable actuation when cable is not fully inside
%only for rectangular cross sections (last modified 04/18/2023, Anup)
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
        nip = Tr.CVTwists{i}(j+1).nip;
         
        if j==Sdiv(i)
            dof_act_start = dof_start;
            dof_act = 0;
            i_sig_l = i_sig;
            dc_l = dc{i}{j}(:,1);
        end
        if j==Ediv(i)
            i_sig_r = i_sig+nip-1;
            dc_r = dc{i}{j}(:,end);
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
    Bq(dof_act_start:dof_act_start+dof_act-1) = J_iL(:,dof_act_start:dof_act_start+dof_act-1)'*dinamico_coAdjoint(g_right)*[0 0 0 u']'*Lscale; %scaling back 

end

end

