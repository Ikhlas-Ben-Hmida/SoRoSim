%Function that calculates the generalized actuation matrix (Bq) at every
%significant point (24.05.2021)

function Bq = ActuationMatrix(S,q)

if isrow(q)
    q=q';
end

if S.Actuated
    J              = S.Jacobian(q);
    nact           = S.nact;
    ndof           = S.ndof;
    Bq             = zeros(ndof,nact);
    n_jact         = S.n_jact;
    
    %revolute, prismatic, helical joints
    Bqj1           = S.Bqj1;
    n_1dof         = size(Bqj1,2);
    Bq(:,1:n_1dof) = S.Bqj1;
    
    %for other joints
    i_jact         = S.i_jact;
    i_tau          = n_1dof+1;
    i_jactq        = S.i_jactq;
    
    for i = i_jact(i_tau:end)
        if S.VLinks(i).jointtype == 'U'
            Bq(i_jactq(i_tau:i_tau+1),i_tau:i_tau+1) = [1 0;0 1];
            i_tau                                    = i_tau+2;
            
        elseif S.VLinks(i).jointtype=='C'
            Bq(i_jactq(i_tau:i_tau+1),i_tau:i_tau+1) = [1 0;0 1];
            i_tau                                    = i_tau+2;
            
        elseif S.VLinks(i).jointtype == 'A'
            i_sig                    = 1;
            f                        = 1;
            
            for ii = 1:i-1
                for jj = 1:S.VLinks(ii).npie-1
                    i_sig            = i_sig+1+S.VLinks(ii).nGauss{jj};
                end
                
                if S.VLinks(ii).linktype == 'r'
                    i_sig            = i_sig+1;
                end
                f                    = f+S.VLinks(ii).npie;
            end
            
            J_here                                   = J((i_sig-1)*6+1:i_sig*6,:);
            S_here                                   = J_here(:,i_jactq(i_tau:i_tau+2));
            B_here                                   = S.Vtwists(f).B;
            
            Bq(i_jactq(i_tau:i_tau+2),i_tau:i_tau+2) = S_here'*B_here;
            i_tau                                    = i_tau+3;
            
        elseif S.VLinks(i).jointtype == 'S'
            i_sig = 1;
            
            for ii = 1:i-1
                for jj = 1:S.VLinks(ii).npie-1
                    i_sig = i_sig+1+S.VLinks(ii).nGauss{jj};
                end
                
                if S.VLinks(ii).linktype == 'r'
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
            for ii = 1:i-1
                for jj = 1:S.VLinks(ii).npie-1
                    i_sig = i_sig+1+S.VLinks(ii).nGauss{jj};
                end
                
                if S.VLinks(ii).linktype == 'r'
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
    n_sact = S.n_sact;
    N      = S.N;
    for ii=1:n_sact
        
         dcii = cell(1,N); dcpii = cell(1,N); Sdivii = cell(1,N); Edivii = cell(1,N);
            
        for i=1:N
            dcii{i}   = S.dc{ii,i};
            dcpii{i}  = S.dcp{ii,i};
            Sdivii{i} = S.Sdiv{ii,i};
            Edivii{i} = S.Ediv{ii,i};
        end
        Insideii     = S.Inside{ii};

        if Insideii
            Bq(:,n_jact+ii) = ComputeCableActuation(S,dcii,dcpii,Sdivii,Edivii,q);
        else
            g = S.FwdKinematics(q);
            Bq(:,n_jact+ii) = ComputeCableActuation2(S,dcii,Sdivii,Edivii,J,g);
        end
        
    end
    %add more type of actuations if needed
    
else
    Bq = 0;
end
end

