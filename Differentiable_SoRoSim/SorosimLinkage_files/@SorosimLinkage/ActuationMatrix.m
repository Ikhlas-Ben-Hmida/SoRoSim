%Function that calculates the generalized actuation matrix (Bq)
%Inspired by dynamicsSolver and updated to match SorosimLinkage properties
function Bq = ActuationMatrix(Linkage,q)
    if isrow(q)
        q = q';
    end
    if Linkage.Actuated
        
        J    = Linkage.Jacobian(q);
        nact = Linkage.nact;
        ndof = Linkage.ndof;
        Bq   = zeros(ndof,nact);
        
        N_jact = Linkage.N_jact; % number of Links whose joints are actuated
        n_jact = Linkage.n_jact; % total number of joint actuators
        
        %% 1. Rigid/Joint Actuation
        % revolute, prismatic, helical joints (1D joints)
        Bj1   = Linkage.Bj1; 
        na_j1 = size(Bj1,2); 
        if na_j1 > 0
            Bq(:,1:na_j1) = Bj1;
        end
        
        % for other multi-DOF joints
        i_jact  = Linkage.i_jact;
        i_u     = na_j1 + 1;
        i_jactq = Linkage.i_jactq;
        
        for ii = na_j1+1:N_jact
            i = i_jact(ii);
            dof_here = Linkage.CVRods{i}(1).dof;
            us_here  = i_u : i_u+dof_here-1;
            dofs_here = i_jactq(us_here);
            
            if Linkage.VLinks(Linkage.LinkIndex(i)).jointtype == 'C'
                Bq(dofs_here, us_here) = eye(2);
            else 
                % Utilizing ActuationPrecompute as done in dynamicsSolver
                i_sig = Linkage.ActuationPrecompute.i_sig_act(ii); 
                
                if Linkage.VLinks(Linkage.LinkIndex(i)).linktype == 'r'
                    gi = Linkage.VLinks(Linkage.LinkIndex(i)).gi;
                else
                    gi = Linkage.VLinks(Linkage.LinkIndex(i)).gi{1};
                end
                
                AdgstepinvS_here = dinamico_Adjoint(gi) * J((i_sig-1)*6+1:i_sig*6, i_jactq(i_u:i_u+dof_here-1));
                Phi_here = Linkage.CVRods{i}(1).Phi;
                Bq(dofs_here, us_here) = AdgstepinvS_here' * Phi_here;
            end
            i_u = i_u + dof_here;
        end
        
        %% 2. Soft/Cable Actuation
        if Linkage.n_sact > 0
            dof_start = 1;
            
            % Setup friction tracking if requested via Linkage properties
            use_friction = isprop(Linkage, 'use_friction') && Linkage.use_friction;
            if use_friction
                mu_array = Linkage.mu;
                
                % Format mu_array to strictly map to the soft actuators (1 to n_sact)
                if isscalar(mu_array)
                    mu_array = mu_array * ones(Linkage.n_sact, 1);
                elseif length(mu_array) == nact
                    % If passed as nact, extract only the soft/cable actuators
                    mu_array = mu_array(n_jact+1 : end);
                elseif length(mu_array) ~= Linkage.n_sact
                    error('ActuationMatrix: Dimension mismatch. Linkage.mu must be a scalar, length nact, or length n_sact.');
                end
                
                gs = Linkage.FwdKinematics(q); % Fetch global poses
                
                phi_sum = zeros(Linkage.n_sact, 1);
                f_prev  = cell(Linkage.n_sact, 1);
                r_prev  = cell(Linkage.n_sact, 1);
                R_prev  = cell(Linkage.n_sact, 1);
                d_prev  = cell(Linkage.n_sact, 1);
            end
            
            i_sig = 1; % Global integration point tracker for gs matrix
            
            for i = 1:Linkage.N
                
                dof_here = Linkage.CVRods{i}(1).dof;
                if ~Linkage.OneBasis
                    dof_start = dof_start + dof_here;
                end
                
                i_sig = i_sig + 1; % Skip Joint DoF frame in gs
                
                for j = 1:Linkage.VLinks(Linkage.LinkIndex(i)).npie-1
                    na_here   = length(Linkage.i_sact{i}{j});
                    dof_here  = Linkage.CVRods{i}(j+1).dof;
                    dofs_here = dof_start : dof_start+dof_here-1;
                    
                    if ~Linkage.OneBasis
                        dof_start = dof_start + dof_here;
                    end
                    
                    Ws  = Linkage.CVRods{i}(j+1).Ws;
                    nip = Linkage.CVRods{i}(j+1).nip; % Defining nip locally!
                    B_here = zeros(dof_here, na_here);
                    
                    % We MUST iterate through all nip points to keep i_sig synced 
                    % with the FwdKinematics gs matrix, even if na_here == 0.
                    for ii = 1:nip
                        if Ws(ii) > 0
                            
                            % 1. Extract current global pose if friction is used
                            if use_friction
                                current_g = gs((i_sig-1)*4 + 1 : i_sig*4, :);
                                r_curr = current_g(1:3, 4);
                                R_curr = current_g(1:3, 1:3);
                            end
                            
                            % 2. Only compute mechanics if there are actuators
                            if na_here > 0
                                Phi = Linkage.CVRods{i}(j+1).Phi((ii-1)*6+1:ii*6,:);
                                xi  = Phi*q(dofs_here) + Linkage.CVRods{i}(j+1).xi_star((ii-1)*6+1:ii*6,1);
                                
                                Phi_a = zeros(6,na_here);
                                xihat_123 = [0 -xi(3) xi(2) xi(4); xi(3) 0 -xi(1) xi(5); -xi(2) xi(1) 0 xi(6)];
                                
                                ia_here = 1;
                                for ia = Linkage.i_sact{i}{j}
                                    dc  = Linkage.dc{ia,i}{j}(:,ii);
                                    dcp = Linkage.dcp{ia,i}{j}(:,ii);
                                    
                                    % Base spatial actuation matrix column
                                    Phi_a_col = SoftActuator_FD(dc, dcp, xihat_123);
                                    
                                    % Friction Calculation (Model 2: Configuration-Dependent)
                                    friction_weight = 1.0;
                                    
                                    % Fetch actuator-specific friction coefficient
                                    if use_friction
                                        mu_ia = mu_array(ia);
                                        
                                        if mu_ia > 0
                                            if isempty(r_prev{ia})
                                                % First point initialization
                                                r_prev{ia} = r_curr;
                                                R_prev{ia} = R_curr;
                                                d_prev{ia} = dc;
                                            else
                                                % Compute cable force direction vector (p)
                                                p_curr = (r_curr - r_prev{ia}) + (R_curr * dc) - (R_prev{ia} * d_prev{ia});
                                                p_norm = norm(p_curr);
                                                
                                                if p_norm > 1e-12
                                                    % Normalize to get direction (f)
                                                    f_curr = p_curr / p_norm;
                                                    
                                                    if ~isempty(f_prev{ia})
                                                        % Compute angle between consecutive directions
                                                        cos_phi = max(min(f_prev{ia}' * f_curr, 1), -1); % Clamp to prevent acos(NaN)
                                                        phi_l = acos(cos_phi);
                                                        phi_sum(ia) = phi_sum(ia) + phi_l;
                                                    end
                                                    f_prev{ia} = f_curr;
                                                end
                                                
                                                % Update trackers for the next Gauss point
                                                r_prev{ia} = r_curr;
                                                R_prev{ia} = R_curr;
                                                d_prev{ia} = dc;
                                            end
                                            
                                            % Apply exponential decay based on cumulative angle
                                            friction_weight = exp(-mu_ia * phi_sum(ia));
                                        end
                                    end
                                    
                                    % Apply friction weight to spatial actuation column
                                    Phi_a(:,ia_here) = Phi_a_col * friction_weight;
                                    ia_here = ia_here + 1;
                                end
                                B_here = B_here + Ws(ii) * Phi' * Phi_a;
                            end
                            
                            % 3. Increment global integration point tracker
                            i_sig = i_sig + 1; 
                        end
                    end
                    
                    % Map B_here into the full Bq matrix using n_jact offset
                    if na_here > 0
                        Bq(dofs_here, n_jact + Linkage.i_sact{i}{j}) = B_here;
                    end
                end
            end
        end
        
    else
        Bq = 0;
    end
end