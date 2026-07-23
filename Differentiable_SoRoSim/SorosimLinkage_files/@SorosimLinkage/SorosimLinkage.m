%Creates the SorosimLinkage class by joining SorosimLinks
%Last modified by Anup Teejo Mathew 05.12.2024

classdef SorosimLinkage
    properties% Independent properties of Linkage class
        %%
        %General Properties

        N            %Total number of Links (not total number of rods)
        ndof         %Total number of DOF
        nsig         %Total number of points at which compulations are performed (significant points): N+N_rigid+sum(N_soft_div*nGauss_div)
        nj           %Total number of virtual joints (rigid joints plus virutal joints of all soft rods)
        VLinks       %Vector of all unique links (obtained from user input) VLinks(
        LinkIndex    %(Nx1) array of indices corresponding to each links. ith Link = S.VLinks(LinkIndex(i))
        CVRods       %Cell element of Rod vectors for each link S.CVRods{i}(j) i: Link index, j=1 for joint, j=2 to ndiv+1 for divisions.

        iLpre        %(Nx1) array corresponding to the Link index of the Link to which the ith Link is connected
        g_ini        %(4Nx4) Fixed initial transformation matrices of Links wrt to the tip of its previous link
        Z_order      %Order of Zannah collocation (2, 4, or 6) default value is 4
        OneBasis     %Enable if every rod is governed by the same q vector. Example a POD basis. (logical 0 by default)

        %Closed Loop Joints (Link A is connect to Link B via a closed loop joint)

        nCLj         %Total number of closed loop joints (default 0)
        iCLA         %(nCLjx1)array corresponding to the index of Link A
        iCLB         %(nCLjx1)array corresponding to the index of Link B
        VRodsCLj     %(nCLjx1)array of Rod vectors corresponding to each closed loop joint
        gACLj        %(nCLjx1)cells of fixed transformation from the tip of Link A to the close loop joint
        gBCLj        %(nCLjx1)cells of fixed transformation from the tip of Link B to the close loop joint
        CLprecompute %Struct element which contains pre-computed Phi_p (cell element of constrain basis), i_sigA (array of significant index corresponding to A), i_sigB (array ofsignificant index corresponding to B), and nCLp (total number of constraints)
        T_BS         %Baumgarte stabilization constant. Lower the value stricter the constrain.
        
        %External Force Properties
        Gravity       %logical 1 if gravity is present, 0 if not
        G             %Value of G

        %Point forces/moments
        np            %Number of point wrenches (default 0)
        LocalWrench   %logical 1 if point force/moment is a follower force (local frame) and 0 if it is wrt global frame
        Fp_loc        %Cell element of Link, Division numbers, and location if along a soft link: [i,j,X]. Corresponding to the point force/moment location (at the end of a soft division/center of mass of rigid link)
        Fp_vec        %Cell element with value of point forces/moments [Mx My Mz Fx Fy Fz]'
        Fp_sig        %Precomputed value of significant point corresponding to point wrenches

        %Custom external force
        CEF          %logical 1 if custom external force is present 0 if not (default value is 0)

        % Underwater Simulation
        UnderWater    %Enable if it is an underwater simulation (default logical 0)
        Rho_water     %Water density by default 1000 kg/m^3
        M_added       %Added mass for fluid simulations (by default zeros(6*nsig,6))
        DL            %Drag-Lift Matrix (by default zeros (6*nsig,6))

        %Actuation Properties
        Actuated      %logical 1 if actuated 0 if not
        nact          %Total number of actuators

        %Joint actuation parameters
        Bj1                 %Pre-computed actuation base for joints with dof = 1
        N_jact              %Total number of links whos joints are actuated
        n_jact              %Total number of joint actuators (eg. 3 for a spherical joint)
        i_jact              %Array of index of links whos joints are actuated
        i_jactq             %(n_jactx1) array of joint coordinate index of all active joints
        WrenchControlled    %(n_jactx1) array of logical 1 if the joint is wrench controlled 0 if joint coordinate controlled
        ActuationPrecompute %struct of precomputed actuation parameters: n_k (total number of constrained DoF), index_q_u, index_q_k, index_u_k, index_u_u

        %Thread-like actuator for soft links
        n_sact            %Number of soft link actuators
        i_sact            %Index of soft actuators present in i th link, jth divisions. Syntax i_sact{i}{j} will give you a row vector with soft actuator index
        dc                %(n_sactxN) cells of local cable position (0, yp, zp) at Gauss quadrature points of all active soft divisions
        dcp               %(n_sactxN) cells of space derivative of the local cable position (0,yp',zp')
        Sdiv              %(n_sactxN) cells of starting division number
        Ediv              %(n_sactxN) cells of ending division  number
        CableFunction     %Struct with cell elements (Cy_fn and Cz_fn) of parameterized functions corresponding to the y and z coodinates of the cable

        %Custom actuation
        CA               %logical 1 if custom actation is present 0 if not (default value: 0)
        CAI              %logical 1 to apply a custom actuator strength (default value: 0)

        %Pre-computed elastic Properties
        K       %Generalized Stiffness matrix
        Damped  %logical 1 if the soft links are elastically damped logical 0 if not (default value is 1)
        D       %Generalized Damping matrix

        CP1     %Custom constant properties of linkage that can be useful
        CP2
        CP3

        %Plotting
        PlotParameters

        %%%PlotParameters is a struct with following elements%%%
        %Lscale                   %Scaling factor for axis
        %CameraPosition           %CameraPosition with respect to origin
        %CameraUpVector           %Orientation of normal
        %CameraTarget             %Target location
        %Light                    %logical 1 if the light is on, 0 if not. (default value: 1)
        %Az_light                 %Light Azimuth wrt camera
        %El_light                 %Light Elevation wrt camera
        %XLim                     %x limit [X_lower_limt X_upper_limit]
        %YLim                     %y limit [Y_lower_limt Y_upper_limit]
        %ZLim                     %z limit [Z_lower_limt Z_upper_limit]
        %FrameRateValue           %FrameRate for dyanmic plot
        %ClosePrevious            %Logical 0 to not close previous image, 1 to close. (default value: 1)
        %CameraRotationSpeed      %For a cinematic rotation of the scene (default 0)
        %VideoResolution          %1 for full screen resolution (higher the better)

        % --- Friction Properties ---
        use_friction = false;
        mu = 0;

    end
    %%
    properties (Dependent = true) %Properties called from Link Class
        linktype
        jointtype
        CS
    end
    %%
    methods

        % Class Constructor Function
        function Linkage = SorosimLinkage(varargin)
            if strcmp(varargin{1},'empty')
                disp('Creating empty SorosimLinkage');

                Linkage.N = 0;
                Linkage.ndof = 0;
                Linkage.nsig = 0;
                Linkage.nCLj = 0;
                Linkage.OneBasis = false;
                Linkage.T_BS = 0.01;
                Linkage.Gravity = false;
                Linkage.G = 0;
                Linkage.np = 0;
                Linkage.UnderWater = false;
                Linkage.Rho_water = 1000;
                Linkage.Actuated = false;
                Linkage.nact = 0;
                Linkage.n_jact = 0;
                Linkage.WrenchControlled = 0;
                Linkage.D = 0;
                Linkage.Damped = true;
                Linkage.K = 0;
                
                Linkage.Z_order = 4; %the value should be 2 or 4
    
                Lscale = 1;
                PlotParameters.Lscale         = Lscale;
                PlotParameters.CameraPosition = [Lscale*2 -Lscale/2 Lscale/2];
                PlotParameters.CameraTarget   = [0 0 0];
                PlotParameters.CameraUpVector = [0 0 1];
                PlotParameters.Light          = true;
                PlotParameters.Az_light       = 0;
                PlotParameters.El_light       = 0;
                PlotParameters.XLim           = [-1.2*Lscale 1.2*Lscale];
                PlotParameters.YLim           = [-1.2*Lscale 1.2*Lscale];
                PlotParameters.ZLim           = [-1.2*Lscale 1.2*Lscale];
                PlotParameters.FrameRateValue = 50;
                PlotParameters.ClosePrevious  = true;
                PlotParameters.Position       = [0.1300 0.1100 0.7750 0.8150]; %default value (normalized)
                PlotParameters.CameraRotationSpeed = 0;
    
                Linkage.PlotParameters = PlotParameters;
            else
    
                % Initializing Links Vector
                VLinks    = SorosimLink.empty(nargin,0);
                LinkNames = cell(nargin,1);
    
                for i=1:nargin
                    L = varargin{i};
                    if ~isa(L, 'SorosimLink')
                        error('Error--->Input must be of SorosimLink class');
                    else
                        VLinks(i)    = L;
                        LinkNames{i} = inputname(i);
                    end
                end
    
                prompt           = {'Enter the total number of Links:'};
                dlgtitle         = 'Number of Links';
                definput         = {num2str(nargin)};
                opts.Interpreter = 'tex';
                ans_act          = inputdlg(prompt,dlgtitle,[1 75],definput,opts);
                if isempty(ans_act)
                    return;
                end

                N                = str2double(ans_act{1});

                Linkage.VLinks = VLinks;
                Linkage.N = 0;
                Linkage.PlotParameters.Lscale = 1;
                Linkage.ndof = 0;
                Linkage.nCLj = 0;
                Linkage.OneBasis = false; %default enable to use POD basis
    
                for i=1:N
                    
                    save('LinkageProgress.mat','Linkage')
    
                    UILinkInput(Linkage,LinkNames,i);
                    load('Temp_LinkageAssembly.mat','g_ini_i','LinkIndex_i','iLpre_i')
    
                    Linkage.LinkIndex = [Linkage.LinkIndex;LinkIndex_i];
                    Linkage.iLpre     = [Linkage.iLpre;iLpre_i];
                    Linkage.g_ini     = [Linkage.g_ini;g_ini_i];
                    Linkage.N         = Linkage.N+1;
    
                    Lscale = 0;
                    for iCLj=1:i
                        Lscale = Lscale+VLinks(Linkage.LinkIndex(iCLj)).L;
                    end
                    if Lscale==0
                        Lscale=1;
                    end
                    Linkage.PlotParameters.Lscale = Lscale;
    
                    ndof       = Linkage.ndof;
                    VRods_i    = SorosimRod.empty(Linkage.VLinks(LinkIndex_i).npie,0);
                    R          = jointRod_pre(Linkage.VLinks(LinkIndex_i),i); %for each joint and rigid link
                    VRods_i(1) = R;
                    ndof       = ndof+R.dof;
    
                    for j=1:Linkage.VLinks(LinkIndex_i).npie-1 %for each of the soft link divisions
    
                        R            = SorosimRod;
                        VRods_i(j+1) = R;
    
                    end
    
                    Linkage.CVRods = [Linkage.CVRods;{VRods_i}];
                    Linkage.ndof     = ndof;
    
                    if Linkage.VLinks(Linkage.LinkIndex(i)).jointtype~='N'
                        close all
                        Linkage.plotq0(i); % Why are we plotting?
                    end
    
                    happy = 0;
                    while ~happy
    
                        R      = jointRod(Linkage.VLinks(Linkage.LinkIndex(i)),i); %for each joint 
    
                        VRods_i(1) = R;
    
                        Linkage.CVRods{i} = VRods_i;
                        for j=1:Linkage.VLinks(Linkage.LinkIndex(i)).npie-1 %for each of the soft link divisions
                            close all
                            Linkage.plotq0(i,j);
                            Linkage.ndof         = Linkage.ndof-VRods_i(j+1).dof;
                            R               = SorosimRod(i,j,Linkage.VLinks(Linkage.LinkIndex(i))); %calling Rod class to generate twist for each piece of soft link
                            VRods_i(j+1)  = R;
                            Linkage.ndof         = Linkage.ndof+R.dof;
                            Linkage.CVRods{i}  = VRods_i;
                        end
    
                        close all
                        Linkage.plotq0;
                        
                        if Linkage.VLinks(Linkage.LinkIndex(i)).linktype=='s'
                            if any(Linkage.VLinks(Linkage.LinkIndex(i)).jointtype=='NSF')
                                quest  = 'Confirm degrees of freedom and reference strain?';
                            elseif any(Linkage.VLinks(Linkage.LinkIndex(i)).jointtype=='RPHCA')
                                quest  = 'Confirm joint definition, degrees of freedom and reference strain?';
                            end
                            FG_ANS = MFquestdlg([0.5, 0.5], quest,['Confirm link ',num2str(i),'?'],'Yes','No','Yes');
                        else
                            if any(Linkage.VLinks(Linkage.LinkIndex(i)).jointtype=='NSF')
                                FG_ANS = 'Yes';
                            elseif any(Linkage.VLinks(Linkage.LinkIndex(i)).jointtype=='RPHUCA')
                                quest  = 'Confirm joint definition?';
                                FG_ANS = MFquestdlg([0.5, 0.5], quest,['Confirm link ',num2str(i),'?'],'Yes','No','Yes');
                            end
                        end
    
                        switch FG_ANS
                            case 'Yes'
                                happy=1;
                            case 'No'
                                happy=0;
                        end
    
                    end
                    
                end
                
                % total number of points at which quantities are computed
                [Linkage.nsig,Linkage.nj]   = find_nsig_nj(Linkage);
                
                save('LinkageProgress.mat','Linkage')
    
                %% closed loop joints
    
                prompt           = {'Enter number of independent closed chain joints:'};
                dlgtitle         = 'Closed Chain';
                definput         = {'0'};
                opts.Interpreter = 'tex';
                ans_act          = inputdlg(prompt,dlgtitle,[1 50],definput,opts);
                nCLj             = str2double(ans_act{1});
    
                if nCLj>0
    
                    Linkage.iCLA       = zeros(nCLj,1);
                    Linkage.iCLA       = zeros(nCLj,1);
                    Linkage.VRodsCLj = SorosimRod.empty(nCLj,0);
    
                    Linkage.gACLj = cell(nCLj,1);
                    Linkage.gBCLj = cell(nCLj,1);
    
                    Phi_p       = cell(nCLj,1);
                    Linkage.nCLj=0;
    
                    for iCL=1:nCLj
                        UIClosedLoop(Linkage,iCL);
                        load('Temp_LinkageClosedJoint.mat','iA','iB','RodAB','gACLj','gBCLj')
    
                        Linkage.nCLj            = Linkage.nCLj+1;
                        Linkage.iCLA(iCL)       = iA;
                        Linkage.iCLB(iCL)       = iB;
                        Linkage.VRodsCLj(iCL)   = RodAB;
                        Linkage.gACLj{iCL}      = gACLj;
                        Linkage.gBCLj{iCL}      = gBCLj;
    
                        Phi_here = RodAB.Phi;
                        dof_here = RodAB.dof;
    
                        Phi_piCL = eye(6); % Start with 6x6 identity matrix

                        for id = 1:dof_here % Loop through each column of Phi_here and remove matching columns from Phi_piCL
                            match_idx = find(all(Phi_piCL == Phi_here(:, id), 1), 1);% Find the matching column index in Phi_piCL
                            if ~isempty(match_idx)
                                % Remove the matching column
                                Phi_piCL(:, match_idx) = [];
                            end
                        end
                        Phi_p{iCL} = Phi_piCL;
    
                    end
    
                    Linkage.CLprecompute.Phi_p = Phi_p; %Perpendicular basis of the closed loop joint
                    Linkage.T_BS               = 0.01; %desired settling time constant (Baumgart)

                    %index corresponding to the start of lambda
    
                end
                close all
                Linkage.plotq0;
    
                save('LinkageProgress.mat','Linkage')
                %% Gravity
    
                quest  = 'Is the system subjected to gravity?';
                FG_ANS = questdlg(quest,'Gravity','Yes','No','Yes');
    
                switch FG_ANS
                    case 'No'
                        Linkage.Gravity = false;
                        Linkage.G       = [0 0 0 0 0 -9.81]';
                    case 'Yes'
                        Linkage.Gravity = true;
                        
                        list  = {'x', 'y', 'z', '-x', '-y','-z'};
                        G_dir = listdlg('PromptString','Direction of Gravity',...
                                'SelectionMode','single',...
                                'ListSize',[160 160],...
                                'ListString',list,...
                                'InitialValue',6);
    
                        if G_dir==1
                            G = [0 0 0 9.81 0 0]';
                        elseif G_dir==2
                            G = [0 0 0 0 9.81 0]';
                        elseif G_dir==3
                            G = [0 0 0 0 0 9.81]';
                        elseif G_dir==4
                            G = [0 0 0 -9.81 0 0]';
                        elseif G_dir==5
                            G = [0 0 0 0 -9.81 0]';
                        elseif G_dir==6
                            G = [0 0 0 0 0 -9.81]';
                        end
                        Linkage.G = G;
                end
    
                close all
                Linkage.plotq0;
                %% Concentrated Point Force Fp 

                prompt           = {'Number of point forces/moments (wrenches)'};
                dlgtitle         = 'Number';
                definput         = {'0'};
                opts.Interpreter = 'tex';
                opts.WindowStyle = 'Normal';
                Fp_ans           = inputdlg(prompt,dlgtitle,[1 50],definput,opts);

                np       = str2double(Fp_ans{1});
                Linkage.Fp_loc = cell(np,1);
                Linkage.Fp_vec = cell(np,1);
                Linkage.np     = 0;
                
                Linkage.LocalWrench = ones(1,np);
                
                for iCLj=1:np

                    quest  = ['Is wrench ', num2str(iCLj),' a local wrench?'];
                    Follow_ANS = questdlg(quest,'Local or Global','Yes','No','No');

                    switch Follow_ANS

                        case 'No'
                            Linkage.LocalWrench(iCLj)= false;
                        case 'Yes'
                            Linkage.LocalWrench(iCLj) = true;
                    end

                    vrgood=false;
                    while ~vrgood
                        prompt           = {'Link number (choose from figure) corresponding to the point of application of force/moment '...
                                            ,'Division number (keep 1 for Rigid Link)','Position [0-1], (only for Soft Link) '...
                                            ,'M_x (Nm) as a function of t (s)','M_y (Nm) as a function of t (s)','M_z (Nm) as a function of t (s)'...
                                            ,'F_x (N) as a function of t (s)','F_y (N) as a function of t (s)','F_z (N) as a function of t (s)'};
                        dlgtitle         = ['Point Wrench',num2str(iCLj)];
                        definput         = {num2str(Linkage.N),'1','1','0','0','0','0','1','0'};
                        opts.Interpreter = 'tex';
                        Fp_ans           = inputdlg(prompt,dlgtitle,[1 75],definput,opts);

                        linknum = eval(Fp_ans{1});
                        divnum  = eval(Fp_ans{2});
                        X1 = eval(Fp_ans{3}); 
                        if linknum<=Linkage.N&&linknum>0
                            if Linkage.VLinks(Linkage.LinkIndex(linknum)).linktype=='r'
                                divmax = 1;
                            else
                                divmax = Linkage.VLinks(Linkage.LinkIndex(linknum)).npie-1;
                            end
                            if divnum<=divmax&&divnum>0&&X1>=0&&X1<=1
                                vrgood=true;
                            end
                        end
                        if ~vrgood
                            uiwait(msgbox('WRONG INPUTS','Error','error'))
                        else
                            if Linkage.VLinks(Linkage.LinkIndex(linknum)).linktype=='s'
                                Linkage.CVRods{linknum}(divnum+1).Xadd = [Linkage.CVRods{linknum}(divnum+1).Xadd;X1];
                            end
                        end

                    end

                    Linkage.Fp_loc{iCLj} = [linknum divnum X1];

                    syms t;
                    Mxs = str2sym(Fp_ans{4});
                    Mys = str2sym(Fp_ans{5});
                    Mzs = str2sym(Fp_ans{6});
                    Fxs = str2sym(Fp_ans{7});
                    Fys = str2sym(Fp_ans{8});
                    Fzs = str2sym(Fp_ans{9});
                    Fps = [Mxs,Mys,Mzs,Fxs,Fys,Fzs]';

                    if any(has(Fps,t))
                        Linkage.Fp_vec{iCLj} = matlabFunction(Fps);
                    else
                        Linkage.Fp_vec{iCLj} = str2func(['@(t) [' num2str(double(Fps')) ']''']);
                    end
                    Linkage.np=Linkage.np+1;

                end
                Linkage.Fp_sig = PointWrenchPoints(Linkage);

                save('LinkageProgress.mat','Linkage')

                %% Actuation

                close all
                Linkage.plotq0;
                
    
                quest  = 'Is the system Actuated?';
                Answer = questdlg(quest,'Actuation','Yes','No','Yes');
    
                switch Answer
                    case 'No'
                        Linkage.Actuated = false;
                        Linkage.nact     = 0;
                        Linkage.n_jact   = 0;
                        Linkage.WrenchControlled  = true;
                    case 'Yes'
                        if exist('CablePoints.mat','file')
                            delete('CablePoints.mat')
                        end
                        n_sact=0;
    
                        [n_jact,i_jact,i_jactq,WrenchControlled,Bj1] = JointActuation(Linkage,false);
    
                        Linkage.n_jact          = n_jact;
                        Linkage.nact            = n_jact; %will update later
                        Linkage.i_jact          = i_jact;
                        Linkage.N_jact          = length(i_jact);
                        Linkage.i_jactq           = i_jactq;
                        Linkage.WrenchControlled  = WrenchControlled;
                        Linkage.Bj1              = Bj1;
                        Linkage.Actuated          = true;

                        [n_k,index_q_u,index_q_k,index_u_k,index_u_u] = JointConstraintPrecompute(Linkage);
                        [ij_act,i_sig_act] = ActuatedJointIndex(Linkage);

                        Linkage.ActuationPrecompute.n_k = n_k; %for length controlled soft also (future version)
                        Linkage.ActuationPrecompute.index_q_u = index_q_u;
                        Linkage.ActuationPrecompute.index_q_k = index_q_k;
                        Linkage.ActuationPrecompute.index_u_k = index_u_k; %change when soft actuators are added
                        Linkage.ActuationPrecompute.index_u_u = index_u_u; %for length controlled soft also (future version)
                        Linkage.ActuationPrecompute.ij_act = ij_act;
                        Linkage.ActuationPrecompute.i_sig_act = i_sig_act;
    
                        if n_jact>0
                            close all
                            Linkage.plotq0;
                        end


                        %Soft link actuation
                        if any(Linkage.linktype=='s')
    
                            if all(Linkage.jointtype=='N')
                                Answer2='Yes';
                            else
                                quest   = 'Is (Are) the soft link(s) actuated?';
                                Answer2 = questdlg(quest,'Actuation','Yes','No','No');
                            end
    
                            switch Answer2
    
                                case 'Yes'
    
                                    [n_sact,dc_fn,dcp_fn,dc,dcp,Sdiv,Ediv] = CableActuation(Linkage);
    
                                    Linkage.dc     = dc;
                                    Linkage.dcp    = dcp;
                                    Linkage.Sdiv   = Sdiv;
                                    Linkage.Ediv   = Ediv;
                                    CableFunction.dc_fn=dc_fn;
                                    CableFunction.dcp_fn=dcp_fn;
                                    Linkage.CableFunction = CableFunction;
                                    
                                    N = Linkage.N; %repeat but may need
                                    i_sact = cell(N, 1); %indexes of soft actuators for every link and division syntax i_sact{i}{j}
                                    for i = 1:N
                                        if Linkage.VLinks(Linkage.LinkIndex(i)).linktype=='s'
                                            ndiv = Linkage.VLinks(Linkage.LinkIndex(i)).npie-1; 
                                            i_sact{i} = cell(1, ndiv);
                                            for j = 1:ndiv
                                                actuator_indices = zeros(1,n_sact); % if present the index otherwise zero
                                                for ia = 1:n_sact
                                                    if ~isempty(dc{ia, i}) && ~isempty(dc{ia, i}{j})
                                                        actuator_indices(ia) = ia;
                                                    end
                                                end
                                                i_sact{i}{j} = nonzeros(actuator_indices)'; % Store the non-empty ia indices for the current i, j
                                            end
                                        end
                                    end
                                    Linkage.i_sact = i_sact;
    
                                case 'No'
                                    n_sact = 0;
                            end
                        end
    
                        Linkage.n_sact = n_sact;
                        Linkage.nact   = n_jact+n_sact;
                        Linkage.ActuationPrecompute.index_u_k = [Linkage.ActuationPrecompute.index_u_k, Linkage.n_jact+1:Linkage.nact]; %including index of known soft actuator input force
                end
                Linkage.CA  = false; %custom internal actuation force
                Linkage.CAI = false; %custom actuation strength

                save('LinkageProgress.mat','Linkage')
                %% Miscellaneous precomputations
                
                %Computational point index precomputations for closed chain joints
                i_sigA   = zeros(Linkage.nCLj,1);
                i_sigB   = zeros(Linkage.nCLj,1);
    
                for iCL=1:Linkage.nCLj
                    i_sigAiCL = 0;
                    iA = Linkage.iCLA(iCL);
                    iB = Linkage.iCLB(iCL);
                    for i=1:iA
                        i_sigAiCL = i_sigAiCL+1;
                        if Linkage.VLinks(Linkage.LinkIndex(i)).linktype=='r'
                            i_sigAiCL = i_sigAiCL+1;
                        end
                        for j=1:Linkage.VLinks(Linkage.LinkIndex(i)).npie-1
                            i_sigAiCL = i_sigAiCL+Linkage.CVRods{i}(j+1).nip;
                        end
                    end
                    i_sigA(iCL) = i_sigAiCL;
    
                    i_sigBiCL = 0;
                    for i=1:iB
                        i_sigBiCL = i_sigBiCL+1;
                        if Linkage.VLinks(Linkage.LinkIndex(i)).linktype=='r'
                            i_sigBiCL = i_sigBiCL+1;
                        end
                        for j=1:Linkage.VLinks(Linkage.LinkIndex(i)).npie-1
                            i_sigBiCL = i_sigBiCL+Linkage.CVRods{i}(j+1).nip;
                        end
                    end
                    i_sigB(iCL) = i_sigBiCL;
                end
                Linkage.CLprecompute.i_sigA = i_sigA;
                Linkage.CLprecompute.i_sigB = i_sigB;
    
    
                % total number of computational points and virtual rigid joints
                [nsig,nj]    = find_nsig_nj(Linkage);
                Linkage.nsig = nsig;
                Linkage.nj   = nj;

                %To avoid redundant constraitns
                Linkage = ClosedLoopRedundancyPrecompute(Linkage); 
                
                Linkage.CEF        = false;  %No custom external force by default
                Linkage.UnderWater = false;  %Enable if it is an underwater simulation (default logical 0)
                Linkage.Rho_water  = 1000;   %Water density by default 1000 kg/m^3
                Linkage.M_added    = zeros(6*nsig,6);  %Added mass for fluid simulations (by default zeros(6*nsig,6))
                Linkage.DL         = zeros(6*nsig,6);  %Drag-Lift Matrix (by default zeros(6*nsig,6))
    
                % Damping matrix (D)
                D   = findD(Linkage);
                Linkage.D = D;
                Linkage.Damped = true;
    
                % Stiffness matrix (K)
                K   = findK(Linkage);
                Linkage.K = K;
                
                Linkage.Z_order = 4; %the value should be 2 or 4

                save('LinkageProgress.mat','Linkage')
    
                %% Plot parameters
                Lscale = Linkage.PlotParameters.Lscale;
                PlotParameters.CameraPosition = [Lscale*2 -Lscale/2 Lscale/2];
                PlotParameters.CameraTarget   = [0 0 0];
                PlotParameters.CameraUpVector = [0 0 1];
                PlotParameters.Light          = true;
                PlotParameters.Az_light       = 0;
                PlotParameters.El_light       = 0;
                PlotParameters.XLim           = [-1.2*Lscale 1.2*Lscale];
                PlotParameters.YLim           = [-1.2*Lscale 1.2*Lscale];
                PlotParameters.ZLim           = [-1.2*Lscale 1.2*Lscale];
                PlotParameters.FrameRateValue = 50;
                PlotParameters.ClosePrevious  = true;
                PlotParameters.Position       = [0.1300 0.1100 0.7750 0.8150]; %default value (normalized)
                PlotParameters.CameraRotationSpeed = 0;
                PlotParameters.VideoResolution = 0.5;
                PlotParameters.Lscale = Lscale;
    
                Linkage.PlotParameters = PlotParameters;
                close all
                Linkage.plotq0;

                save('LinkageProgress.mat','Linkage')
    
                %%
                if exist('Temp_LinkageAssembly.mat', 'file')
                    delete('Temp_LinkageAssembly.mat')
                end
                if exist('Temp_LinkageClosedJoint.mat', 'file')
                    delete('Temp_LinkageClosedJoint.mat')
                end
                if exist('Base_properties.mat', 'file')
                    delete('Base_properties.mat')
                end
                if exist('cableactuation.mat', 'file')
                    delete('cableactuation.mat')
                end
    
            end
        end %Class constructor

        %% Methods and set-get functions
        g       = FwdKinematics(Linkage,q,i,j)           %to get the transformation matrix at every significant points (arranged as column array) i: link, j: division (j=0 for joints)
        J       = Jacobian(Linkage,q,i,j)                %to get the Jacobian at every significant points (arranged as column array)
        Jd      = Jacobiandot(Linkage,q,qd,i,j)          %to get the derivative of Jacobian at every significant points (arranged as column array)
        xi      = ScrewStrain(Linkage,q,i,j)              %to get the screw strain at every significant points (arranged as column array)
        eta     = ScrewVelocity(Linkage,q,qd,i,j)        %to get the screw velocity at every significant points (arranged as column array)
        D       = findD(Linkage);                    %to compute and get the generalized damping matrix
        K       = findK(Linkage)                        %to compute and get the generalized stiffness matrix
        Bq      = ActuationMatrix(Linkage,q);             %to get the generalized actuation matrix (custom actuation not included)
        M       = GeneralizedMassMatrix(Linkage,q)        %to get the generalized mass matrix
        C       = GeneralizedCoriolisMatrix(Linkage,q,qd) %to get the generalized coriolis matrix
        F       = GeneralizedExternalForce(Linkage,q,qd,t)  %to get the generalized external force matrix
        [t,qqd] = dynamics(Linkage,x0,dynamicActionopen,dynamicsOptions)            %for dynamic simulation
        [q,u,lambda] = statics(Linkage,x0,action,staticsOptions)            %for static simulation

        plotq0(Linkage,Lh,Dh,CLh)     %to plot the free body diagram of the linkage
        plotq(Linkage,q)              %to plot the state of the linkage for a given q
        plotqt(Linkage,t,qqd,options) %to get dynamic simulation video output for a given t (time array) and qqd (array of joint coordinates and their time derivatives)

        %--------------------------------------------------------------------------
        %GET FUNCTIONS FOR DEPENDENT PROPERTIES: Connect the properties of
        %SorosimLink and Rod class to SorosimLinkage and allows calling values of the
        %properties
        %--------------------------------------------------------------------------

        % Link get functions

        function v = get.linktype(S)
            v = [S.VLinks.linktype];
        end

        function v = get.jointtype(S)
            v=[S.VLinks.jointtype];
        end

        function v = get.CS(S)
            v=[S.VLinks.CS];
        end

        %------------------------------------------------------------------
    end %methods

    %% Set and property update functions

    methods % Property Update functions

        %Changing nCLj:
        function S=set.VLinks(S, value)
            S.VLinks = value;
            S = S.Update();
        end
        function S=set.CVRods(S, value)
            S.CVRods = value;
            S = S.Update();
        end
        
        function Linkage = Update(Linkage)
            if isempty(Linkage.Z_order) %updated last
                return
            end
            Linkage.ndof=0;
            for i=1:Linkage.N
                Linkage.ndof = Linkage.ndof+Linkage.CVRods{i}(1).dof;
                for j=1:Linkage.VLinks(Linkage.LinkIndex(i)).npie-1
                    Linkage.ndof = Linkage.ndof+Linkage.CVRods{i}(j+1).dof;
                end
            end
            Linkage.D = findD(Linkage);
            Linkage.K = findK(Linkage);

            [Linkage.nsig,Linkage.nj] = find_nsig_nj(Linkage);
            [~,~,i_jactq_new,~,Bj1_new] = JointActuation(Linkage,true);
            Linkage.i_jactq = i_jactq_new;
            Linkage.Bj1    = Bj1_new;

            [n_k,index_q_u,index_q_k,index_u_k,index_u_u] = JointConstraintPrecompute(Linkage);
            [ij_act,i_sig_act] = ActuatedJointIndex(Linkage);

            Linkage.ActuationPrecompute.n_k = n_k;
            Linkage.ActuationPrecompute.index_q_u = index_q_u;
            Linkage.ActuationPrecompute.index_q_k = index_q_k;
            Linkage.ActuationPrecompute.index_u_k = index_u_k;
            Linkage.ActuationPrecompute.index_u_u = index_u_u;
            Linkage.ActuationPrecompute.ij_act = ij_act;
            Linkage.ActuationPrecompute.i_sig_act = i_sig_act;

            Linkage.Fp_sig=PointWrenchPoints(Linkage);

            i_sigA   = zeros(Linkage.nCLj,1);
            i_sigB   = zeros(Linkage.nCLj,1);

            for iCL=1:Linkage.nCLj
                i_sigAiCL = 0;
                ia = Linkage.iCLA(iCL);
                iB = Linkage.iCLB(iCL);
                for i=1:ia
                    i_sigAiCL = i_sigAiCL+1;
                    if Linkage.VLinks(Linkage.LinkIndex(i)).linktype=='r'
                        i_sigAiCL = i_sigAiCL+1;
                    end
                    for j=1:Linkage.VLinks(Linkage.LinkIndex(i)).npie-1
                        i_sigAiCL = i_sigAiCL+Linkage.CVRods{i}(j+1).nip;
                    end
                end
                i_sigA(iCL) = i_sigAiCL;

                i_sigBiCL = 0;
                for i=1:iB
                    i_sigBiCL = i_sigBiCL+1;
                    if Linkage.VLinks(Linkage.LinkIndex(i)).linktype=='r'
                        i_sigBiCL = i_sigBiCL+1;
                    end
                    for j=1:Linkage.VLinks(Linkage.LinkIndex(i)).npie-1
                        i_sigBiCL = i_sigBiCL+Linkage.CVRods{i}(j+1).nip;
                    end
                end
                i_sigB(iCL) = i_sigBiCL;
            end
            Linkage.CLprecompute.i_sigA = i_sigA;
            Linkage.CLprecompute.i_sigB = i_sigB;

            %To avoid redundant constraitns
            Linkage = ClosedLoopRedundancyPrecompute(Linkage); 
            
            for ia=1:Linkage.n_sact
                for i=1:Linkage.N

                    dc_fni  = Linkage.CableFunction.dc_fn{ia,i};
                    dcp_fni = Linkage.CableFunction.dcp_fn{ia,i};
                    Xc = 0;
                    if Linkage.Ediv(ia,i)>0
                        for j=Linkage.Sdiv(ia,i):Linkage.Ediv(ia,i)
                            ld = Linkage.VLinks(Linkage.LinkIndex(i)).ld{j};
                            nip = Linkage.CVRods{i}(j+1).nip;
                            Xs = Linkage.CVRods{i}(j+1).Xs;
                            dcj  = zeros(3,nip);
                            dcpj = zeros(3,nip);
                            
                            dcj(:,1) = dc_fni(Xc);
                            if ~isempty(dcp_fni)
                                dcpj(:,1) = dcp_fni(Xc);
                            end
                            
                            for ii=2:nip
                                Xc = Xc+(Xs(ii)-Xs(ii-1))*ld;
                                dcj(:,ii) = dc_fni(Xc);
                                if ~isempty(dcp_fni)
                                    dcpj(:,ii) = dcp_fni(Xc);
                                end
                            end
                            Linkage.dc{ia,i}{j}=dcj;
                            Linkage.dcp{ia,i}{j}=dcpj;
                        end
                        
                    end
                end
            end
        end
        

    end % Property update methods
end %classdef
