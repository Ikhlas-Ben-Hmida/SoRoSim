%Creates the Linkage class by joining Links
%Last modified by Anup Teejo Mathew 21/05/2021
classdef Linkage
    properties% Independent properties of Linkage class
        %%
        %General Properties

        N              %Total number of Links in series
        ntot           %Total number of pieces
        ndof           %Total number of DOF
        n_sig          %Total number of points at which quantites are evaluated (significant points): N+N_rigid+sum(N_soft_div*nGauss_div)
        VLinks         %Link vector from base to tip
        Vtwists        %Twist vector from base to tip
        Ltot           %Total length of linkage
        g_ini          %Initial transformation matrix
        g0             %Cell element of the initial transformation matrix of a piece
        q_scale        %(ndofx1) array of multiplier for each joint coordinate

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %External Force Properties
        Gravity               %1 if gravity is present 0 if not
        G                     %Value of G
        
        %Point forces/moments
        PointForce            %1 if point force/moment is present 0 if not
        np                    %Number of point forces
        Fp_loc                %Cell element of piece numbers corresponding to the point force/moment location (at the tip of soft link/center of mass of rigid link)
        Fp_vec                %Cell element with value of point forces/moments [Mx My Mz Fx Fy Fz]'

        %Custom external force
        CEFP                  %1 if custom external force is present 0 if not (default value is 0)
        M_added               %For external force that depend on etadot

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Actuation Properties
        Actuated          %1 if actuated 0 if not
        nact              %Total number of actuators
        
        %Joint actuation parameters
        Bqj1              %Actuation base for joints with dof = 1
        n_jact            %total number of joint actuators
        i_jact            %index of Links whos joints are actuated
        i_jactq           %q index of all active joints
        WrenchControlled  %1 if the dof is wrench controlled 0 if q controlled
        
        %Thread-like actuator for soft links
        n_sact            %Number of soft link actuators
        dc                %Cell element of local cable position (0, yp, zp) at Gauss quadrature points of all active soft divisions
        dcp               %Cell element of space derivative of the local cable position (0,yp',zp')
        Sdiv              %Starting division
        Ediv              %Ending division
        Inside            %1 if the actuator is fully inside the linkage, 0 if not
        
        %Custom actuation
        CAP               %1 if custom actation is present 0 if not (default value is 0)
        CAS               %1 to apply a custom actuator strength (default value: 0)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Elastic Properties
        K               %Stiffness matrix Qspace
        Damped          %1 if the soft links are elastically damped 0 if not (default value is 1)
        D               %Damping matrix in Qspace
        
        CP1       %Custom constant properties of linkage
        CP2
        CP3
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        CablePoints       %For displaying the soft actuator path 
        PlotParameters
        
        %%%PlotParameters is a struct with following elements%%%

        %CameraPosition     CameraPosition with respect to origin
        %CameraUpVector     Orientation of normal
        %CameraTarget       Target location 
        %Az_light           Light Azimuth wrt camera
        %El_light           Light Elevation wrt camera
        %X_lim              X limit [X_lower_limt X_upper_limit]
        %Y_lim              Y limit [Y_lower_limt Y_upper_limit]
        %Z_lim              Z limit [Z_lower_limt Z_upper_limit]
        %FrameRateValue     FrameRate for dyanmic plot
        %ClosePrevious      0 to not close previous image 1 to close
  
        
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
        function S = Linkage(varargin)
            
            % Initializing Links Vector
            VLinks = Link.empty(nargin,0);
            for i=1:nargin
                L = varargin{i};
                if isa(L, 'Link')
                    VLinks(i) = L;
                else
                    disp('Error--->Input must be of Link class');
                end
            end
            % Assigning values to properties of the class
            S.VLinks = VLinks;
            N        = length(S.VLinks);
            S.N      = N;
            
            ntot = 0;
            Ltot = 0;
            for i=1:N
                ntot = ntot+S.VLinks(i).npie;
                Ltot = Ltot+S.VLinks(i).L;
            end
            S.ntot = ntot;
            S.Ltot = Ltot;
            
            % total number of points at which quantities are computed
            n_sig = N; %all joints
            for i=1:N
                for j=1:S.VLinks(i).npie-1
                    n_sig = n_sig+S.VLinks(i).nGauss{j};
                end
                if S.VLinks(i).linktype=='r'
                    n_sig = n_sig+1;
                end
            end
            S.n_sig = n_sig;

            %Pre-rotation of linkage:
            
            g_ini=eye(4);  
            S.g_ini = g_ini;
 
            %% Twist

            g0       = findg0(S);
            S.g0     = g0;
            Vtwists  = Twist.empty(S.ntot,0);
            
            happy=0;

            % To build the twist vector (Vtwists)
            while ~happy
                q_scale  = [];
                f       = 1; %index of piece 1 to ntot
                ndof    = 0;

                for i=1:N %for each link

                    T          = jointtwist_pre(S.VLinks(i),[]); %for each joint and rigid link
                    Vtwists(f) = T;
                    ndof       = ndof+T.dof;

                    f=f+1;

                    for j=1:S.VLinks(i).npie-1 %for each of the soft link divisions
                        
                        T          = Twist;
                        Vtwists(f) = T;
                        ndof       = ndof+T.dof;
                        f          = f+1;

                    end

                end %twist vector loop

                S.ndof    = ndof;
                S.Vtwists = Vtwists;
                
                f = 1; %index of piece 1 to ntot

                for i=1:N %for each link
                    
                    if S.VLinks(i).jointtype~='N'
                        close all
                        S.plotq0(f);
                    end
                    
                    T          = jointtwist(S.VLinks(i),i); %for each joint and rigid link
                    Vtwists(f) = T;
                    q_scale    = [q_scale; ones(T.dof,1)]; 
                    S.Vtwists  = Vtwists;

                    f=f+1;

                    for j=1:S.VLinks(i).npie-1 %for each of the soft link divisions

                        close all
                        S.plotq0(f);
                        T          = Twist(i,j,S.VLinks(i).Xs{j},S.VLinks(i).lp{j}); %calling Twist class to generate twist for each piece of soft link
                        Vtwists(f) = T;
                        ndof       = ndof+T.dof;

                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        Lscale = S.VLinks(i).lp{j};
                        for ii=1:6
                            if ii<4
                                if T.Bdof(ii)
                                    for kk=1:T.Bodr(ii)+1
                                        q_scale = [q_scale; 1/Lscale^kk];
                                    end
                                end
                            else
                                if T.Bdof(ii)
                                    for kk=1:T.Bodr(ii)+1
                                        q_scale = [q_scale; 1/Lscale^(kk-1)];
                                    end
                                end
                            end
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                        f=f+1;
                        S.ndof     = ndof;
                        S.Vtwists  = Vtwists;
                    end

                end

                close all
                S.plotq0;

                quest  = 'Confirm linkage geometry?';
%                 Options.WindowStyle='Normal';
                FG_ANS = MFquestdlg([0.5, 0.5], quest,'Confirmation','Yes','No','Yes');

                switch FG_ANS
                    case 'Yes'
                        happy=1;
                    case 'No'
                        happy=0;
                        close all
                end
            end
            S.q_scale = q_scale;          

            %% External Force Vector F(q)
            
            % Gravity
            
            quest  = 'Is the system subjected to gravity?';
            FG_ANS = questdlg(quest,'Gravity','Yes','No','Yes');
            
            switch FG_ANS
                case 'No'
                    S.Gravity = 0;
                    S.G       = 0;
                case 'Yes'
                    S.Gravity        = 1;
                    prompt           = {'Direction of gravity: (x, y, z, -x, -y or -z)'};
                    dlgtitle         = 'Gravitational Force';
                    definput         = {'-y'};
                    opts.Interpreter = 'tex';
                    G_d              = inputdlg(prompt,dlgtitle,[1 50],definput,opts);

                    G_dir = G_d{1};
                    
                    if G_dir=='x'
                        G = [0 0 0 9.81 0 0]';
                    elseif G_dir=='y'
                        G = [0 0 0 0 9.81 0]';
                    elseif G_dir=='z'
                        G = [0 0 0 0 0 9.81]';
                    elseif all(G_dir=='-x')
                        G = [0 0 0 -9.81 0 0]';
                    elseif all(G_dir=='-y')
                        G = [0 0 0 0 -9.81 0]';
                    elseif all(G_dir=='-z')
                        G = [0 0 0 0 0 -9.81]';
                    end
                    S.G = G;
            end
            
            close all
            S.plotq0;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Concentrated Point Force Fp
            
            quest  = 'Is the system subjected to an external point force/moment?';
            FP_ANS = questdlg(quest,'Point Force','Yes','No','No');
            
            switch FP_ANS
                
                case 'No'
                    S.PointForce = 0;
                case 'Yes'
                    
                    S.PointForce = 1;
                    
                    prompt           = {'Number of point forces/moments'};
                    dlgtitle         = 'Number';
                    definput         = {'1'};
                    opts.Interpreter = 'tex';
                    opts.WindowStyle = 'Normal';
                    Fp_ans           = inputdlg(prompt,dlgtitle,[1 50],definput,opts);

                    np     = str2num(Fp_ans{1});
                    Fp_loc = cell(1,np);
                    Fp_vec = cell(1,np);
                    S.np   = np;

                    for ii=1:np
                        goodpiece=0;
                        while ~goodpiece
                            prompt           = {['Piece number (choose from figure) corresponding to the point of application of force/moment ',num2str(ii)]...
                                                ,'M_x (Nm) as a function of t (s)','M_y (Nm) as a function of t (s)','M_z (Nm) as a function of t (s)','F_x (N) as a function of t (s)','F_y (N) as a function of t (s)','F_z (N) as a function of t (s)'};
                            dlgtitle         = 'Point Force/Moment';
                            definput         = {num2str(S.ntot),'0','0','0','0','1','0'};
                            opts.Interpreter = 'tex';
                            Fp_ans           = inputdlg(prompt,dlgtitle,[1 50],definput,opts);

                            piecenum   = str2num(Fp_ans{1});
                            f=1;

                            for i=1:N
                                if S.VLinks(i).linktype=='r'
                                    if piecenum==f
                                        goodpiece=1;
                                        break
                                    end
                                end
                                f=f+1;
                                for j=1:S.VLinks(i).npie-1
                                    if piecenum==f
                                        goodpiece=1;
                                        break
                                    end
                                    f=f+1;
                                end
                            end
                            if ~goodpiece
                                uiwait(msgbox('WRONG PIECE NUMBER','Error','error'))
                            end
                                
                        end

                        Fp_loc{ii} = piecenum;
                        
                        syms t;
                        Mxs = str2sym(Fp_ans{2});
                        Mys = str2sym(Fp_ans{3});
                        Mzs = str2sym(Fp_ans{4});
                        Fxs = str2sym(Fp_ans{5});
                        Fys = str2sym(Fp_ans{6});
                        Fzs = str2sym(Fp_ans{7});
                        Fps = [Mxs,Mys,Mzs,Fxs,Fys,Fzs]';

                        if any(has(Fps,t))
                            Fp_vec{ii} = matlabFunction(Fps);
                        else 
                            Fp_vec{ii} = str2func(['@(t) [' num2str(double(Fps')) ']''']);
                        end
                        
                        
                        
                        
                        
                        
%                         Mx         = str2num(Fp_ans{2});
%                         My         = str2num(Fp_ans{3});
%                         Mz         = str2num(Fp_ans{4});
%                         Fx         = str2num(Fp_ans{5});
%                         Fy         = str2num(Fp_ans{6});
%                         Fz         = str2num(Fp_ans{7});
% 
%                         Fp_vec{ii} = [Mx My Mz Fx Fy Fz]';

                    end
                    
                    S.Fp_loc = Fp_loc;
                    S.Fp_vec = Fp_vec;
                    close all
                    S.plotq0;
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Custom external force
                    
            S.CEFP     = 0;
            S.M_added  = zeros((n_sig-N)*6,6);
            
            %% Actuation 
            
            quest  = 'Is the system Actuated?';
            Answer = questdlg(quest,'Actuation','Yes','No','Yes');
            
            switch Answer
                case 'No'
                    S.Actuated = 0;
                    S.nact     = 0;
                    S.n_jact   = 0;
                case 'Yes'
                    if exist('CablePoints.mat','file')
                        delete('CablePoints.mat')   
                    end
                    n_sact=0;
                    
                    [n_jact,i_jact,i_jactq,WrenchControlled,Bqj1] = JointActuation(S);

                    S.n_jact            = n_jact;
                    S.i_jact            = i_jact;
                    S.i_jactq           = i_jactq;
                    S.WrenchControlled  = WrenchControlled;
                    S.Bqj1              = Bqj1;
                    S.Actuated=1;
                    
                    
                    if n_jact>0
                        close all
                        S.plotq0;
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Soft link actuation 
                    if any(S.linktype=='s')
                        
                        if all(S.jointtype=='N')
                            Answer2='Yes';
                        else
                            quest   = 'Is (Are) the soft link(s) actuated?';
                            Answer2 = questdlg(quest,'Actuation','Yes','No','No');
                        end

                        switch Answer2
                            
                            case 'Yes'

                                [n_sact,dc,dcp,Sdiv,Ediv,Inside] = CableActuation(S);
                                
                                S.dc     = dc;
                                S.dcp    = dcp;
                                S.Sdiv   = Sdiv;
                                S.Ediv   = Ediv;
                                S.Inside = Inside;
                                CablePoints   = load('CablePoints.mat');
                                S.CablePoints = CablePoints;
                                
                            case 'No'
                                n_sact = 0;
                        end
                    end
                    
                    S.n_sact = n_sact;
                    S.nact   = n_jact+n_sact;
            end
            S.CAP = 0;
            S.CAS = 0;
            %% Constant coefficients
            % Damping matrix (K) estimation
            D     = findD(S);
            S.D   = D;
            S.Damped = 1;
            % Stiffness matrix (K) estimation
            K    = findK(S);
            S.K  = K;
                     
            %% Plot parameters
            
            PlotParameters.CameraPosition = [Ltot Ltot/2 Ltot];
            PlotParameters.CameraTarget   = [0 0 0];
            PlotParameters.CameraUpVector = [0 1 0];
            PlotParameters.Az_light       = 0;
            PlotParameters.El_light       = 0;
            PlotParameters.X_lim          = [-1.2*Ltot 1.2*Ltot];
            PlotParameters.Y_lim          = [-1.2*Ltot 1.2*Ltot];
            PlotParameters.Z_lim          = [-1.2*Ltot 1.2*Ltot];
            PlotParameters.FrameRateValue = 100;
            PlotParameters.ClosePrevious  = 1;
            
            S.PlotParameters = PlotParameters;
            
            S.plotq0;

        end %Class constructor
        
        %% Methods and set-get functions
        g0      = findg0(S);                          %to get the initial transformation matrix of all pieces
        g       = FwdKinematics(S,q);                 %to get the transformation matrix at every significant points (arranged as column array)
        J       = Jacobian(S,q);                      %to get the Jacobian at every significant points (arranged as column array)
        Jd      = Jacobiandot(S,q,qd);                %to get the derivative of Jacobian at every significant points (arranged as column array)
        eta     = ScrewVelocity(S,q,qd);              %to get the screw velocity at every significant points (arranged as column array)
        D       = findD(S);                           %to compute and get the generalized damping matrix
        K       = findK(S)                            %to compute and get the generalized stiffness matrix
        Bq      = ActuationMatrix(S,q);               %to get the generalized actuation matrix (custom actuation not included)
        M       = GeneralizedMassMatrix(S,q)          %to get the generalized mass matrix
        C       = GeneralizedCoriolisMatrix(S,q,qd)   %to get the generalized coriolis matrix
        F       = GeneralizedExternalForce(S,q)       %to get the generalized external force matrix
        [t,qqd] = dynamics(S);                        %for dynamic simulation
        [q,u]   = statics(S,qu0)                      %for static simulation
        
        plotq0(S,f);       %to plot the free body diagram of the linkage
        plotq(S,q);        %to plot the state of the linkage for a given q
        plotqqd(S,t,qqd);  %to get dynamic simulation video output for a given t (time array) and qqd (array of joint coordinates and their time derivatives)
        
        %--------------------------------------------------------------------------
        %GET FUNCTIONS FOR DEPENDENT PROPERTIES: Connect the properties of
        %link and twist class to Linkage and allows calling values of the
        %properties
        %--------------------------------------------------------------------------

        % Link get functions

        function v = get.linktype(L1)
            v = [L1.VLinks.linktype];
        end
        
        function v = get.jointtype(L1)
            v=[L1.VLinks.jointtype];
        end
        
        function v = get.CS(L1)
            v=[L1.VLinks.CS];
        end
        
        %------------------------------------------------------------------
    end %methods
end %classdef
