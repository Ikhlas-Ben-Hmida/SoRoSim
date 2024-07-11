%Creates the SorosimLinkage class by joining SorosimLinks
%Last modified by Anup Teejo Mathew 02.03.2022
classdef SorosimLinkage
    properties% Independent properties of Linkage class
        %%
        %General Properties

        N            %Total number of Links
        ndof         %Total number of DOF
        nsig         %Total number of points at which compulations are performed (significant points): N+N_rigid+sum(N_soft_div*nGauss_div)
        VLinks       %Vector of all unique links (obtained from user input)
        LinkIndex    %(Nx1) array of indices corresponding to each links. ith Link = Tr.VLinks(LinkIndex(i))
        CVTwists     %Cell element of Twist vectors for each link

        iLpre        %(Nx1) array corresponding to the Link index of the Link to which the ith Link is connected
        g_ini        %(4Nx4) Fixed initial transformation matrices of Links wrt to the tip of its previous link
        Z_order      %order of Zannah collocation (2, 4, or 6) default value is 4

        %Closed Loop Joints (Link A is connect to Link B via a closed loop joint)

        nCLj         %Total number of closed loop joints
        iACL         %(nCLjx1)array corresponding to the index of Link A
        iCLB         %(nCLjx1)array corresponding to the index of Link B
        VTwistsCLj   %(nCLjx1)array of Twist vectors corresponding to each closed loop joint
        gACLj        %(nCLjx1)cells of fixed transformation from the tip of Link A to the close loop joint
        gBCLj        %(nCLjx1)cells of fixed transformation from the tip of Link B to the close loop joint
        CLprecompute %Struct element which contains pre-computed BpCLj (cell element of constrain basis), i_sigA (array of significant index corresponding to A), i_sigB (array ofsignificant index corresponding to B), and nCLp (total number of constraints)
        T_BS         %Baumgarte stabilization constant. Lower the value stricter the constrain.
        
        %External Force Properties
        Gravity       %logical 1 if gravity is present, 0 if not
        G             %Value of G

        %Point forces/moments
        PointForce    %logical 1 if point force/moment is present 0 if not
        FollowerForce %logical 1 if point force/moment is a follower force (local frame) and 0 if it is wrt global frame
        np            %Number of point forces
        Fp_loc        %Cell element of Link and Division numbers [i,j] corresponding to the point force/moment location (at the end of a soft division/center of mass of rigid link)
        Fp_vec        %Cell element with value of point forces/moments [Mx My Mz Fx Fy Fz]'

        %Custom external force
        CEFP          %logical 1 if custom external force is present 0 if not (default value is 0)
        M_added       %For the computation of external force that depend on etadot

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Actuation Properties
        Actuated          %logical 1 if actuated 0 if not
        nact              %Total number of actuators

        %Joint actuation parameters
        Bqj1              %Pre-computed actuation base for joints with dof = 1
        n_jact            %Total number of joint actuators
        i_jact            %(n_jactx1) array of index of links whos joints are actuated
        i_jactq           %(n_jactx1) array of joint coordinate index of all active joints
        WrenchControlled  %(n_jactx1) array of logical 1 if the joint is wrench controlled 0 if joint coordinate controlled

        %Thread-like actuator for soft links
        n_sact            %Number of soft link actuators
        dc                %(n_sactxN) cells of local cable position (0, yp, zp) at Gauss quadrature points of all active soft divisions
        dcp               %(n_sactxN) cells of space derivative of the local cable position (0,yp',zp')
        Sdiv              %(n_sactxN) cells of starting division number
        Ediv              %(n_sactxN) cells of ending division  number
        Inside            %(1xn_sact) cells with logical 1 if the actuator is fully inside the linkage, 0 if not
        CableFunction     %Struct with cell elements (Cy_fn and Cz_fn) of parameterized functions corresponding to the y and z coodinates of the cable

        %Custom actuation
        CAP               %logical 1 if custom actation is present 0 if not (default value: 0)
        CAS               %logical 1 to apply a custom actuator strength (default value: 0)

        %Pre-computed elastic Properties
        K       %Generalized Stiffness matrix
        Damped  %1 if the soft links are elastically damped 0 if not (default value is 1)
        D       %Generalized Damping matrix

        CP1     %Custom constant properties of linkage
        CP2
        CP3

        %Plotting
        PlotParameters

        %%%PlotParameters is a struct with following elements%%%
        %Lscale
        %CameraPosition     CameraPosition with respect to origin
        %CameraUpVector     Orientation of normal
        %CameraTarget       Target location
        %Light              logical 1 if the light is on, 0 if not. (default value: 1)
        %Az_light           Light Azimuth wrt camera
        %El_light           Light Elevation wrt camera
        %X_lim              X limit [X_lower_limt X_upper_limit]
        %Y_lim              Y limit [Y_lower_limt Y_upper_limit]
        %Z_lim              Z limit [Z_lower_limt Z_upper_limit]
        %FrameRateValue     FrameRate for dyanmic plot
        %ClosePrevious      logical 0 to not close previous image, 1 to close. (default value: 1)


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
        function Tr = SorosimLinkage(varargin)


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
            ans_act          = inputdlg(prompt,dlgtitle,[1 50],definput,opts);
            N                = str2num(ans_act{1});

            % Assigning values to properties of the class
            Tr.VLinks = VLinks;
            Tr.N = 0;
            Tr.PlotParameters.Lscale = 1;
            Tr.ndof = 0;

            for i=1:N

                UILinkInput(Tr,LinkNames,i);
                load('Temp_LinkageAssembly.mat','g_ini_i','LinkIndex_i','iLpre_i')

                Tr.LinkIndex = [Tr.LinkIndex;LinkIndex_i];
                Tr.iLpre     = [Tr.iLpre;iLpre_i];
                Tr.g_ini     = [Tr.g_ini;g_ini_i];
                Tr.N         = Tr.N+1;

                Lscale = 0;
                for ii=1:i
                    Lscale = Lscale+VLinks(Tr.LinkIndex(ii)).L;
                end
                if Lscale==0
                    Lscale=1;
                end
                Tr.PlotParameters.Lscale = Lscale;

                ndof         = Tr.ndof;
                VTwists_i    = SorosimTwist.empty(Tr.VLinks(LinkIndex_i).npie,0);
                T            = jointtwist_pre(Tr.VLinks(LinkIndex_i),i); %for each joint and rigid link
                VTwists_i(1) = T;
                ndof         = ndof+T.dof;

                for j=1:Tr.VLinks(LinkIndex_i).npie-1 %for each of the soft link divisions

                    T            = SorosimTwist;
                    VTwists_i(j+1) = T;

                end

                Tr.CVTwists = [Tr.CVTwists;{VTwists_i}];
                Tr.ndof     = ndof;

                if Tr.VLinks(Tr.LinkIndex(i)).jointtype~='N'
                    close all
                    Tr.plotq0(i); % Why are we plotting?
                end

                happy = 0;
                while ~happy

                    T            = jointtwist(Tr.VLinks(Tr.LinkIndex(i)),i); %for each joint and rigid link
                    VTwists_i(1) = T;

                    Tr.CVTwists{i} = VTwists_i;
                    for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie-1 %for each of the soft link divisions\
                        close all
                        Tr.plotq0(i,j);
                        Tr.ndof         = Tr.ndof-VTwists_i(j+1).dof;
                        T               = SorosimTwist(i,j,Tr.VLinks(Tr.LinkIndex(i))); %calling Twist class to generate twist for each piece of soft link
                        VTwists_i(j+1)  = T;
                        Tr.ndof         = Tr.ndof+T.dof;
                        Tr.CVTwists{i}  = VTwists_i;
                    end

                    close all
                    Tr.plotq0;
                    
                    if Tr.VLinks(Tr.LinkIndex(i)).linktype=='s'
                        if any(Tr.VLinks(Tr.LinkIndex(i)).jointtype=='NSF')
                            quest  = 'Confirm degrees of freedom and reference strain?';
                        elseif any(Tr.VLinks(Tr.LinkIndex(i)).jointtype=='RPHCA')
                            quest  = 'Confirm joint definition, degrees of freedom and reference strain?';
                        end
                        FG_ANS = MFquestdlg([0.5, 0.5], quest,['Confirm link ',num2str(i),'?'],'Yes','No','Yes');
                    else
                        if any(Tr.VLinks(Tr.LinkIndex(i)).jointtype=='NSF')
                            FG_ANS = 'Yes';
                        elseif any(Tr.VLinks(Tr.LinkIndex(i)).jointtype=='RPHUCA')
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
            nsig    = findnsig(Tr);
            Tr.nsig = nsig;
            


            %% closed loop

            prompt           = {'Enter number of closed loop joints:'};
            dlgtitle         = 'Closed loop joints';
            definput         = {'0'};
            opts.Interpreter = 'tex';
            ans_act          = inputdlg(prompt,dlgtitle,[1 50],definput,opts);
            nCLj             = str2num(ans_act{1});

            if nCLj>0

                Tr.iACL       = zeros(nCLj,1);
                Tr.iCLB       = zeros(nCLj,1);
                Tr.VTwistsCLj = SorosimTwist.empty(nCLj,0);

                Tr.gACLj = cell(nCLj,1);
                Tr.gBCLj = cell(nCLj,1);

                Bp       = cell(nCLj,1);
                Tr.nCLj=0;

                for iCL=1:nCLj
                    UIClosedLoop(Tr,iCL);
                    load('Temp_LinkageClosedJoint.mat','iA','iB','TwistAB','gACLj','gBCLj')

                    Tr.nCLj            = Tr.nCLj+1;
                    Tr.iACL(iCL)       = iA;
                    Tr.iCLB(iCL)       = iB;
                    Tr.VTwistsCLj(iCL) = TwistAB;
                    Tr.gACLj{iCL}      = gACLj;
                    Tr.gBCLj{iCL}      = gBCLj;

                    B_here   = TwistAB.B;
                    dof_here = TwistAB.dof;

                    BpiCL = eye(6);
                    for id=1:dof_here
                        flag=1;
                        k=1;
                        while flag
                            if isequal(B_here(:,id),BpiCL(:,k))
                                BpiCL(:,k)=[];
                                flag = 0;
                            end
                            k=k+1;
                        end
                    end
                    Bp{iCL} = BpiCL;

                end

                nCLp = 0;
                for ii=1:Tr.nCLj
                    nCLp = nCLp+size(Bp{ii},2); % can pre compute
                end
                Tr.CLprecompute.Bp     = Bp; %Perpendicular basis of the closed loop joint
                Tr.CLprecompute.nCLp   = nCLp; % ndof of Bp
                %sig point of tip of A and sig point of tip of B will be computed later after point load
                Tr.T_BS                = 0.01; %desired settling time constant

            end
            close all
            Tr.plotq0;
            %% External Force Vector F(q)

            % Gravity

            quest  = 'Is the system subjected to gravity?';
            FG_ANS = questdlg(quest,'Gravity','Yes','No','Yes');

            switch FG_ANS
                case 'No'
                    Tr.Gravity = false;
                    Tr.G       = 0;
                case 'Yes'
                    Tr.Gravity = true;
                    
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
                    Tr.G = G;
            end

            close all
            Tr.plotq0;
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Concentrated Point Force Fp

            quest  = 'Is the system subjected to an external point force/moment?';
            FP_ANS = questdlg(quest,'Point Force','Yes','No','No');

            switch FP_ANS

                case 'No'
                    Tr.PointForce = false;
                case 'Yes'

                    Tr.PointForce = true;

                    prompt           = {'Number of point forces/moments (wrenches)'};
                    dlgtitle         = 'Number';
                    definput         = {'1'};
                    opts.Interpreter = 'tex';
                    opts.WindowStyle = 'Normal';
                    Fp_ans           = inputdlg(prompt,dlgtitle,[1 50],definput,opts);

                    np        = str2num(Fp_ans{1});
                    Tr.Fp_loc = cell(np,1);
                    Tr.Fp_vec = cell(np,1);
                    Tr.np     = 0;
                    
                    Tr.FollowerForce = cell(1,np);
                    
                    for ii=1:np

                        quest  = ['Is wrench ', num2str(ii),' a follower wrench?'];
                        Follow_ANS = questdlg(quest,'Follower Force','Yes','No','No');

                        switch Follow_ANS

                            case 'No'
                                Tr.FollowerForce{ii} = false;
                            case 'Yes'
                                Tr.FollowerForce{ii} = true;
                        end

                        vrgood=false;
                        while ~vrgood
                            prompt           = {'Link number (choose from figure) corresponding to the point of application of force/moment '...
                                ,'Division number (keep 1 for Rigid Link)','Position [0-1], (only for Soft Link) '...
                                ,'M_x (Nm) as a function of t (s)','M_y (Nm) as a function of t (s)','M_z (Nm) as a function of t (s)','F_x (N) as a function of t (s)','F_y (N) as a function of t (s)','F_z (N) as a function of t (s)'};
                            dlgtitle         = ['Point Force/Moment (in global frame) ',num2str(ii)];
                            definput         = {num2str(Tr.N),'1','1','0','0','0','0','1','0'};
                            opts.Interpreter = 'tex';
                            Fp_ans           = inputdlg(prompt,dlgtitle,[1 75],definput,opts);

                            linknum = str2num(Fp_ans{1});
                            divnum  = str2num(Fp_ans{2});
                            X1 = str2num(Fp_ans{3}); 
                            if linknum<=Tr.N&&linknum>0
                                if Tr.VLinks(Tr.LinkIndex(linknum)).linktype=='r'
                                    divmax = 1;
                                else
                                    divmax = Tr.VLinks(Tr.LinkIndex(linknum)).npie-1;
                                end
                                if divnum<=divmax&&divnum>0&&X1>=0&&X1<=1
                                    vrgood=true;
                                end
                            end
                            if ~vrgood
                                uiwait(msgbox('WRONG INPUTS','Error','error'))
                            else
                                if Tr.VLinks(Tr.LinkIndex(linknum)).linktype=='s'
                                    Tr.CVTwists{linknum}(divnum+1).Xadd = [Tr.CVTwists{linknum}(divnum+1).Xadd;X1];
                                end
                            end

                        end

                        Tr.Fp_loc{ii} = [linknum divnum X1];

                        syms t;
                        Mxs = str2sym(Fp_ans{4});
                        Mys = str2sym(Fp_ans{5});
                        Mzs = str2sym(Fp_ans{6});
                        Fxs = str2sym(Fp_ans{7});
                        Fys = str2sym(Fp_ans{8});
                        Fzs = str2sym(Fp_ans{9});
                        Fps = [Mxs,Mys,Mzs,Fxs,Fys,Fzs]';

                        if any(has(Fps,t))
                            Tr.Fp_vec{ii} = matlabFunction(Fps);
                        else
                            Tr.Fp_vec{ii} = str2func(['@(t) [' num2str(double(Fps')) ']''']);
                        end
                        Tr.np=Tr.np+1;

                    end
                    close all
                    Tr.plotq0;
            end
            i_sigA   = zeros(nCLj,1);
            i_sigB   = zeros(nCLj,1);

            for iCL=1:nCLj
                i_sigAiCL = 0;
                iA = Tr.iACL(iCL);
                iB = Tr.iCLB(iCL);
                for i=1:iA
                    i_sigAiCL = i_sigAiCL+1;
                    if Tr.VLinks(Tr.LinkIndex(i)).linktype=='r'
                        i_sigAiCL = i_sigAiCL+1;
                    end
                    for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie-1
                        i_sigAiCL = i_sigAiCL+Tr.CVTwists{i}(j+1).nip;
                    end
                end
                i_sigA(iCL) = i_sigAiCL;

                i_sigBiCL = 0;
                for i=1:iB
                    i_sigBiCL = i_sigBiCL+1;
                    if Tr.VLinks(Tr.LinkIndex(i)).linktype=='r'
                        i_sigBiCL = i_sigBiCL+1;
                    end
                    for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie-1
                        i_sigBiCL = i_sigBiCL+Tr.CVTwists{i}(j+1).nip;
                    end
                end
                i_sigB(iCL) = i_sigBiCL;
            end
            Tr.CLprecompute.i_sigA = i_sigA;
            Tr.CLprecompute.i_sigB = i_sigB;
            % total number of points at which quantities are computed
            nsig    = findnsig(Tr);
            Tr.nsig = nsig;
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Custom external force

            Tr.CEFP     = false;

            %% Actuation

            quest  = 'Is the system Actuated?';
            Answer = questdlg(quest,'Actuation','Yes','No','Yes');

            switch Answer
                case 'No'
                    Tr.Actuated = false;
                    Tr.nact     = 0;
                    Tr.n_jact   = 0;
                case 'Yes'
                    if exist('CablePoints.mat','file')
                        delete('CablePoints.mat')
                    end
                    n_sact=0;

                    [n_jact,i_jact,i_jactq,WrenchControlled,Bqj1] = JointActuation(Tr,false);

                    Tr.n_jact            = n_jact;
                    Tr.i_jact            = i_jact;
                    Tr.i_jactq           = i_jactq;
                    Tr.WrenchControlled  = WrenchControlled;
                    Tr.Bqj1              = Bqj1;
                    Tr.Actuated          = true;


                    if n_jact>0
                        close all
                        Tr.plotq0;
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Soft link actuation
                    if any(Tr.linktype=='s')

                        if all(Tr.jointtype=='N')
                            Answer2='Yes';
                        else
                            quest   = 'Is (Are) the soft link(s) actuated?';
                            Answer2 = questdlg(quest,'Actuation','Yes','No','No');
                        end

                        switch Answer2

                            case 'Yes'

                                [n_sact,dc_fn,dcp_fn,dc,dcp,Sdiv,Ediv,Inside] = CableActuation(Tr);

                                Tr.dc     = dc;
                                Tr.dcp    = dcp;
                                Tr.Sdiv   = Sdiv;
                                Tr.Ediv   = Ediv;
                                Tr.Inside = Inside;
                                CableFunction.dc_fn=dc_fn;
                                CableFunction.dcp_fn=dcp_fn;
                                Tr.CableFunction = CableFunction;

                            case 'No'
                                n_sact = 0;
                        end
                    end

                    Tr.n_sact = n_sact;
                    Tr.nact   = n_jact+n_sact;
            end
            Tr.CAP = false;
            Tr.CAS = false;
            %% Constant coefficients
            
            % Damping matrix (K) estimation
            D     = findD(Tr);
            Tr.D   = D;
            Tr.Damped = true;
            % Stiffness matrix (K) estimation
            K    = findK(Tr);
            Tr.K  = K;

            
            Tr.Z_order   = 4; %the value should be 2 or 4
            Tr.nsig      = findnsig(Tr);
            %% Plot parameters
            PlotParameters.Lscale         = Lscale;
            PlotParameters.CameraPosition = [Lscale*2 -Lscale/2 Lscale/2];
            PlotParameters.CameraTarget   = [0 0 0];
            PlotParameters.CameraUpVector = [0 0 1];
            PlotParameters.Light          = true;
            PlotParameters.Az_light       = 0;
            PlotParameters.El_light       = 0;
            PlotParameters.X_lim          = [-1.2*Lscale 1.2*Lscale];
            PlotParameters.Y_lim          = [-1.2*Lscale 1.2*Lscale];
            PlotParameters.Z_lim          = [-1.2*Lscale 1.2*Lscale];
            PlotParameters.FrameRateValue = 50;
            PlotParameters.ClosePrevious  = true;
            PlotParameters.Position       = [0.1300 0.1100 0.7750 0.8150]; %default value (normalized)

            Tr.PlotParameters = PlotParameters;
            close all
            Tr.plotq0;

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


        end %Class constructor

        %% Methods and set-get functions
        g       = FwdKinematics(Tr,q,t,i,j);           %to get the transformation matrix at every significant points (arranged as column array) i: link, j: division (j=0 for joints)
        J       = Jacobian(Tr,q,t,i,j);                  %to get the Jacobian at every significant points (arranged as column array)
        Jd      = Jacobiandot(Tr,q,qd,t,i,j);          %to get the derivative of Jacobian at every significant points (arranged as column array)
        xi      = ScrewStrain(Tr,q,t,i,j)              %to get the screw strain at every significant points (arranged as column array)
        eta     = ScrewVelocity(Tr,q,qd,t,i,j);        %to get the screw velocity at every significant points (arranged as column array)
        D       = findD(Tr,q,qd,t);                    %to compute and get the generalized damping matrix
        K       = findK(Tr,q,t)                        %to compute and get the generalized stiffness matrix
        Bq      = ActuationMatrix(Tr,q,t);             %to get the generalized actuation matrix (custom actuation not included)
        M       = GeneralizedMassMatrix(Tr,q,t)        %to get the generalized mass matrix
        C       = GeneralizedCoriolisMatrix(Tr,q,qd,t) %to get the generalized coriolis matrix
        F       = GeneralizedExternalForce(Tr,q,qd,t)  %to get the generalized external force matrix
        [t,qqd] = dynamics(Tr,qqd0,odetype,dt);             %for dynamic simulation
        [q,u]   = statics(Tr,qu0,magnifier)            %for static simulation

        plotq0(Tr,Lh,Dh,CLh);       %to plot the free body diagram of the linkage
        plotq(Tr,q,t);              %to plot the state of the linkage for a given q
        plotqqd(Tr,t,qqd);          %to get dynamic simulation video output for a given t (time array) and qqd (array of joint coordinates and their time derivatives)

        %--------------------------------------------------------------------------
        %GET FUNCTIONS FOR DEPENDENT PROPERTIES: Connect the properties of
        %SorosimLink and Twist class to SorosimLinkage and allows calling values of the
        %properties
        %--------------------------------------------------------------------------

        % Link get functions

        function v = get.linktype(Tr)
            v = [Tr.VLinks.linktype];
        end

        function v = get.jointtype(Tr)
            v=[Tr.VLinks.jointtype];
        end

        function v = get.CS(Tr)
            v=[Tr.VLinks.CS];
        end

        %------------------------------------------------------------------
    end %methods

    %% Set and property update functions

    methods % Property Update functions

        %Changing nCLj:
        function Tr=set.VLinks(Tr, value)
            Tr.VLinks = value;
            Tr = Tr.Update();
        end
        function Tr=set.CVTwists(Tr, value)
            Tr.CVTwists = value;
            Tr = Tr.Update();
        end
        
        function Tr = Update(Tr)
            if isempty(Tr.Z_order)
                return
            end
            Tr.ndof=0;
            for i=1:Tr.N
                Tr.ndof = Tr.ndof+Tr.CVTwists{i}(1).dof;
                for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie-1
                    Tr.ndof = Tr.ndof+Tr.CVTwists{i}(j+1).dof;
                end
            end
            Tr.D       = findD(Tr);
            Tr.K       = findK(Tr);
            Tr.nsig    = findnsig(Tr);
            [~,~,i_jactq_new,~,Bqj1_new] = JointActuation(Tr,true);
            Tr.i_jactq           = i_jactq_new;
            Tr.Bqj1              = Bqj1_new;

            i_sigA   = zeros(Tr.nCLj,1);
            i_sigB   = zeros(Tr.nCLj,1);

            for iCL=1:Tr.nCLj
                i_sigAiCL = 0;
                ia = Tr.iACL(iCL);
                iB = Tr.iCLB(iCL);
                for i=1:ia
                    i_sigAiCL = i_sigAiCL+1;
                    if Tr.VLinks(Tr.LinkIndex(i)).linktype=='r'
                        i_sigAiCL = i_sigAiCL+1;
                    end
                    for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie-1
                        i_sigAiCL = i_sigAiCL+Tr.CVTwists{i}(j+1).nip;
                    end
                end
                i_sigA(iCL) = i_sigAiCL;

                i_sigBiCL = 0;
                for i=1:iB
                    i_sigBiCL = i_sigBiCL+1;
                    if Tr.VLinks(Tr.LinkIndex(i)).linktype=='r'
                        i_sigBiCL = i_sigBiCL+1;
                    end
                    for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie-1
                        i_sigBiCL = i_sigBiCL+Tr.CVTwists{i}(j+1).nip;
                    end
                end
                i_sigB(iCL) = i_sigBiCL;
            end
            Tr.CLprecompute.i_sigA = i_sigA;
            Tr.CLprecompute.i_sigB = i_sigB;
            
            for ia=1:Tr.n_sact
                for i=1:Tr.N

                    dc_fni  = Tr.CableFunction.dc_fn{ia,i};
                    dcp_fni = Tr.CableFunction.dcp_fn{ia,i};
                    Xc = 0;
                    if Tr.Ediv(ia,i)>0
                        for j=Tr.Sdiv(ia,i):Tr.Ediv(ia,i)
                            lp = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
                            nip = Tr.CVTwists{i}(j+1).nip;
                            Xs = Tr.CVTwists{i}(j+1).Xs;
                            dcj  = zeros(3,nip);
                            dcpj = zeros(3,nip);
                            
                            dcj(:,1) = dc_fni(Xc);
                            if ~isempty(dcp_fni)
                                dcpj(:,1) = dcp_fni(Xc);
                            end
                            
                            for ii=2:nip
                                Xc = Xc+(Xs(ii)-Xs(ii-1))*lp;
                                dcj(:,ii) = dc_fni(Xc);
                                if ~isempty(dcp_fni)
                                    dcpj(:,ii) = dcp_fni(Xc);
                                end
                            end
                            Tr.dc{ia,i}{j}=dcj;
                            Tr.dcp{ia,i}{j}=dcpj;
                        end
                        
                    end
                end
            end
        end
        

    end % Property update methods
end %classdef
