%This toolbox allows the modelling and simulation of static and dynamic
%soft robot open chains
%Run this file to initialize the Rigid-Soft Robotics Toolbox
%Last modified by Anup Teejo Mathew - 20/05/2021

%% Basic operation instructions:
%You can start using this toolbox by defining links, in order to create a link
%define a variable as shown below:

% LinkName = Link

%You will then have to answer questions in the form of dialougue boxes
%In order to define the link's gemetric and material properties.

%After your link is created you can create a linkage which connects any
%links you've created previously. A linkage is a chain of 1 or more links 
%connected in series to form an open chain and can be defied as shown below:

% LinkageName = Linkage(LinkName1,LinkName2.....LinkNameN)

%Dialougue boxes will then be used to collect input that is used to define
%the linkage properties

%You can access any link or linkage properties by calling the property name
%as shown below:

% LinkName.LinkPropertyName

% LinkageName.LinkagePropertyName

%Similarly methods can also be used to evaluate properties outside of the 
%classconstructor after the object is created 

%LinkageName.MethodName(Input1...InputN)

%You can run a static simulation by using the below command:
 
% LinkageName.statics

%the static simulation's output is a vector of joint angles

%A dynamic simulation can be performed using:

%LinkageName.dynamics

%The two outputs of a dynamic simulation are the time vector (t) and a
%matrix containing the configuration (q) and velocity (qd) at every time 
%element of (t). It is a good practice to save the outputs of the dynamic
%simulation for easy access.

%The examples folder of the toolbox contains some saved linkages and links
%you can run simulations for 

%%%%%%%%%%%%%%%%%%%%%%%% CLASSES AND THEIR PROPERTIES AND METHODS%%%%%%%%%%%%%%%%%%%%%
%% 1. Link Class

% %General Properties
% 
% jointtype       %Type of joint used to connect the link (lumped DoF). (R) for Revolute,(P) for Prismatic, (H) for Helical, (U) for Universal, (C) for Cylindrical, (A) for Planar, (S) for Spherical, (F) for Free motion and (N) for Fixed
% linktype        %'s' for soft or 'r' for rigid
% CS              %Cross section: 'R' for rectangular 'C' for circular
% npie            %Number of pieces. For a rigid link, npie = 1. For a soft link, npie=1+number of divisions 
% nGauss          %Number of Gaussian Quadrature points + 2
% Xs              %Gaussian Quadrature points inc 0 and 1
% Ws              %Gaussian Quadrature weights
% 
% 
% %Geometric Properties
% 
% lp        %Length of each divisions of the link (soft link) [m]
% L         %Total length of the link [m]
% r         %Radius as function of X1 (X1=X/L, X1 varies from 0 to 1) [m]
% h         %Height as function of X1 (X1=X/L, X1 varies from 0 to 1) [m]
% w         %Width as function of X1 (X1=X/L, X1 varies from 0 to 1) [m]
% cx        %x coordinate of origin wrt. previous frame [m]
% cy        %y coordinate of origin wrt. previous frame [m]
% cz        %z coordinate of origin wrt. previous frame [m]
% 
% %Material
% 
% E         %Young's modulus [Pa]
% Poi       %Poisson's ratio [-]
% G         %Shear modulus [Pa]
% Eta       %Material Damping [Pa.s]
% Rho       %Density [kg/m^3]
% 
% Ms        %Screw inertia
% Es        %Screw stiffness
% Gs        %Screw Damping
% 
% %Plot Properties
% 
% n_l       %Number of cross sections per division.
% n_r       %Number of radial points if the cross section is circular
% color     %color of link

%% 2. Twist

% Bdof         %(6x1) array specifying the allowable DoFs of a soft piece. 1 if allowed 0 if not.
% Bodr         %(6x1) array specifying the order of allowed DoF (0: constant, 1: linear, 2: quadratic,...)
% B            %(6xdof) Base matrix calculated at lumped joints or ((6xnGauss)xdof) base matrices computed at every significant points of a soft division
% B_Z1         %Base calculated at first Zanna point (Xs+Z1*(delta(Xs)))
% B_Z2         %Base calculated at the second Zanna point (Xs+Z2*(delta(Xs)))
% xi_starfn    %Reference strain vector as a function of X
% xi_star      %(6x1) reference strain vector at the lumped joint or ((6xnGauss)x3) reference strain vectors computed at Gauss quadrature and Zannah collocation points
% dof          %Total degrees of freedom of the piece


%% 3. Linkage

%General Properties

% N              %Total number of Links in series
% ntot           %Total number of pieces
% ndof           %Total number of DOF
% n_sig          %Total number of points at which quantites are evaluated (significant points): N+N_rigid+sum(N_soft_div*nGauss_div)
% VLinks         %Link vector from base to tip
% Vtwists        %Twist vector from base to tip
% Ltot           %Total length of linkage
% g_ini          %Initial transformation matrix
% g0             %Cell element of the initial transformation matrix of a piece
% q_scale        %(ndofx1) array of multiplier for each joint coordinate
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %External Force Properties
% Gravity               %1 if gravity is present 0 if not
% G                     %Value of G
% 
% %Point forces/moments
% PointForce            %1 if point force/moment is present 0 if not
% np                    %Number of point forces
% Fp_loc                %Cell element of piece numbers corresponding to the point force/moment location (at the tip of soft link/center of mass of rigid link)
% Fp_vec                %Cell element with value of point forces/moments [Mx My Mz Fx Fy Fz]'
% 
% %Custom external force
% CEFP                  %1 if custom external force is present 0 if not (default value is 0)
% M_added               %For external force that depend on etadot
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Actuation Properties
% Actuated          %1 if actuated 0 if not
% nact              %Total number of actuators
% 
% %Joint actuation parameters
% Bqj1              %Actuation base for joints with dof = 1
% n_jact            %total number of joint actuators
% i_jact            %index of Links whos joints are actuated
% i_jactq           %q index of all active joints
% WrenchControlled  %1 if the dof is wrench controlled 0 if q controlled
% 
% %Thread-like actuator for soft links
% n_sact            %Number of soft link actuators
% dc                %Cell element of local cable position (0, yp, zp) at Gauss quadrature points of all active soft divisions
% dcp               %Cell element of space derivative of the local cable position (0,yp',zp')
% Sdiv              %Starting division
% Ediv              %Ending division
% Inside            %1 if the actuator is fully inside the linkage, 0 if not
% 
% %Custom actuation
% CAP               %1 if custom actation is present 0 if not (default value is 0)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Elastic Properties
% K               %Stiffness matrix Qspace
% Damped          %1 if the soft links are elastically damped 0 if not (default value is 1)
% D               %Damping matrix in Qspace
% 
% CP1       %Custom constant properties of linkage
% CP2
% CP3
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CablePoints       %For displaying the soft actuator path
% PlotParameters
% 
% %%%PlotParameters is a struct with following elements%%%
% 
% %CameraPosition     CameraPosition with respect to origin
% %CameraUpVector     Orientation of normal
% %CameraTarget       Target location 
% %Az_light           Light Azimuth wrt camera
% %El_light           Light Elevation wrt camera
% %X_lim              X limit [X_lower_limt X_upper_limit]
% %Y_lim              Y limit [Y_lower_limt Y_upper_limit]
% %Z_lim              Z limit [Z_lower_limt Z_upper_limit]
% %FrameRateValue     FrameRate for dyanmic plot
% %ClosePrevious      0 to not close previous image 1 to close


%%%%METHODS%%%%
% g0      = findg0(S);                          %to get the initial transformation matrix of all pieces
% g       = FwdKinematics(S,q);                 %to get the transformation matrix at every significant points (arranged as column array)
% J       = Jacobian(S,q);                      %to get the Jacobian at every significant points (arranged as column array)
% Jd      = Jacobiandot(S,q,qd);                %to get the derivative of Jacobian at every significant points (arranged as column array)
% eta     = ScrewVelocity(S,q,qd);              %to get the screw velocity at every significant points (arranged as column array)
% D       = findD(S);                           %to compute and get the generalized damping matrix
% K       = findK(S)                            %to compute and get the generalized stiffness matrix
% Bq      = ActuationMatrix(S,q);               %to get the generalized actuation matrix (custom actuation not included)
% M       = GeneralizedMassMatrix(S,q)          %to get the generalized mass matrix
% C       = GeneralizedCoriolisMatrix(S,q,qd)   %to get the generalized coriolis matrix
% F       = GeneralizedExternalForce(S,q)       %to get the generalized external force matrix
% [t,qqd] = dynamics(S);                        %for dynamic simulation
% [q,u]   = statics(S,qu0)                      %for static simulation
% 
% plotq0(S,f);       %to plot the free body diagram of the linkage
% plotq(S,q);        %to plot the state of the linkage for a given q
% plotqqd(S,t,qqd);  %to get dynamic simulation video output for a given t (time array) and qqd (array of joint coordinates and their time derivatives)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

%% startup
clc
clear variables

addpath('Basic functions')
addpath('Actuation')
addpath('Link files')
addpath('Linkage files')
addpath('Custom')
addpath('Static and dynamic files')


if exist('.\Basis_properties.mat','file')
    delete('Basis_properties.mat')
end
if exist('.\cableactuation.mat','file')
    delete('cableactuation.mat')
end
if exist('.\CablePoints.mat','file')
    delete('CablePoints.mat')
end


disp('Type LinkName=Link to create the links (joint and body) and follow the instructions')
disp('Type LinkageName=Linkage(LinkName1,LinkName2,...,LinkNameN) to create linkages (combination of links 1 to N) and follow the instructions')
disp('For static equilibrium problem type [q,u]=Linkagename.statics')
disp('For dynamics problem type [t,qqd] = Linkagename.dynamics')
