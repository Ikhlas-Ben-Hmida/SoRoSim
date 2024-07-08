%This toolbox allows the modelling and simulation of static and dynamic
%Run this file to initialize the Rigid-Soft Robotics Toolbox
%Last modified by Anup Teejo Mathew - 18/01/2022

%% Basic operation instructions:
%You can start using this toolbox by defining links, in order to create a link
%define a variable as shown below:

% LinkName = SorosimLink

%You will then have to answer questions in the form of dialougue boxes
%In order to define the link's gemetric and material properties.

%After your link is created you can create a linkage which connects any
%links you've created previously. A linkage is a chain of 1 or more links.
%Syntax for creating a Linkage is given below:

% LinkageName = SorosimLinkage(LinkName1,LinkName2.....LinkNameN)

%Dialougue boxes will then be used to collect input that is used to define
%the linkage properties. You can access any link or linkage properties by calling the property name
%as shown below:

% LinkName.LinkPropertyName
% LinkageName.LinkagePropertyName

%Similarly methods can also be used to evaluate properties outside of the 
%classconstructor after the object is created 

%LinkageName.MethodName(Input1...InputN)

%You can run a static simulation by using the below command:
 
% LinkageName.statics
%The static simulation's output is a vector of joint angles

%A dynamic simulation can be performed using:

%LinkageName.dynamics
%The two outputs of a dynamic simulation are the time vector (t) and a
%matrix containing the configuration (q) and velocity (qd) at every time 
%element of (t). It is a good practice to save the outputs of the dynamic
%simulation for easy access.

%The examples folder of the toolbox contains some saved Linkages and Links
%you can run simulations for 

%%%%%%%%%%%%%%%%%%%%%%%% CLASSES AND THEIR PROPERTIES AND METHODS%%%%%%%%%%%%%%%%%%%%%
%% 1. SorosimLink Class

% %General Properties
        
% jointtype    %Type of joint used to connect the link (lumped DoF). (R) for Revolute,(P) for Prismatic, (H) for Helical, (C) for Cylindrical, (A) for Planar, (S) for Spherical, (F) for Free motion and (N) for Fixed
% linktype     %'s' for soft or 'r' for rigid
% CS           %Cross section: 'R' for rectangular 'C' for circular 'E' for elliptical
% npie         %Number of pieces. For a rigid link, npie = 1. For a soft link, npie=1+number of divisions 
% 
% %Geometric Properties
% 
% lp        %Length of each divisions of the link (soft link) [m]
% L         %Total length of the link [m]
% r         %Radius as a function of X1 (X1=X/L, X1 varies from 0 to 1) [m]
% h         %Height as a function of X1 (X1=X/L, X1 varies from 0 to 1) [m]
% w         %Width as a function of X1 (X1=X/L, X1 varies from 0 to 1) [m]
% a         %Semi-major axis as a function of X1 (X1=X/L, X1 varies from 0 to 1) [m]
% b         %Semi-minor axis as a function of X1 (X1=X/L, X1 varies from 0 to 1) [m]
% gi        %Transformation from joint to center of mass for ridig link to center of area for soft link
% gf        %Transformation to joint from center of mass for ridig link from center of area for soft link
% 
% %Material
% 
% E         %Young's modulus [Pa]
% Poi       %Poisson's ratio [-]
% G         %Shear modulus [Pa]
% Eta       %Material Damping [Pa.s]
% Rho       %Density [kg/m^3]
% Kj        %Joint Stiffness Matrix
% Dj        %Joint Damping Matrix
% 
% M        %Inertia matrix (only for rigid body)
% 
% %Plot Properties
% 
% n_l       %Number of cross sections per division. (default value: 10)
% n_r       %Number of radial points if the cross section is circular or ellipsoidal (default value: 18)
% color     %Color of link (random by default)
% 
% CPF       %Custom plot function for rigid bodies (logical 1 or 0)
% PlotFn    %Handle of function to plot the geometry (for rigid link)
% Lscale    %Scaling factor for plotting symbols or axes

%% 2. SorosimTwist

% Type         %Base type (Monomial, Lagrange Polynomial, Linear Interpolation, Gaussian, Custom, Non-linear Gaussian)
% SubClass     %Now only for FEM Like basis (linear,quadratic,cubic)
% Bdof         %(6x1) array specifying the allowable DoFs of a soft piece. 1 if allowed 0 if not.
% Bodr         %(6x1) array specifying the order of allowed DoF (0: constant, 1: linear, 2: quadratic,...)
% dof          %degress of freedom of each base
% 
% Bh           %Function handle for base
% B            %(6xdof) Base matrix calculated at lumped joints or ((6xnGauss)xdof) base matrices computed at every significant points of a soft division
% B_Z1         %Base calculated at 4th order first Zanna point (Xs+Z1*(delta(Xs)))
% B_Z2         %Base calculated at 4th order second Zanna point (Xs+Z2*(delta(Xs)))
% B_Z          %Base calculated at 2nd order Zanna point 
% 
% xi_starfn    %Reference strain vector as a function of X
% xi_star      %(6x1) reference strain vector at the lumped joint or ((6xnGauss)x4) reference strain vectors computed at Gauss quadrature and Zannah collocation points
% 
% Link         %Link associated with this twist only for soft link
% div          %Division associated with this twist
% nip           %number of integration point including boundaries
% Xs           %integration points 
% Ws           %weights of integration point
% Ms           %Inertia matrix of cross-section (6nip x 6) matrix
% Es           %Stiffness matrix (6nip x 6) matrix
% Gs           %Damping matrix (6nip x 6) matrix
% 
% Xadd         %additional integration points (nx1) vector 
% CI           %logical 0 by default 1 if custom integration is enabled
% CIFn         %function handle for custom integration

%%%%Methods%%%%
% Updatexi_star(T)
% UpdateMEG(T)
% UpdateAll(T)

%% 3. SorosimLinkage

% %General Properties
% 
% N            %Total number of Links
% ndof         %Total number of DOF
% nsig         %Total number of points at which compulations are performed (significant points): N+N_rigid+sum(N_soft_div*nGauss_div)
% VLinks       %Vector of all unique links (obtained from user input)
% LinkIndex    %(Nx1) array of indices corresponding to each links. ith Link = Tr.VLinks(LinkIndex(i))
% CVTwists     %Cell element of Twist vectors for each link
% 
% iLpre        %(Nx1) array corresponding to the Link index of the Link to which the ith Link is connected
% g_ini        %(4Nx4) Fixed initial transformation matrices of Links wrt to the tip of its previous link
% Z_order      %order of Zannah collocation (2, 4, or 6) default value is 4
% 
% %Closed Loop Joints (Link A is connect to Link B via a closed loop joint)
% 
% nCLj         %Total number of closed loop joints
% iACL         %(nCLjx1)array corresponding to the index of Link A
% iCLB         %(nCLjx1)array corresponding to the index of Link B
% VTwistsCLj   %(nCLjx1)array of Twist vectors corresponding to each closed loop joint
% gACLj        %(nCLjx1)cells of fixed transformation from the tip of Link A to the close loop joint
% gBCLj        %(nCLjx1)cells of fixed transformation from the tip of Link B to the close loop joint
% CLprecompute %Struct element which contains pre-computed BpCLj (cell element of constrain basis), i_sigA (array of significant index corresponding to A), i_sigB (array ofsignificant index corresponding to B), and nCLp (total number of constraints)
% T_BS         %Baumgarte stabilization constant. Lower the value stricter the constrain.
% 
% %External Force Properties
% Gravity       %logical 1 if gravity is present, 0 if not
% G             %Value of G
% 
% %Point forces/moments
% PointForce    %logical 1 if point force/moment is present 0 if not
% FollowerForce %logical 1 if point force/moment is a follower force (local frame) and 0 if it is wrt global frame
% np            %Number of point forces
% Fp_loc        %Cell element of Link and Division numbers [i,j] corresponding to the point force/moment location (at the end of a soft division/center of mass of rigid link)
% Fp_vec        %Cell element with value of point forces/moments [Mx My Mz Fx Fy Fz]'
% 
% %Custom external force
% CEFP          %logical 1 if custom external force is present 0 if not (default value is 0)
% M_added       %For the computation of external force that depend on etadot
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Actuation Properties
% Actuated          %logical 1 if actuated 0 if not
% nact              %Total number of actuators
% 
% %Joint actuation parameters
% Bqj1              %Pre-computed actuation base for joints with dof = 1
% n_jact            %Total number of joint actuators
% i_jact            %(n_jactx1) array of index of links whos joints are actuated
% i_jactq           %(n_jactx1) array of joint coordinate index of all active joints
% WrenchControlled  %(n_jactx1) array of logical 1 if the joint is wrench controlled 0 if joint coordinate controlled
% 
% %Thread-like actuator for soft links
% n_sact            %Number of soft link actuators
% dc                %(n_sactxN) cells of local cable position (0, yp, zp) at Gauss quadrature points of all active soft divisions
% dcp               %(n_sactxN) cells of space derivative of the local cable position (0,yp',zp')
% Sdiv              %(n_sactxN) cells of starting division number
% Ediv              %(n_sactxN) cells of ending division  number
% Inside            %(1xn_sact) cells with logical 1 if the actuator is fully inside the linkage, 0 if not
% CableFunction     %Struct with cell elements (Cy_fn and Cz_fn) of parameterized functions corresponding to the y and z coodinates of the cable
% 
% %Custom actuation
% CAP               %logical 1 if custom actation is present 0 if not (default value: 0)
% CAS               %logical 1 to apply a custom actuator strength (default value: 0)
% 
% %Pre-computed elastic Properties
% K       %Generalized Stiffness matrix
% Damped  %1 if the soft links are elastically damped 0 if not (default value is 1)
% D       %Generalized Damping matrix
% 
% CP1     %Custom constant properties of linkage
% CP2
% CP3
% 
% %Plotting
% PlotParameters
% 
% %%%PlotParameters is a struct with following elements%%%
% %Lscale
% %CameraPosition     CameraPosition with respect to origin
% %CameraUpVector     Orientation of normal
% %CameraTarget       Target location
% %Light              logical 1 if the light is on, 0 if not. (default value: 1)
% %Az_light           Light Azimuth wrt camera
% %El_light           Light Elevation wrt camera
% %X_lim              X limit [X_lower_limt X_upper_limit]
% %Y_lim              Y limit [Y_lower_limt Y_upper_limit]
% %Z_lim              Z limit [Z_lower_limt Z_upper_limit]
% %FrameRateValue     FrameRate for dyanmic plot
% %ClosePrevious      logical 0 to not close previous image, 1 to close. (default value: 1)


%%%%METHODS%%%%

% g       = FwdKinematics(Tr,q,t,i,j);           %to get the transformation matrix at every significant points (arranged as column array) i: link, j: division (j=0 for joints)
% J       = Jacobian(Tr,q,t,i,j);                %to get the Jacobian at every significant points (arranged as column array)
% Jd      = Jacobiandot(Tr,q,qd,t,i,j);          %to get the derivative of Jacobian at every significant points (arranged as column array)
% xi      = ScrewStrain(Tr,q,t,i,j)              %to get the screw strain at every significant points (arranged as column array)
% eta     = ScrewVelocity(Tr,q,qd,t,i,j);        %to get the screw velocity at every significant points (arranged as column array)
% D       = findD(Tr,q,qd,t);                    %to compute and get the generalized damping matrix
% K       = findK(Tr,q,t)                        %to compute and get the generalized stiffness matrix
% Bq      = ActuationMatrix(Tr,q,t);             %to get the generalized actuation matrix (custom actuation not included)
% M       = GeneralizedMassMatrix(Tr,q,t)        %to get the generalized mass matrix
% C       = GeneralizedCoriolisMatrix(Tr,q,qd,t) %to get the generalized coriolis matrix
% F       = GeneralizedExternalForce(Tr,q,qd,t)  %to get the generalized external force matrix
% [t,qqd] = dynamics(Tr,qqd0,odetype,dt);             %for dynamic simulation
% [q,u]   = statics(Tr,qu0,magnifier)            %for static simulation
% 
% plotq0(Tr,Lh,Dh,CLh);       %to plot the free body diagram of the linkage
% plotq(Tr,q,t);              %to plot the state of the linkage for a given q
% plotqqd(Tr,t,qqd);          %to get dynamic simulation video output for a given t (time array) and qqd (array of joint coordinates and their time derivatives)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

%% startup
clc
clear variables

addpath(genpath('Basic functions'))
addpath('Custom')
addpath('SorosimLink files')
addpath('SorosimTwist files')
addpath(genpath('SorosimLinkage files')) %include subfolders


if exist('.\Basis_properties.mat','file')
    delete('Basis_properties.mat')
end
if exist('.\cableactuation.mat','file')
    delete('cableactuation.mat')
end
if exist('.\CablePoints.mat','file')
    delete('CablePoints.mat')
end

disp('Welcome to SoRoSim Toolbox')
disp('Type LinkName=SorosimLink to create the links (joint and body)')
disp('Type LinkageName=SorosimLinkage(LinkName1,LinkName2,...,LinkNameN) to create linkages by combining links')
disp('For static equilibrium problem type [q,u]=LinkageName.statics')
disp('For dynamics problem type [t,qqd] = LinkageName.dynamics')
