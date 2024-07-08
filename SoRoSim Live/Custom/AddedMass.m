%Function that allows the user to define a custom external force, which is
%a function of eta dot, to be added to the default linkage mass matrix.

%M_added is computed at all significant point except at joints, Inertia
%metrix for rigid links, distributed for soft links

%In order to define the added mass use the comand: LinkageName.M_added = AddedMass(LinkageName)

%The default value of the added mass is []
function M_added = AddedMass(Tr)

%-------------------------------------------------------------------------
% dynamic input environment
nsig_nj  = Tr.nsig-Tr.N; %everything except rigid joints
M_added  = zeros(6*nsig_nj,6);

%change from here
i_sig_nj = 1;

Rho_water    = 1000;         % [Kg/m^3] water density
%-------------------------------------------------------------------------
Rb = 0.09;
Lcyl = 0.03;

Afr       =pi*Rb^2;                          % [m^2] frontale
Vsp       =4*pi*Rb^3/3;                      % [m^3] volume sfera
Vcyl      =Afr*Lcyl;                       % [m^3] volume cilindro
% drag lift body
a           = Rb+Lcyl/2;                         % [m] length ellipsoid (from the center)
b           = Rb;                                % [m] length ellipsoid (from the center)
ecce        = 1-(b/a)^2;                         % [-] eccentricita ellipsoid
alpha       = 2*(1-ecce^2)*...
              (0.5*log((1+ecce)/(1-ecce))-ecce)/...
              (ecce^3);                                  % [-] added mass term
beta        = 1/ecce^2-(1-ecce^2)*log((1+ecce)/(1-ecce))/...
              (2*ecce^3);                                % [-] added mass term

gamma       = ((b^2-a^2)^2*(alpha-beta))/...             % [m^2] added mass term
               5*(2*(b^2-a^2)+(b^2+a^2)*(beta-alpha));
M_added((i_sig_nj-1)*6+1:i_sig_nj*6,:) = Rho_water*(Vcyl+Vsp)*diag([gamma 0 gamma beta/(2-beta) alpha/(2-alpha) beta/(2-beta)]);     % massa da aggiungere
i_sig_nj = i_sig_nj+1;

%-------------------------------------------------------------------------
% Geometrical input of shaft

for i=2:5
    Rs       = Tr.VLinks(Tr.LinkIndex(i)).r(0);                % [m] Radius
    Ls       = Tr.VLinks(Tr.LinkIndex(i)).L;                 % [m] Length
    As       = pi*Rs^2;                       % [m^2] Area
    Vs       = As*Ls;                         % [m^3] volume

    % Dynamic parameters
    Blsy     = 1.5;                                    % [-] Added mass coef. in y
    Blsz     = 1.5;                                    % [-] Added mass coef. in z

    M_added((i_sig_nj-1)*6+1:i_sig_nj*6,:) = Rho_water*Vs*diag([0 0 0 0 Blsy Blsz]);  % addedd mass matrix
    i_sig_nj = i_sig_nj+1;
end

%% Soft Links
Bly     = 0.6;                                   % [-] added mass coeff.
Blz     = 0.6; 
%same coefficients for hook and filament
for i=6:Tr.N-1 % 6 to 9 is hook, 10 to 13 is filament, i=N is Rod
    for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie-1
        r   = Tr.VLinks(Tr.LinkIndex(i)).r{j}; % [m] Radius
        nip = Tr.CVTwists{i}(j+1).nip;
        Xs  = Tr.CVTwists{i}(j+1).Xs;
        for jj=1:nip
            r_here = r(Xs(jj));
            A_here = pi*r_here^2;
            M_added((i_sig_nj-1)*6+1:i_sig_nj*6,:) = Rho_water*diag([0 0 0 0 A_here*Bly A_here*Blz]);
            i_sig_nj = i_sig_nj+1;
        end
    end
end

%% Geometrical input of rod

M_added((i_sig_nj-1)*6+1:i_sig_nj*6,:)    = zeros(6,6);% added mass   
end

