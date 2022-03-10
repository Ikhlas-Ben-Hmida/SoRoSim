%Function that allows the user to define a custom external force, which is
%a function of eta dot, to be added to the default linkage mass matrix.
%In order to define the added mass use the comand: LinkageName.M_added = AddedMass(LinkageName)
%The default value of the added mass is equal to zero, if no custm added mass is defined
function M_added = AddedMass(Tr)

%-------------------------------------------------------------------------
% dynamic input environment

ro_water    = 1000;         % [Kg/m^3] water density

%-------------------------------------------------------------------------
% Geometrical input of body

Rb       = Tr.VLinks(1).r(0);             % [m] Radius
Lb       = Tr.VLinks(1).L;                % [m] Length
Vol      = pi*Rb^2*Lb;                    % [m^3] volume cilinder

% Dynamic parameters
ecce     = 3/4;
alpha    = 2*(1-ecce^2)*...
    (0.5*log((1+ecce)/(1-ecce))-ecce)/...
    (ecce^3);                      % [-] added mass term
beta     = 1/ecce^2-(1-ecce^2)*log((1+ecce)/(1-ecce))/...
    (2*ecce^3);                    % [-] added mass term
gamma    = 9*Rb^2*(alpha-beta)/...
    (5*(5*(beta-alpha)-6));        % [m^2] added mass term

Maddb    = ro_water*Vol*diag([gamma 0 gamma beta/(2-beta) alpha/(2-alpha) beta/(2-beta)]);% added mass

%-------------------------------------------------------------------------
% Geometrical input of shaft

Rs       = Tr.VLinks(2).r(0);                % [m] Radius
Ls       = Tr.VLinks(2).L;                 % [m] Length
As       = pi*Rs^2;                       % [m^2] Area
Vs       = As*Ls;                         % [m^3] volume

% Dynamic parameters
Blsy     = 1.5;                                    % [-] Added mass coef. in y
Blsz     = 1.5;                                    % [-] Added mass coef. in z
Madds    = ro_water*Vs*diag([0 0 0 0 Blsy Blsz]);  % addedd mass matrix

%-------------------------------------------------------------------------
% Geometrical input of Disk

Rs       = Tr.VLinks(3).r(0);                % [m] Radius
Ls       = Tr.VLinks(3).L;                 % [m] Length
As       = pi*Rs^2;                       % [m^2] Area
Vs       = As*Ls;                         % [m^3] volume

% Dynamic parameters
Blsy     = 1.5;                                    % [-] Added mass coef. in y
Blsz     = 1.5;                                    % [-] Added mass coef. in z
Maddd    = ro_water*Vs*diag([0 0 0 0 Blsy Blsz]);  % addedd mass matrix

%-------------------------------------------------------------------------
% Geometrical input  Hook1
Rh       = Tr.VLinks(4).r{1}(0);             % [m] Radius
nGauss   = Tr.VLinks(4).nGauss{1};
Ah       = pi*Rh^2;                       % [m^2]

% Dynamic parameters
Blhy     = 0.6;                                         % [-] added mass coeff
Blhz     = Blhy;                                        % [-] added mass coeff
Maddh    = ro_water*diag([0 0 0 0 Ah*Blhy Ah*Blhz]);    % added mass
Maddh    = repmat(Maddh,nGauss,1);

%-------------------------------------------------------------------------
% Geometrical input  Hook2
Rh       = Tr.VLinks(4).r{1}(0);             % [m] Radius
nGauss   = Tr.VLinks(4).nGauss{1};
Ah       = pi*Rh^2;                       % [m^2]

% Dynamic parameters
Blhy     = 0.6;                                         % [-] added mass coeff
Blhz     = Blhy;                                        % [-] added mass coeff
Maddh2    = ro_water*diag([0 0 0 0 Ah*Blhy Ah*Blhz]);    % added mass
Maddh2    = repmat(Maddh2,nGauss,1);

%-------------------------------------------------------------------------
% Geometrical input  Hook3
Rh       = Tr.VLinks(4).r{1}(0);             % [m] Radius
nGauss   = Tr.VLinks(4).nGauss{1};
Ah       = pi*Rh^2;                       % [m^2]

% Dynamic parameters
Blhy     = 0.6;                                         % [-] added mass coeff
Blhz     = Blhy;                                        % [-] added mass coeff
Maddh3    = ro_water*diag([0 0 0 0 Ah*Blhy Ah*Blhz]);    % added mass
Maddh3    = repmat(Maddh3,nGauss,1);

%-------------------------------------------------------------------------
% Geometrical input of filament1

if Tr.VLinks(4).CS=='C'
    r        = Tr.VLinks(5).r{1};
    nGauss   = Tr.VLinks(5).nGauss{1};
    Xs       = Tr.VLinks(5).Xs{1};
    Blfy     = Blhy;                                   % [-] added mass coeff.
    Blfz     = Blhz;                                   % [-] added mass coeff.
    Maddf    = zeros(6*nGauss,6);

    for jj=1:nGauss
        Rf_here                               = r(Xs(jj));
        Af_here                               = pi*Rf_here^2;
        Maddf(6*(jj-1)+1:6*(jj-1)+6,:)        = ro_water*diag([0 0 0 0 Af_here*Blfy Af_here*Blfz]);
    end
else
    a        = Tr.VLinks(5).a{1};
    b        = Tr.VLinks(5).b{1};
    nGauss   = Tr.VLinks(5).nGauss{1};
    Xs       = Tr.VLinks(5).Xs{1};
    Blfy     = Blhy;                                   % [-] added mass coeff.
    Blfz     = Blhz;                                   % [-] added mass coeff.
    Maddf    = zeros(6*nGauss,6);

    for jj=1:nGauss
        Af_here                               = pi*a(Xs(jj))*b(Xs(jj));
        Maddf(6*(jj-1)+1:6*(jj-1)+6,:)        = ro_water*diag([0 0 0 0 Af_here*Blfy Af_here*Blfz]);
    end
end


%-------------------------------------------------------------------------
% Geometrical input of filament2

if Tr.VLinks(4).CS=='C'
    r        = Tr.VLinks(5).r{1};
    nGauss   = Tr.VLinks(5).nGauss{1};
    Xs       = Tr.VLinks(5).Xs{1};
    Blfy     = Blhy;                                   % [-] added mass coeff.
    Blfz     = Blhz;                                   % [-] added mass coeff.
    Maddf2    = zeros(6*nGauss,6);

    for jj=1:nGauss
        Rf_here                               = r(Xs(jj));
        Af_here                               = pi*Rf_here^2;
        Maddf2(6*(jj-1)+1:6*(jj-1)+6,:)        = ro_water*diag([0 0 0 0 Af_here*Blfy Af_here*Blfz]);
    end
else
    a        = Tr.VLinks(5).a{1};
    b        = Tr.VLinks(5).b{1};
    nGauss   = Tr.VLinks(5).nGauss{1};
    Xs       = Tr.VLinks(5).Xs{1};
    Blfy     = Blhy;                                   % [-] added mass coeff.
    Blfz     = Blhz;                                   % [-] added mass coeff.
    Maddf2   = zeros(6*nGauss,6);

    for jj=1:nGauss
        Af_here                               = pi*a(Xs(jj))*b(Xs(jj));
        Maddf2(6*(jj-1)+1:6*(jj-1)+6,:)        = ro_water*diag([0 0 0 0 Af_here*Blfy Af_here*Blfz]);
    end
end

%-------------------------------------------------------------------------
% Geometrical input of filament3

if Tr.VLinks(4).CS=='C'
    r        = Tr.VLinks(5).r{1};
    nGauss   = Tr.VLinks(5).nGauss{1};
    Xs       = Tr.VLinks(5).Xs{1};
    Blfy     = Blhy;                                   % [-] added mass coeff.
    Blfz     = Blhz;                                   % [-] added mass coeff.
    Maddf3    = zeros(6*nGauss,6);

    for jj=1:nGauss
        Rf_here                               = r(Xs(jj));
        Af_here                               = pi*Rf_here^2;
        Maddf3(6*(jj-1)+1:6*(jj-1)+6,:)        = ro_water*diag([0 0 0 0 Af_here*Blfy Af_here*Blfz]);
    end
else
    a        = Tr.VLinks(5).a{1};
    b        = Tr.VLinks(5).b{1};
    nGauss   = Tr.VLinks(5).nGauss{1};
    Xs       = Tr.VLinks(5).Xs{1};
    Blfy     = Blhy;                                   % [-] added mass coeff.
    Blfz     = Blhz;                                   % [-] added mass coeff.
    Maddf3    = zeros(6*nGauss,6);

    for jj=1:nGauss
        Af_here                               = pi*a(Xs(jj))*b(Xs(jj));
        Maddf3(6*(jj-1)+1:6*(jj-1)+6,:)        = ro_water*diag([0 0 0 0 Af_here*Blfy Af_here*Blfz]);
    end
end

M_added  = [Maddb;Madds;Maddd;Maddh;Maddh2;Maddh3;Maddf;Maddf2;Maddf3];

end

