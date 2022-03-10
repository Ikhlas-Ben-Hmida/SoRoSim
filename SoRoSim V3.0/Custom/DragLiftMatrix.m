%Function that allows the user to define a constant drag/Lift matrix
%to be saved as a constant property and used for external force
%calculations.
%In order to define the DL matrix the user can use the command:
% LinkageName.CustomProperty1=DragLiftMatrix(LinkageName)
%The user can specify different custom properties by changing or adding
%functions that follow the format of this file
%Last modified by Costanza Armanini 28/06/2021
function DL = DragLiftMatrix(Tr)

%-------------------------------------------------------------------------
% dynamic input environment

ro_water    = 1000;         % [Kg/m^3] water density

%-------------------------------------------------------------------------
% Geometrical input of body

Rb       = Tr.VLinks(1).r(0);                % [m] Radius
Lb       = Tr.VLinks(1).L;                 % [m] Length
Afr      = pi*Rb^2;                       % [m^2] frontale
Ala      = 2*Rb*Lb;                       % [m^2] laterale
Asu      = 2*pi*Rb^2+2*pi*Rb*Lb;          % [m^2] superficiale

Cmbx     = 2.5;                             % [-] coeff. viscosità longitudinale
Cmby     = 2.5;                             % [-] coeff. viscosità trasversale
Cmbz     = 2.5;                             % [-] coeff. viscosità trasversale
Clbx     = 2.5;                             % [-] coeff. viscosità longitudinale
Clby     = 2.5;                             % [-] coeff. viscosità trasversale
Clbz     = 2.5;                             % [-] coeff. viscosità trasversale
Db       = 0.5*ro_water*diag([Cmbx*Asu*Rb^3 Cmby*Ala*Rb^3 Cmbz*Ala*Rb^3 Clbx*Afr Clby*Ala Clbz*Ala]); % drag coef matrix


%-------------------------------------------------------------------------
% Geometrical input shaft

Rs       =Tr.VLinks(2).r(0);                % [m] Radius
Ls       =Tr.VLinks(2).L;                 % [m] Length

Clsx        =0.01;                                   % [-] Lift coeff. in x
Clsy        =2.5;                                    % [-] Lift coeff. in y
Clsz        =Clsy;                                    % [-] Lift coeff. in z
Ds          =0.5*ro_water*diag([0 0 0 2*pi*Rs*Ls*Clsx Rs*Ls*Clsy Rs*Ls*Clsz]);    % drag coef matrix


%-------------------------------------------------------------------------
% Geometrical input of Disk

Rb       = Tr.VLinks(3).r(0);                % [m] Radius
Lb       = Tr.VLinks(3).L;                 % [m] Length
Afr      = pi*Rb^2;                       % [m^2] frontale
Ala      = 2*Rb*Lb;                       % [m^2] laterale
Asu      = 2*pi*Rb^2+2*pi*Rb*Lb;          % [m^2] superficiale

Cmbx     = 2.5;                             % [-] coeff. viscosità longitudinale
Cmby     = 2.5;                             % [-] coeff. viscosità trasversale
Cmbz     = 2.5;                             % [-] coeff. viscosità trasversale
Clbx     = 2.5;                             % [-] coeff. viscosità longitudinale
Clby     = 2.5;                             % [-] coeff. viscosità trasversale
Clbz     = 2.5;                             % [-] coeff. viscosità trasversale
Dd       = 0.5*ro_water*diag([Cmbx*Asu*Rb^3 Cmby*Ala*Rb^3 Cmbz*Ala*Rb^3 Clbx*Afr Clby*Ala Clbz*Ala]); % drag coef matrix



%-------------------------------------------------------------------------
% Geometrical input  Hook

Rh       = Tr.VLinks(4).r{1}(0);                             % [m] Radius
nGauss   = Tr.VLinks(4).nGauss{1};

Cdhx        =0;                                     % [-] drag coeff. viscosità longitudinale
Cdhy        =1.1;                                    % [-] drag coeff. viscosità trasversale
Cdhz        =Cdhy;                                  % [-] drag coeff. viscosità trasversale
Clhx        =-0.1;                                    % [-] lift coeff.
Clhy        =0;                                     % [-] lift coeff.
Clhz        =0;                                     % [-] lift coeff.
Dh          =ro_water*([zeros(3,3) zeros(3,3); zeros(3,3) Rh*[0.5*pi*Cdhx -Clhz Clhy; Clhz Cdhy -Clhx; -Clhy Clhx Cdhz]]);
Dh          =repmat(Dh,nGauss,1);



%-------------------------------------------------------------------------
% Geometrical input  Hook2

Rh       = Tr.VLinks(4).r{1}(0);                             % [m] Radius
nGauss   = Tr.VLinks(4).nGauss{1};

Cdhx        =0;                                     % [-] drag coeff. viscosità longitudinale
Cdhy        =1.1;                                    % [-] drag coeff. viscosità trasversale
Cdhz        =Cdhy;                                  % [-] drag coeff. viscosità trasversale
Clhx        =-0.1;                                    % [-] lift coeff.
Clhy        =0;                                     % [-] lift coeff.
Clhz        =0;                                     % [-] lift coeff.
Dh2          =ro_water*([zeros(3,3) zeros(3,3); zeros(3,3) Rh*[0.5*pi*Cdhx -Clhz Clhy; Clhz Cdhy -Clhx; -Clhy Clhx Cdhz]]);
Dh2          =repmat(Dh2,nGauss,1);


%-------------------------------------------------------------------------
% Geometrical input  Hook3

Rh       = Tr.VLinks(4).r{1}(0);                             % [m] Radius
nGauss   = Tr.VLinks(4).nGauss{1};

Cdhx        =0;                                     % [-] drag coeff. viscosità longitudinale
Cdhy        =1.1;                                    % [-] drag coeff. viscosità trasversale
Cdhz        =Cdhy;                                  % [-] drag coeff. viscosità trasversale
Clhx        =-0.1;                                    % [-] lift coeff.
Clhy        =0;                                     % [-] lift coeff.
Clhz        =0;                                     % [-] lift coeff.
Dh3          =ro_water*([zeros(3,3) zeros(3,3); zeros(3,3) Rh*[0.5*pi*Cdhx -Clhz Clhy; Clhz Cdhy -Clhx; -Clhy Clhx Cdhz]]);
Dh3          =repmat(Dh3,nGauss,1);


%-------------------------------------------------------------------------
% Geometrical input of filament1
if Tr.VLinks(4).CS=='C'
    r      =Tr.VLinks(5).r{1};

    Cdfx        =Cdhx;                                      % [-] drag coeff. viscosità longitudinale
    Cdfy        =Cdhy;                                      % [-] drag coeff. viscosità trasversale
    Cdfz        =Cdhz ;                                     % [-] drag coeff. viscosità trasversale
    Clfx        =Clhx;                                      % [-] lift coeff.
    Clfy        =Clhy ;                                     % [-] lift coeff.
    Clfz        =Clhz;                                      % [-] lift coeff.
    Df          =zeros(6*nGauss,6);
    Xs     = Tr.VLinks(5).Xs{1};
    nGauss = Tr.VLinks(5).nGauss{1};

    for jj=1:nGauss
        Rf_here                              =r(Xs(jj));       % [m] radius definition based on equivalent segment areas
        Df(6*(jj-1)+1:6*(jj-1)+6,:)          =ro_water*([zeros(3,3) zeros(3,3); zeros(3,3) Rf_here*[0.5*pi*Cdfx -Clfz Clfy; Clfz Cdfy -Clfx; -Clfy Clfx Cdfz]]);
    end
else
    a      =Tr.VLinks(5).a{1};
    b      =Tr.VLinks(5).b{1};

    Cdfx        =Cdhx;                                      % [-] drag coeff. viscosità longitudinale
    Cdfy        =Cdhy;                                      % [-] drag coeff. viscosità trasversale
    Cdfz        =Cdhz ;                                     % [-] drag coeff. viscosità trasversale
    Clfx        =Clhx;                                      % [-] lift coeff.
    Clfy        =Clhy ;                                     % [-] lift coeff.
    Clfz        =Clhz;                                      % [-] lift coeff.
    Df          =zeros(6*nGauss,6);
    Xs     = Tr.VLinks(5).Xs{1};
    nGauss = Tr.VLinks(5).nGauss{1};

    for jj=1:nGauss
        Rf_here                              =sqrt((a(Xs(jj))^2+b(Xs(jj))^2)/2);       % [m] rms value used
        Df(6*(jj-1)+1:6*(jj-1)+6,:)          =ro_water*([zeros(3,3) zeros(3,3); zeros(3,3) Rf_here*[0.5*pi*Cdfx -Clfz Clfy; Clfz Cdfy -Clfx; -Clfy Clfx Cdfz]]);
    end
end


%-------------------------------------------------------------------------
% Geometrical input of filament2

if Tr.VLinks(4).CS=='C'
    r      =Tr.VLinks(5).r{1};

    Cdfx        =Cdhx;                                      % [-] drag coeff. viscosità longitudinale
    Cdfy        =Cdhy;                                      % [-] drag coeff. viscosità trasversale
    Cdfz        =Cdhz ;                                     % [-] drag coeff. viscosità trasversale
    Clfx        =Clhx;                                      % [-] lift coeff.
    Clfy        =Clhy ;                                     % [-] lift coeff.
    Clfz        =Clhz;                                      % [-] lift coeff.
    Df2          =zeros(6*nGauss,6);
    Xs     = Tr.VLinks(5).Xs{1};
    nGauss = Tr.VLinks(5).nGauss{1};

    for jj=1:nGauss
        Rf_here                              =r(Xs(jj));       % [m] radius definition based on equivalent segment areas
        Df2(6*(jj-1)+1:6*(jj-1)+6,:)          =ro_water*([zeros(3,3) zeros(3,3); zeros(3,3) Rf_here*[0.5*pi*Cdfx -Clfz Clfy; Clfz Cdfy -Clfx; -Clfy Clfx Cdfz]]);
    end
else
    a      =Tr.VLinks(5).a{1};
    b      =Tr.VLinks(5).b{1};

    Cdfx        =Cdhx;                                      % [-] drag coeff. viscosità longitudinale
    Cdfy        =Cdhy;                                      % [-] drag coeff. viscosità trasversale
    Cdfz        =Cdhz ;                                     % [-] drag coeff. viscosità trasversale
    Clfx        =Clhx;                                      % [-] lift coeff.
    Clfy        =Clhy ;                                     % [-] lift coeff.
    Clfz        =Clhz;                                      % [-] lift coeff.
    Df2          =zeros(6*nGauss,6);
    Xs     = Tr.VLinks(5).Xs{1};
    nGauss = Tr.VLinks(5).nGauss{1};

    for jj=1:nGauss
        Rf_here                              =sqrt((a(Xs(jj))^2+b(Xs(jj))^2)/2);       % [m] rms value used
        Df2(6*(jj-1)+1:6*(jj-1)+6,:)          =ro_water*([zeros(3,3) zeros(3,3); zeros(3,3) Rf_here*[0.5*pi*Cdfx -Clfz Clfy; Clfz Cdfy -Clfx; -Clfy Clfx Cdfz]]);
    end
end


%-------------------------------------------------------------------------
% Geometrical input of filament3

if Tr.VLinks(4).CS=='C'
    r      =Tr.VLinks(5).r{1};

    Cdfx        =Cdhx;                                      % [-] drag coeff. viscosità longitudinale
    Cdfy        =Cdhy;                                      % [-] drag coeff. viscosità trasversale
    Cdfz        =Cdhz ;                                     % [-] drag coeff. viscosità trasversale
    Clfx        =Clhx;                                      % [-] lift coeff.
    Clfy        =Clhy ;                                     % [-] lift coeff.
    Clfz        =Clhz;                                      % [-] lift coeff.
    Df3          =zeros(6*nGauss,6);
    Xs     = Tr.VLinks(5).Xs{1};
    nGauss = Tr.VLinks(5).nGauss{1};

    for jj=1:nGauss
        Rf_here                              =r(Xs(jj));       % [m] radius definition based on equivalent segment areas
        Df3(6*(jj-1)+1:6*(jj-1)+6,:)          =ro_water*([zeros(3,3) zeros(3,3); zeros(3,3) Rf_here*[0.5*pi*Cdfx -Clfz Clfy; Clfz Cdfy -Clfx; -Clfy Clfx Cdfz]]);
    end
else
    a      =Tr.VLinks(5).a{1};
    b      =Tr.VLinks(5).b{1};

    Cdfx        =Cdhx;                                      % [-] drag coeff. viscosità longitudinale
    Cdfy        =Cdhy;                                      % [-] drag coeff. viscosità trasversale
    Cdfz        =Cdhz ;                                     % [-] drag coeff. viscosità trasversale
    Clfx        =Clhx;                                      % [-] lift coeff.
    Clfy        =Clhy ;                                     % [-] lift coeff.
    Clfz        =Clhz;                                      % [-] lift coeff.
    Df3          =zeros(6*nGauss,6);
    Xs     = Tr.VLinks(5).Xs{1};
    nGauss = Tr.VLinks(5).nGauss{1};

    for jj=1:nGauss
        Rf_here                              =sqrt((a(Xs(jj))^2+b(Xs(jj))^2)/2);       % [m] rms value used
        Df3(6*(jj-1)+1:6*(jj-1)+6,:)          =ro_water*([zeros(3,3) zeros(3,3); zeros(3,3) Rf_here*[0.5*pi*Cdfx -Clfz Clfy; Clfz Cdfy -Clfx; -Clfy Clfx Cdfz]]);
    end
end


DL=[Db;Ds;Dd;Dh;Dh2;Dh3;Df;Df2;Df3];

end

