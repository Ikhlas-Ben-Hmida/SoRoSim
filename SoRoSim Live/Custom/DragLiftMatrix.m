%Function that allows the user to define a constant drag/Lift matrix
%to be saved as a constant property and used for external force
%calculations.

%DL is computed at all significant point except at joints

%In order to define the DL matrix the user can use the command:
%LinkageName.CP1=DragLiftMatrix(LinkageName)

%The user can specify different custom properties by changing or adding
%functions that follow the format of this file
%Last modified by Anup Teejo Mathew 09/02/2023
function DL = DragLiftMatrix(Tr)

n    = Tr.nsig-Tr.N; %everything except rigid joints
DL   = zeros(6*n,6);

%change from here
i = 1;
i_sig_nj = 1;

Rho_water    = 1000;         % [Kg/m^3] water density
%-------------------------------------------------------------------------
% Geometrical input of body
Rb = 0.09;

Afr       =pi*Rb^2;                          % [m^2] frontale
Ala       =(pi+4)*Rb^2;                      % [m^2] laterale
Asu       =8*pi*Rb^2;                        % [m^2] superficiale

Cmbx        =2.5;                                      % [-] coeff. viscosità 
Cmby        =0;                                      % [-] coeff. viscosità 
Cmbz        =Cmby;                                      % [-] coeff. viscosità 
Clbx        =1;                                      % [-] coeff. viscosità 
Clby        =1;                                      % [-] coeff. viscosità moto
Clbz        =Clbx;                                      % [-] coeff. viscosità 

DL((i_sig_nj-1)*6+1:i_sig_nj*6,:) = 0.5*Rho_water*diag([Cmbx*Ala*Rb^3 Cmby*Asu*Rb^3 Cmbz*Ala*Rb^3 Clbx*Ala Clby*Afr Clbz*Ala]); % drag coef matrix
i_sig_nj = i_sig_nj+1;

%-------------------------------------------------------------------------
% Geometrical input shaft

for i=2:5

    Rs = Tr.VLinks(Tr.LinkIndex(i)).r(0);              % [m] Radius
    Ls = Tr.VLinks(Tr.LinkIndex(i)).L;                 % [m] Length

    Clsx = 0.01;                                   % [-] Lift coeff. in x
    Clsy = 2.5;                                    % [-] Lift coeff. in y
    Clsz = Clsy;                                    % [-] Lift coeff. in z

    DL((i_sig_nj-1)*6+1:i_sig_nj*6,:) = 0.5*Rho_water*diag([0 0 0 2*pi*Rs*Ls*Clsx Rs*Ls*Clsy Rs*Ls*Clsz]);    % drag coef matrix
    i_sig_nj = i_sig_nj+1;
    
end

%% Soft Links

Cdx        =0;                                     % [-] drag coeff. viscosità longitudinale
Cdy        =1.1;                                   % [-] drag coeff. viscosità trasversale
Cdz        =Cdy;                                  % [-] drag coeff. viscosità trasversale
Clx        =-0.1;                                  % [-] lift coeff.
Cly        =0;                                     % [-] lift coeff.
Clz        =0;                                     % [-] lift coeff.

%same coefficients for hook and filament

for i=6:Tr.N-1 % 6 to 9 is hook, 10 to 13 is filament, i=N is Rod
    for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie-1
        r   = Tr.VLinks(Tr.LinkIndex(i)).r{j}; % [m] Radius
        nip = Tr.CVTwists{i}(j+1).nip;
        Xs  = Tr.CVTwists{i}(j+1).Xs;
        for jj=1:nip
            r_here =r(Xs(jj));       % [m] radius definition based on equivalent segment areas
            DL((i_sig_nj-1)*6+1:i_sig_nj*6,:) =Rho_water*([zeros(3,3) zeros(3,3); zeros(3,3) r_here*[0.5*pi*Cdx -Clz Cly; Clz Cdy -Clx; -Cly Clx Cdz]]);
            i_sig_nj = i_sig_nj+1;
        end
    end
end

%% Geometrical input of rod

DL((i_sig_nj-1)*6+1:i_sig_nj*6,:)  = zeros(6,6); % drag coef matrix

end

