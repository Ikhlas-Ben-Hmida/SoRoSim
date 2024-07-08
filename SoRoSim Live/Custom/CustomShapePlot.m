function CustomShapePlot(g)

%Use this function to plot a custom shape for any rigid link. g is the
%reference frame of the rigid body. Make sure to enable, the property CPF
%(Custom Plot Function) of the Link to true (1)

N_face = 12;
a=0.1; % zodiac edge length
phi=(1+sqrt(5))/2; %Golden ratio
Ra = a*phi/2; %Rc=Ra*sqrt(3) Rc is radius to corner
Rf = a*phi^2/(2*sqrt(3-phi)); % radius to frame

Points = zeros(3,20);
Points_New = zeros(3,20);
Points(:,1) = [1 1 -1];
Points(:,2) = [1 1 1];
Points(:,3) = [1 -1 -1];
Points(:,4) = [1 -1 1];
Points(:,5) = [-1 1 -1];
Points(:,6) = [-1 1 1];
Points(:,7) = [-1 -1 -1];
Points(:,8) = [-1 -1 1];
Points(:,9) = [0 phi 1/phi];
Points(:,10) = [0 phi -1/phi];
Points(:,11) = [0 -phi 1/phi];
Points(:,12) = [0 -phi -1/phi];
Points(:,13) = [1/phi 0 phi];
Points(:,14) = [1/phi 0 -phi];
Points(:,15) = [-1/phi 0 phi];
Points(:,16) = [-1/phi 0 -phi];
Points(:,17) = [phi 1/phi 0];
Points(:,18) = [phi -1/phi 0];
Points(:,19) = [-phi 1/phi 0];
Points(:,20) = [-phi -1/phi 0];
X = Points(1,:);
Y = Points(2,:);
Z = Points(3,:);
Points = Ra*Points;

Face1 = [5 10 1 14 16];
Face2 = [4 11 8 15 13];
Face3 = [12 7 16 14 3];
Face4 = [2 13 15 6 9];
Face5 = [5 19 20 7 16];
Face6 = [17 18 4 13 2];
Face7 = [5 19 6 9 10];
Face8 = [3 12 11 4 18];
Face9 = [10 9 2 17 1];
Face10 = [12 7 20 8 11];
Face11 = [1 17 18 3 14];
Face12 = [6 19 20 8 15];

%% ROTATION OF THE DODECAEDRON IN ORDER TO OBTAIN THE REFRERENCE FRAME THAT WE WANT

FaceX = Face2;
P1=[X(FaceX(1)) Y(FaceX(1)) Z(FaceX(1))];
P2=[X(FaceX(2)) Y(FaceX(2)) Z(FaceX(2))];
P3=[X(FaceX(3)) Y(FaceX(3)) Z(FaceX(3))];
v1=P1-P2;
v2=P1-P3;
C = cross(v1,v2);
u = C/sqrt(sum(C.^2));
th = acos(u(3)/(u(2)*u(2)+u(3)*u(3)));
R = [cos(th) sin(th);-sin(th) cos(th)];

T = eye(3);
T(2:3,2:3) = inv(R);

for i=1:length(Points)
    Points_New(:,i) = T*Points(:,i);
end
X = Points_New(1,:);
Y = Points_New(2,:);
Z = Points_New(3,:);
Pc=[mean(X(Face1)),mean(Y(Face1)),mean(Z(Face1));
    mean(X(Face2)),mean(Y(Face2)),mean(Z(Face2));
    mean(X(Face3)),mean(Y(Face3)),mean(Z(Face3));
    mean(X(Face4)),mean(Y(Face4)),mean(Z(Face4));
    mean(X(Face5)),mean(Y(Face5)),mean(Z(Face5));
    mean(X(Face6)),mean(Y(Face6)),mean(Z(Face6));
    mean(X(Face7)),mean(Y(Face7)),mean(Z(Face7));
    mean(X(Face8)),mean(Y(Face8)),mean(Z(Face8));
    mean(X(Face9)),mean(Y(Face9)),mean(Z(Face9));
    mean(X(Face10)),mean(Y(Face10)),mean(Z(Face10));
    mean(X(Face11)),mean(Y(Face11)),mean(Z(Face11));
    mean(X(Face12)),mean(Y(Face12)),mean(Z(Face12))];

lp = size(Points_New,2);
PointsH = [Points_New; ones(1,lp)];
Points_New = g*PointsH;
X = Points_New(1,:);
Y = Points_New(2,:);
Z = Points_New(3,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

patch(X([Face1 Face1(1)]),Y([Face1 Face1(1)]),Z([Face1 Face1(1)]),[.7 .7 .7],'EdgeColor','none')
hold on;patch(X([Face2 Face2(1)]),Y([Face2 Face2(1)]),Z([Face2 Face2(1)]),[.7 .7 .7],'EdgeColor','none')
patch(X([Face3 Face3(1)]),Y([Face3 Face3(1)]),Z([Face3 Face3(1)]),[.7 .7 .7],'EdgeColor','none')
patch(X([Face4 Face4(1)]),Y([Face4 Face4(1)]),Z([Face4 Face4(1)]),[.7 .7 .7],'EdgeColor','none')
patch(X([Face5 Face5(1)]),Y([Face5 Face5(1)]),Z([Face5 Face5(1)]),[.7 .7 .7],'EdgeColor','none')
patch(X([Face6 Face6(1)]),Y([Face6 Face6(1)]),Z([Face6 Face6(1)]),[.7 .7 .7],'EdgeColor','none')
patch(X([Face7 Face7(1)]),Y([Face7 Face7(1)]),Z([Face7 Face7(1)]),[.7 .7 .7],'EdgeColor','none')
patch(X([Face8 Face8(1)]),Y([Face8 Face8(1)]),Z([Face8 Face8(1)]),[.7 .7 .7],'EdgeColor','none')
patch(X([Face9 Face9(1)]),Y([Face9 Face9(1)]),Z([Face9 Face9(1)]),[.7 .7 .7],'EdgeColor','none')
patch(X([Face10 Face10(1)]),Y([Face10 Face10(1)]),Z([Face10 Face10(1)]),[.7 .7 .7],'EdgeColor','none')
patch(X([Face11 Face11(1)]),Y([Face11 Face11(1)]),Z([Face11 Face11(1)]),[.7 .7 .7],'EdgeColor','none')
patch(X([Face12 Face12(1)]),Y([Face12 Face12(1)]),Z([Face12 Face12(1)]),[.7 .7 .7],'EdgeColor','none')

light

% xlabel('x (mm)');ylabel('y (mm)');zlabel('z (mm)')
axis equal
% grid on;


%% OBTENGO LOS VECTORES UNITARIOS PERPENDICULARES A CADA CARA

% CARAS = [Face1;Face2;Face3;Face4;Face5;Face6;Face7;Face8;Face9;Face10;Face11;Face12]; 
Ui = zeros(N_face,3);
euli = zeros(N_face,3);
n_r=18;
n_l=3;
Xpatch  = zeros(n_r,(n_r-1)*(n_l-2)+2);
Ypatch  = zeros(n_r,(n_r-1)*(n_l-2)+2);
Zpatch  = zeros(n_r,(n_r-1)*(n_l-2)+2);


r     = 0.022;
theta = linspace(0,2*pi,n_r);
x     = zeros(1,n_r);
y     = r*sin(theta);
z     = r*cos(theta);
pos   = [x;y;z;ones(1,n_r)];


for i=1:N_face

Ui(i,:) = Pc(i,:)/norm(Pc(i,:));


if i==1
    eul = [0 -pi/2 0];
elseif i==2
    eul = [0 pi/2 0];
else
    eul = [0 0 0];
end
if i>2
    nn = floor((i+1)/2)-1;
    eul(3) = pi/2+(nn-1)*2*pi/5+pi*(rem(i,2)-1);
    eul(2) = atan(0.5)*(-1)^i;
end
euli(i,:) = eul;

xit1 = [eul(1) 0 0 0 0 0]';
xit2 = [0 eul(2) 0 0 0 0]';
xit3 = [0 0 eul(3) 0 0 0]';
xi2 = [ 0 0 0 Rf 0 0 ]';
g_here = g*variable_expmap_g(xit1)*variable_expmap_g(xit3)*variable_expmap_g(xit2)*variable_expmap_g(xi2); %note the value of eul (euler angles) do the sequence of transformations to find the frame of shafts like done here

%rotate first about z axis by eul(3) then about y axis by eul(2)

i_patch = 1;
pos_here = g_here*pos;
x_here   = pos_here(1,:);
y_here   = pos_here(2,:);
z_here   = pos_here(3,:);

Xpatch(:,i_patch) = x_here';
Ypatch(:,i_patch) = y_here';
Zpatch(:,i_patch) = z_here';
i_patch           = i_patch+1;

x_pre    = x_here;
y_pre    = y_here;
z_pre    = z_here;

g_here   = g_here*[eye(3) [0.05;0;0];[0 0 0 1]];
pos_here = g_here*pos;
x_here   = pos_here(1,:);
y_here   = pos_here(2,:);
z_here   = pos_here(3,:);

for jj=1:n_r-1

    Xpatch(1:5,i_patch)   = [x_pre(jj) x_here(jj) x_here(jj+1) x_pre(jj+1) x_pre(jj)]';
    Xpatch(6:end,i_patch) = x_pre(jj)*ones(n_r-5,1);
    Ypatch(1:5,i_patch)   = [y_pre(jj) y_here(jj) y_here(jj+1) y_pre(jj+1) y_pre(jj)]';
    Ypatch(6:end,i_patch) = y_pre(jj)*ones(n_r-5,1);
    Zpatch(1:5,i_patch)   = [z_pre(jj) z_here(jj) z_here(jj+1) z_pre(jj+1) z_pre(jj)]';
    Zpatch(6:end,i_patch) = z_pre(jj)*ones(n_r-5,1);
    i_patch = i_patch+1;

end

Xpatch(:,i_patch) = x_here';
Ypatch(:,i_patch) = y_here';
Zpatch(:,i_patch) = z_here';

patch(Xpatch,Ypatch,Zpatch,[0.3 0.3 0.3],'EdgeColor','none');

end

% close all
% save('EulerAngles','euli')
% for assembly, rotate first about z axis by eul(3) then about y axis by eul(2) thats all            

end