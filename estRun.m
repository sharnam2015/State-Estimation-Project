function [x,y,theta,internalStateOut] = estRun(time, dt,...
    internalStateIn, steeringAngle, pedalSpeed, measurement)

x = internalStateIn.x;
y = internalStateIn.y;
theta = internalStateIn.theta;
color = internalStateIn.color;
Pp = internalStateIn.P;

%%%%%%%%%%%%%%%% UKF %%%%%%%%%%%%%%%%%%%%%%
% Variables for modelling
r = 0.425;
B = 0.8;
v = 5*r*pedalSpeed;
s = zeros(3,6);

%% Noise settings
W = [1.0893,1.5333;1.5333,2.9880];
%%Defining variables
T = chol(Pp);
X = [x;y;theta];
%V = zeros(3,3);
sigmaV = 0.1;
Xm = zeros(3,1);
Pm = zeros(3,1);
delB = 0.1;
delr = 0.05;
delv = 5*pedalSpeed*delr;
deltheta = (delv*tan(steeringAngle)/B)-(v*tan(steeringAngle)*delB/(B*B));
delx = (delv*cos(theta))-(v*sin(theta)*deltheta);
dely = (delv*sin(theta))+(v*cos(theta)*deltheta);
vn = [delx;dely;deltheta];
%%Using different V matrices
%V = [delx*delx,delx*dely,delx*deltheta;
%    dely*delx,dely*dely,dely*deltheta;
%    delx*deltheta,deltheta*dely,deltheta*deltheta];
%V = [delx,0,0;0,dely,0;0,0,deltheta];
%vn
%V = zeros(3)
%sigmaV = 0.02125^2;
%sigmaB = 0.07968^2;
%V =dt*[sigmaV, 0, 0;
%    0, sigmaV, 0;
%    0, 0, sigmaV + sigmaB];

b = 0.002;
V = [b,0,0;0,b,0;0,0,b];
%V = [delx,0,0;0,dely,0;0,0,deltheta];
%V = zeros(3);

%%Initializing Sigma Points
s1 = X+(sqrt(3)*T(:,1));
s4 = X-(sqrt(3)*T(:,1));
s2 = X+(sqrt(3)*T(:,2));
s5 = X-(sqrt(3)*T(:,2));
s3 = X+(sqrt(3)*T(:,3));
s6 = X-(sqrt(3)*T(:,3));
s = [s1,s2,s3,s4,s5,s6];
sx = zeros(6,1);
sy = zeros(6,1);
stheta = zeros(6,1);

%Sending the points through the model
for i=1:6
    sx(i) = s(1,i) + (v*cos(theta))*dt;
    sy(i) = s(2,i) + (v*sin(theta))*dt;
    stheta(i) = s(3,i) + ((v*tan(steeringAngle))*dt/B);
end

usx = mean(sx);
usy = mean(sy);
utheta = mean(stheta);
usm = [usx;usy;utheta];
%sVr = [var(sx);var(sy);var(stheta)];
sus = [sx,sy,stheta];

for j = 1:6
    Pp = Pp + ((sus(j,:)')-usm)*transpose((sus(j,:)')-usm);
end
Pp = Pp/6;
Pp = Pp+V;

xp = usx;
yp = usy;
thetap = utheta;
szx = zeros(6,1);
szy = zeros(6,1);
sz = zeros(6,2);
Xp = [xp;yp;thetap];

%%Checking for measurement and if measurement is available then updating
if ~isnan(measurement(1)) & ~isnan(measurement(2))
    measurex = measurement(1);
    measurey = measurement(2);
    
    % Measurement passing sigma points through the z model
    for i = 1:1:6
        szx(i)= sx(i)+ 0.5*B*cos(stheta(i));
        szy(i) = sy(i) + 0.5*B*sin(stheta(i));
    end
    sz = [szx,szy];
    zex = [mean(szx);mean(szy)];
    Pzz = zeros(2,2);
    Pxz = zeros(3,2);
    Pxs = zeros(3,2);
    Pzs = zeros(2,2);
    
    % Calculating Pzz
    for j = 1:6
        Pzz = Pzs + (((sz(j,:)')-zex)*transpose(((sz(j,:)')-zex)));
        Pzs = Pzz;
        Pxz = Pxz + ((([sx(j);sy(j);stheta(j)])-Xp)*transpose(((sz(j,:)')-zex)));
        Pxs = Pxz;
    end
    Pzz = ((1/6)*(Pzz))+W;
    Pxz = ((1/6)*Pxz);
    K = zeros(3,2);
    K = Pxz*inv(Pzz);
    Z = zeros(2,1);
    Z = [measurex;measurey];
    Xm = Xp+(K*(Z-zex));
    Pm = Pp - (K*Pzz*(transpose(K)));
else
    Pm = Pp;
    Xm = Xp;
end

x = Xm(1);
y = Xm(2);
theta = Xm(3);
% we're unreliable about our favourite colour:
if strcmp(color, 'green')
    color = 'red';
else
    color = 'green';
end

%% OUTPUTS %%
% Update the internal state (will be passed as an argument to the function
% at next run), must obviously be compatible with the format of
internalStateOut.x = x;
internalStateOut.y = y;
internalStateOut.theta = theta;
internalStateOut.P = Pm;
internalStateOut.color = color;
end
