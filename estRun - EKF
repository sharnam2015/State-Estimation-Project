function [x,y,theta,internalStateOut] = estRun(time,...
    dt, internalStateIn, steeringAngle, pedalSpeed, measurement)
% In this function you implement your estimator. The function arguments
% are:
%  time: current time in [s]
%  dt: current time step [s]
%  internalStateIn: the estimator internal state, definition up to you.
%  steeringAngle: the steering angle of the bike, gamma, [rad]
%  pedalSpeed: the rotational speed of the pedal, omega, [rad/s]
%  measurement: the position measurement valid at the current time step
%
% Note: the measurement is a 2D vector, of x-y position measurement.
%  The measurement sensor may fail to return data, in which case the
%  measurement is given as NaN (not a number).
%
% The function has four outputs:
%  est_x: your current best estimate for the bicycle's x-position
%  est_y: your current best estimate for the bicycle's y-position
%  est_theta: your current best estimate for the bicycle's rotation theta
%  internalState: the estimator's internal state, in a format that can...
% be understood by the next call to this function

% Assigning variables from previous time step
x = internalStateIn.x;
y = internalStateIn.y;
theta = internalStateIn.theta;

%%%%%%%%%%%%%%%% EKF %%%%%%%%%%%%%%%%%%%%%%

% Variables for modelling
r = 0.425;
B = 0.8;
v = 5*r*pedalSpeed;

%%%%%%%%%%%%%%%% Noise %%%%%%%%%%%%%%%%%%%%

% xm(0)
Pm = internalStateIn.P; %This value is intialized in estInitialize

% Process Noise [E[v] = 0]
varR = 0.02125^2;
varB = 0.07968^2;
V = dt*[varR, 0, 0;
    0, varR, 0;
    0, 0, varR + varB];


% Measurement Noise [E[w] = 0]
W =[1.0893, 1.5333;
    1.5333, 2.9880];

%%%%%%%%%%%%% Prediction Update %%%%%%%%%%%

% Prediction Model is not needed as A and L are calculated by hand.
% Linearized Matrices
A =[1, 0, -(v*sin(theta))*dt;
    0, 1, (v*cos(theta))*dt;
    0, 0, 1];

L = eye(3);

% Prediction Step
xp = x + (v*cos(theta))*dt;
yp = y + (v*sin(theta))*dt;
thetap = theta + (v*tan(steeringAngle))*dt/B;
Xp = [xp; yp; thetap];

Pp = A*Pm*A' + L*V*L';

%%%%%%%%%%%%% Measurement Update %%%%%%%%%%%

% Measurement Model
h =[xp + 0.5*B*cos(thetap);
    yp + 0.5*B*sin(thetap)];

% Linearized Matrices
H =[1, 0, -0.5*B*sin(thetap);
    0, 1, 0.5*B*cos(thetap)];

M = eye(2);

% Measurement Step
K = Pp*H'*inv(H*Pp*H' + M*W*M'); %Gain

if ~isnan(measurement)
    xm = Xp + K*(measurement' - h);
    %         Pm = (eye(3) - K*H)*Pp; %This method is less accurate
    Pm = (eye(3) - K*H)*Pp*(eye(3) - K*H)' + K*W*K';
else
    xm = Xp;
    Pm = Pp;
end

% Storing new values before returning the function
x = xm(1);
y = xm(2);
theta = xm(3);


%% OUTPUTS %%
% Update the internal state (will be passed as an argument to the function
% at next run), must obviously be compatible with the format of
% internalStateIn:

internalStateOut.x = x;
internalStateOut.y = y;
internalStateOut.theta = theta;
internalStateOut.P = Pm;

end
