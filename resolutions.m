% Define physical constants.
hbar = 1.0546e-34; % J*s
m0 = 9.1094e-31;    % kg
c = 3e8;            % m/s
eV = 1.6022e-19;    % J = kg*m^2/s^2

% Specify electron beam energy and FWHM.
E0 = 20e3;          % eV
delE = 2;           % eV

% Specify incident beam angle and RHEED screen distance.
theta = 3;          % degrees
phi = 0;            % degrees
d = 30;             % cm

% Specify the RHEED screen coordinates (xd, yd) and corresponding pixel
% sizes delxd and delyd.
xd = 0;                         % cm
yd = d*tan(theta*pi/180);       % cm
delxd = 10e-4;                  % cm
delyd = 10e-4;                  % cm

% Calculate the derivative dK0/dE.
dK0dE = (1/hbar)*(m0 + (E0*eV)/(c^2))/sqrt(2*m0*E0*eV + (E0*eV/c)^2);   % m/J

% Calculate the reciprocal space resolution corresponding to primary beam 
% FWHM = delE.
K0 = calck0(E0);                    % Angstroms^-1
delK0 = dK0dE*delE*eV*1e-10;        % Angstroms^-1

% Calculate the scattering vector derivatives with respect to the screen
% position coordinates. Assume small scattering vector magnitude S and
% small inner potential such that the magnitude of the scattered wavevector
% K and incident wavevector K0 are approximately equal.
dSzdyd = K0*(d^2)*(d^2 + yd^2)^(-3/2);          % Angstroms^-1 cm^-1

% Partial derivatives w.r.t. xd. Assume that Sx ~= K ~= K0 and that Sy ~=
% 0.
Sx = K0; Sy = 0;                                % Angstroms^-1
dSydxd = (1/d)*(Sx + K0*cos(phi*pi/180));       % Angstroms^-1 cm^-1
dSxdxd = -d/(xd^2)*(Sy + K0*sin(phi*pi/180));   % Angstroms^-1 cm^-1

% Calculate the reciprocal space resolution corresponding to pixel size
% delxd, delyd.
delSz = dSzdyd*delyd;               % Angstroms^-1
delSx = dSxdxd*delxd;               % Angstroms^-1
delSy = dSydxd*delxd;               % Angstroms^-1

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assume xd ~= 0, yd = d*tan(theta), Sy is small, neglect inner potential.
% Assume azimuth phi = 0.

% Define physical constants.
hbar = 1.0546e-34; % J*s
m0 = 9.1094e-31;    % kg
c = 3e8;            % m/s
eV = 1.6022e-19;    % J = kg*m^2/s^2

% Specify electron beam energy and FWHM.
E0 = 20e3;          % eV
delE = 2;           % eV

% Specify incident beam angle and RHEED screen distance.
theta = 3;              % degrees
theta = theta*pi/180;   % radians
phi = 0;                % degrees
d = 30;                 % cm

% Specify the RHEED screen coordinates (xd, yd) and corresponding pixel
% sizes delxd and delyd.
xd = 0;                         % cm
yd = d*tan(theta*pi/180);       % cm
% delxd = 10e-4;                  % cm
% delyd = 10e-4;                  % cm

% Calculate the derivative dK0/dE.
dK0dE = (1/hbar)*(m0 + (E0*eV)/(c^2))/sqrt(2*m0*E0*eV + (E0*eV/c)^2);   % m^-1*J^-1

% Calculate the reciprocal space resolution corresponding to primary beam 
% FWHM = delE.
K0 = calck0(E0);                    % Angstroms^-1
delK0 = dK0dE*delE*eV*1e-10;        % Angstroms^-1

% Calculate the scattering vector components. Assume Sy is small.
Sz = 2*K0*sin(theta);                                       % Angstroms^-1
Sx = sqrt(K0^2 - (Sz - K0*sin(theta))^2) - K0*cos(theta);   % Angstroms^-1
Sy = 0;                                                     % Angstroms^-1

% Calculate the partial derivative of scattering vector w.r.t.
% incident beam wavevector.
dSxdK0 = (K0*(1 - sin(theta)^2) + Sz*sin(theta))/sqrt(K0^2 - (Sz - K0*sin(theta))^2) - cos(theta);     % Dimensionless
dSzdK0 = (K0*(1 - cos(theta)^2) - Sx*cos(theta))/sqrt(K0^2 - (Sx + K0*cos(theta))^2) + sin(theta);     % Dimensionless

% Calculate the partial derivative of RHEED screen coordinate xd w.r.t. 
% scattering vector x-component Sx.
dxddSx = -d*Sy/(K0^2)*(1 + Sx/K0)^(-2);                     % cm*Angstroms
dyddSz = (d/K0)*(1 - (Sz/K0 - sin(theta))^2)^(-3/2);        % cm*Angstroms

% Approximate the RHEED screen resolution by the product of partial
% derivatives.
delxd = dxddSx*dSxdK0*dK0dE*(1e-10)*(delE*eV);              % cm
delyd = dyddSz*dSzdK0*dK0dE*(1e-10)*(delE*eV);              % cm
delSx = dSxdK0*delK0;                                       % Angstroms^-1
delSz = dSzdK0*delK0;                                       % Angstroms^-1