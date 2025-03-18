function [ x, y ] = ProjScreen( k, K0, Vi, thetai, psii, d )
% ProjScreen
% Project the internal scattered wavevectors k onto the RHEED screen.
% Account for the inner potential Vi.
%
% Inputs:
% k             Nx3 array of internal scattered wavevectors in Cartesian
%               coordinates. In Angstroms^-1
% K0            Incident electron beam wavevector, in Angstroms^-1
% Vi            Inner potential, in eV
% thetaii       Incident beam angle, in radians
% psii          Incident beam azimuth, in radians
% d             Distance from sample to RHEED screen, in cm
%
% Outputs:
% x, y          Screen coordinates, in cm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define some physical constants.
hbar = 1.0546e-34; % J*s
m0 = 9.1094e-31;    % kg
eV = 1.6022e-19;    % J = kg*m^2/s^2

% Calculate the x, y, and z components of the incident electron beam.
K0x = K0*cos(thetai)*cos(psii);
K0y = K0*cos(thetai)*sin(psii);
K0z = -K0*sin(thetai);

% Calculate the external (vacuum) scattering vector S. The perpendicular
% component is modified by the inner potential Vi.
Sx = k(:,1) - K0x;
Sy = k(:,2) - K0y;
Kz = sqrt(k(:,3).^2 - 2*m0*eV*Vi/(hbar^2)*(1e-20));
Sz = Kz - K0z;

% Convert the scattering vectors to scattering angles.
thetaf = asin(Sz/K0 - sin(thetai));
psif = atan((Sy/K0 + sin(psii))./(Sx/K0 + cos(psii)));

% Project the scattering angles onto the RHEED screen.
x = d*tan(psif);
y = d*tan(thetaf);

end