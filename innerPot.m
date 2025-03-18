function [ s, k0z ] = innerPot( K0, theta, psi, S, Vi )
% innerPot
% Calculate the internal scattering vectors after applying boundary
% conditions to the incident electron beam.
%
% Inputs:
% K0        Electron beam wavevector in vacuum, in Ang^-1
% theta     Incident beam angle, in radians
% psi       Incident beam azimuth, in radians
% S         Array of vacuum scattering vector triplets, in Ang^-1
% Vi        Inner potential in crystal, in eV
%
% Outputs:
% s         Array of internal scattering vector triplets, in Ang^-1
% k0z       z-component (transverse) of incident wavevector, in Ang^-1

% Define some physical constants.
hbar = 1.0546e-34; % J*s
m0 = 9.1094e-31;    % kg
% c = 3e8;            % m/s
eV = 1.6022e-19;    % J = kg*m^2/s^2

% Unpack the vacuum scattering vector.
Sx = squeeze(S(:,:,1));
Sy = squeeze(S(:,:,2));
Sz = squeeze(S(:,:,3));

% Calculate the x, y, and z-components of the vacuum wavevector.
% K0x = K0*cos(theta)*cos(psi);   % Angstroms^-1
% K0y = K0*cos(theta)*sin(psi);   % Angstroms^-1
K0z = -K0*sin(theta);           % Angstroms^-1

% Calculate the z-component of the scattered electron wavevector in vacuum.
Kz = Sz + K0z;                  % Angstroms^-1

% Calculate the z-component of the scattered wavevector in the crystal.
kz = sqrt((Kz*1e10).^2 + 2*m0*eV*Vi/(hbar^2))*(1e-10);       % Angstroms^-1
k0z = sqrt((K0z*1e10).^2 + 2*m0*eV*Vi/(hbar^2))*(1e-10);     % Angstroms^-1

% Calculate the z-component of the internal scattering vector.
sz = kz - k0z;                  % Angstroms^-1

% Pack up the internal scattering vector. The transverse x- and
% y-components are unchanged from vacuum.
s = zeros(size(S));
s(:,:,1) = Sx;
s(:,:,2) = Sy;
s(:,:,3) = sz;

% Pack up the internal scattered wavevector. The transverse x- and
% y-components are unchanged from vacuum.
% k = zeros(size(S));
% k(:,:,1) = Sx + K0x;
% k(:,:,2) = Sy + K0y;
% k(:,:,3) = kz;

end