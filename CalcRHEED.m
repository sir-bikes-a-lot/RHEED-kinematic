function [r, xd, yd, S, I, crystal, xk1, yk1, xk2, yk2] = CalcRHEED(filename, theta, psi, T, radius, d, hkl, E0)
% Development function to calculate RHEED patterns at different
% azimuths.
addpath('Crystal models');
addpath('Libraries');
addpath('vasplab');

% Import parameters for each atom from library.
Ga = elements('Ga', 'GaN');
N = elements('N', 'GaN');
ElementsLib = struct('Ga', Ga, 'N', N);

% Calculate electron beam wavevector.
K0 = calck0(E0);      % Ang^-1
lambda = 2*pi./K0;   % Angstroms

% Specify number of mesh points to use for RHEED screen coordinates.
Nx = 201;
Ny = 201;

% Set up RHEED screen coordinate mesh and corresponding reciprocal space
% mesh.
[xd, yd, S] = calcSmesh(radius, d, theta, 0, K0, Nx, Ny);

% DEBUG: Set up a 1D line.
% Sz = linspace(0,10,1001);
% Sx = 0*Sz + 0*pi/5.189;
% Sy = 0*Sz + 0*pi/5.189;
% S = transpose([Sx; Sy; Sz]);

% TO-DO: Implement a way to calculate the inner potential, Vi.
Vi = 10;                                        % eV

% Convert the vacuum scattering vectors S to the internal scattering 
% vectors s.
[s, k0z] = innerPot(K0, theta, psi, S, Vi);     % Angstroms^-1

% Calculate magnitude of the internal scattering vectors.
smag = sqrt(s(:,:,1).^2 + s(:,:,2).^2 + s(:,:,3).^2)/(4*pi);    % Angstroms^-1

% Import crystal model.
crystal = LoadCrystal(filename, hkl, ElementsLib, smag, T);

% Loop through all the atoms in the list.
M = length(crystal.atoms);
r = zeros(M, 3);

for m=1:M
    % Rotate crystal by azimuthal angle psi.
    r(m,:) = basistrans(transpose(crystal.r(m,:)), pi/2, psi);       % Angstroms
end

% Calculate the mean absorption coefficient. Assume a fixed inelastic
% scattering energy loss, delE.
delE = 15;                                          % eV
mu0 = calcmu(crystal, ElementsLib, K0, delE, E0);    % Ang^-1
mu = mu0./sin(theta);                                   % Ang^-1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEBUG: Turn off Kikuchi line calculation.

% Calculate the Kikuchi line locations.
[kk1, kk2] = Kikuchi(crystal.UClattice, K0, Vi, 2);

% Project the Kikuchi line wavevectors onto the RHEED screen.
sz = size(kk1);
M = sz(1);
xk1 = zeros(M, sz(2)); yk1 = xk1;
xk2 = zeros(M, sz(2)); yk2 = xk2;

for m=1:M
    [xk1(m,:), yk1(m,:)] = ProjScreen(squeeze(kk1(m,:,:)), K0, Vi, theta, psi, d);
    [xk2(m,:), yk2(m,:)] = ProjScreen(squeeze(kk2(m,:,:)), K0, Vi, theta, psi, d);
end

% xk1 = 0; yk1 = 0; xk2 = 0; yk2 = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DEBUG
% figure;
% plot(xk1(1,:), yk1(1,:), 'k-'); hold on;
% plot(xk2(1,:), yk2(1,:), 'k-');
% for m=2:M
%     plot(xk1(m,:), yk1(m,:), 'k-');
%     plot(xk2(m,:), yk2(m,:), 'k-');
% end
% hold off;
% axis([min(min(xd)), max(max(xd)), min(min(yd)), max(max(yd))]);

% Calculate diffracted intensity.
A = calcF(crystal, mu, s, r, k0z);      % Dimensionless
I = A.*conj(A);                         % Dimensionless

end