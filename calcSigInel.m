function [SigInel] = calcSigInel(atom, elements, K0, delE, E)
% calcSigInel
% Calculate the inelastic scattering cross-section of the atom using the
% Compton incoherent scattering functions.
%
% Inputs:
% atom      String containing the atom name
% elements  Library of elemental parameters
% K0        Incident electron beam wavevector, in Angstroms^-1
% delE      Energy loss between incident and scattered electrons, in eV
% E         Incident electron beam energy, in eV
%
% Outputs:
% SigInel   Inelastic scattering cross-section of the atom, in Angstroms^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fetch the incoherent scattering function from the database of elements.
element = elements.(atom);
Incoherent = element.Incoherent;

% Define some physical constants.
hbar = 1.0546e-34;  % J*s
m0 = 9.1094e-31;    % kg
eV = 1.6022e-19;    % J = kg*m^2/s^2
aH = 5.29177e-11;   % m

% Calculate integration limits Smin and Kn + K0.
Smin = (K0/2)*(delE/E);                         % Angstroms^-1
Kn = sqrt(K0^2 - 2*m0*delE*eV/(hbar^2)*1e-20);  % Angstroms^-1

% Set up a 1-D vector of scattering parameters to integrate over.
N = 101;
Sq = logspace(log10(Smin), log10(K0 + Kn), N);            % Angstroms^-1

% Interpolate this vector between the incoherent scattering function lookup
% table values.
Incq = interp1(Incoherent.S, Incoherent.Inc, Sq);   % Dimensionless

% Calculate the differential inelastic scattering cross-section.
dSigInel = ((2/(aH*1e10))^2)*Incq./(Sq.^4);         % Angstroms^2

% Integrate the differential inelastic scattering cross-section to obtain
% the total inelastic scattering cross-section.
SigInel = trapz(Sq, Sq.*dSigInel);          % dimensionless
SigInel = (2*3.14159/(K0^2))*SigInel;       % Angstroms^2
end