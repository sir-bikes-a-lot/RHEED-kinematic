function [k0] = calck0(E)
% calck0
% Calculate electron beam wavevector k0.
% Inputs:
% E         Electron beam energy, in eV
%
% Ouputs:
% k0        Electron beam wavevector, in Ang^-1

hbar = 1.0546e-34; % J*s
m0 = 9.1094e-31;    % kg
c = 3e8;            % m/s
eV = 1.6022e-19;    % J = kg*m^2/s^2
k0 = 1/hbar*sqrt(2*m0*E*eV + ((E*eV).^2)/(c^2))*1e-10;       % Ang^-1

end