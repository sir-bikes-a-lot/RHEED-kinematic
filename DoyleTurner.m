function [ f ] = DoyleTurner( DT, B, s )
% DoyleTurner
% Calculate the complex electron scattering factor for scattering angle
% theta and electron wavelength lambda. Calculate the Debye-Waller
% correction factor in another function.
%
% Inputs:
% DoyleTurner   Structure containing the Doyle-Turner parameters for the
%               complex electron scattering factor.
% B             Debye-Waller factor, in Ang^2
% s             Array of scattering vectors, in Angstroms^-1
%
% Outputs:
% f             Complex electron scattering factor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unpack the Doyle-Turner parameters.
aRe = [DT.a1, DT.a2, DT.a3, DT.a4, DT.a5];                 % Angstroms
bRe = [DT.b1, DT.b2, DT.b3, DT.b4, DT.b5];                 % Angstroms^2
aIm = [DT.aTDS1, DT.aTDS2, DT.aTDS3, DT.aTDS4, DT.aTDS5];  % Angstroms
bIm = [DT.bTDS1, DT.bTDS2, DT.bTDS3, DT.bTDS4, DT.bTDS5];  % Angstroms^2

% DEBUG
% E = 20e3;           % eV
% K0 = calck0(E);      % Ang^-1
% lambda = 2*pi/K0;    % Angstroms
% s = sin(3*pi/180)/lambda;

% Calculate the complex electron scattering factor.
f = zeros([length(aRe), size(s)]);

for i=1:length(aRe)
    f(i,:,:) = aRe(i)*exp(-bRe(i)*s.^2) + 1i*aIm(i)*exp(-(bIm(i) - B/2)*s.^2);
    % f(i,:,:) = aRe(i)*exp(-bRe(i)*s.^2);
end

f = squeeze(sum(f,1)).*exp(-B*s.^2);

end