function [A] = calcF(crystal, mu, S, r, k0z)
% calcF()
% Calculate the kinematically diffracted amplitude based on scattering
% factors, reciprocal space vectors, and atomic positions.
% Assume z = 0 corresponds to a "mean plane" surface, with atomic positions 
% z < 0 below the surface.
% Inputs:
% f         Atomic scattering factors. Vector of length m
% mu        Effective absorption coefficient, in Angstroms^-1
% S         Reciprocal space vectors. 3D array
% r         Atomic positions. Array with length m and height 3
% k0z       z-component (transverse) of incident wavevector, in Ang^-1
%
% Outputs:
% A         Electron diffraction amplitude. 3D array

% Initialize loop variables.
[m,n] = size(squeeze(S(:,:,1)));
A = zeros(m, n);
M = length(crystal.atoms);
temp = zeros(1, M);

for i=1:m
    for j=1:n
        % Get the (Sx, Sy, Sz) triplet.
        s = squeeze(S(i,j,:));

        % If Sz is less than -k0z, the beam is totally internally
        % reflected.
        if (s(3) <= -k0z)
            temp = 0*temp;

        % Otherwise, calculate the diffracted amplitude.
        else
            for k=1:M
                % Assign the atomic scattering factors.
                atom = crystal.atoms(k);
                f = crystal.f.(atom)(i,j);
                % DEBUG
                % f = crystal.f.(atom)(46,101);
                % DEBUG
                % f = 1;

                % Calculate the diffracted amplitude.
                temp(k) = f*exp(1i*(s(1)*r(k,1) + s(2)*r(k,2) + s(3)*r(k,3)))*exp(mu*r(k,3)/2);
            end
        end

        A(i,j) = sum(temp);
    end
end

end