function [ G, Gmag ] = GetG( lattice, hkl )
% GetG
% Calculate the reciprocal lattice vectors G for each h,k,l triplet.
%
% Inputs:
% lattice       3x3 array with the unit cell lattice vectors, in Cartesian
%               coordinates. In Angstroms
% hkl           Nx3 array of N different h,k,l triplets
%
% Outputs:
% G             Nx3 array of reciprocal lattice vectors in Cartesian
%               coordinates. In Angstroms^-1
% Gmag          Length N vector of magnitude of reciprocal lattice vectors.

% Calculate the reciprocal lattice vectors.
a = lattice(1,:);                   % Angstroms
b = lattice(2,:);                   % Angstroms
c = lattice(3,:);                   % Angstroms
omega = dot(a, cross(b, c));        % Angstroms^3

a1 = 2*pi/omega*cross(b,c);         % Angstroms^-1
b1 = 2*pi/omega*cross(c,a);         % Angstroms^-1
c1 = 2*pi/omega*cross(a,b);         % Angstroms^-1

% Make sure the h,k,l triplets are in the shape Nx3.
% sz = size(hkl);
% 
% if sz(2) ~= 3
%     hkl = transpose(hkl);
% end

% Generate the array of reciprocal lattice vectors Ghkl.
sz = size(hkl);
N = sz(1);
Ghkl = zeros(N, 3, 3);

for n=1:N
    Ghkl(n,1,:) = a1*hkl(n,1);
    Ghkl(n,2,:) = b1*hkl(n,2);
    Ghkl(n,3,:) = c1*hkl(n,3);
end

% Flatten Ghkl into an Nx3 array of Cartesian coordinates. Calculate the
% magnitude of G.
G = zeros(N, 3);
Gmag = zeros(N, 1);

for n=1:N
    G(n,:) = squeeze(sum(Ghkl(n,:,:), 2));
    Gmag(n) = sqrt(dot(G(n,:), G(n,:)));
end

end