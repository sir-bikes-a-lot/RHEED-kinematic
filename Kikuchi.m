function [ kk1, kk2 ] = Kikuchi( lattice, K0, Vi, n)
% Kikuchi
% Solve for the location of the Kikuchi lines.
%
% Inputs:
% lattice       3x3 array with the unit cell lattice vectors, in Cartesian
%               coordinates. In Angstroms
% K0            Incident electron beam wavevector, in Angstroms^-1
% Vi            Inner potential, in eV
% n             Max index for h,k,l triplets.
%
% Outputs:
% kk            N x Nx x 3 array of scattered wavevector triplets
%               corresponding to the Kikuchi lines.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the array of h,k,l triplets up to n.
v = -n:n;                               % Integers between -n and n
p = [v, v, v];                          % 3 of each integer in this vector
hkl = unique(nchoosek(p, 3), 'rows');   % Permute the integer triplets and keep the unique rows
sz = size(hkl);                         % Get the size
N = sz(1);                              % This is the length of the h,k,l triplet vector

% Get the array of reciprocal lattice vectors.
[G, Gmag] = GetG(lattice, hkl);

% Define some physical constants.
hbar = 1.0546e-34; % J*s
m0 = 9.1094e-31;    % kg
eV = 1.6022e-19;    % J = kg*m^2/s^2

% The wavevector magnitude in the crystal is increased by the inner
% potential Vi.
k = sqrt(K0^2 + 2*m0*eV*Vi/(hbar^2)*(1e-20));       % Angstroms^-1

% Set up a vector of internal scattered wavevectors kx to solve for the
% Kikuchi lines.
kx = linspace(0, K0, 10001);          % Angstroms^-1

% Solve for the Kikuchi lines, which satisfy k.G = (|G|^2)/2.
kk1 = zeros(N, length(kx), 3);
kk2 = kk1;

for m=1:N
    % Define auxiliary variables a, b, and c to solve the quadratic
    % equation.
    a = G(m,2)^2 + G(m,3)^2;
    b = -2*G(m,2)*((Gmag(m)^2)/2 - kx.*G(m,1));
    c = ((Gmag(m)^2)/2 - kx*G(m,1)).^2 - (G(m,3)^2)*(k^2 - kx.^2);
    ky1 = (-b + sqrt(b.^2 - 4*a.*c))./(2*a);
    ky2 = (-b - sqrt(b.^2 - 4*a.*c))./(2*a);

    % Switch solution branch depending on the sign of b.
    % if b < 0
    %     ky = (-b - sqrt(b.^2 - 4*a.*c))./(2*a);
    % else
    %     ky = (-b + sqrt(b.^2 - 4*a.*c))./(2*a);
    % end

    % If the value of ky is imaginary, assign it the value 'NaN'.
    for i=1:length(kx)
        if ~isreal(ky1(i))
            ky1(i) = NaN;
        end
        if ~isreal(ky2(i))
            ky2(i) = NaN;
        end
    end

    % Solve for kz.
    kz1 = sqrt(k^2 - kx.^2 - ky1.^2);
    kz2 = sqrt(k^2 - kx.^2 - ky2.^2);

    % Pack up the scattered wavevector.
    kk1(m,:,1) = kx;
    kk1(m,:,2) = ky1;
    kk1(m,:,3) = kz1;

    kk2(m,:,1) = kx;
    kk2(m,:,2) = ky2;
    kk2(m,:,3) = kz2;
end

end