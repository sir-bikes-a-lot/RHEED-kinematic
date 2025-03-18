function [ crystal ] = LoadCrystal( filename, hkl, elements, s, T )
% LoadCrystal
% Load a VASP POSCAR file containing the crystal information, including the
% atom IDs and coordinates. Uses the VASPLAB package: https://github.com/max-radin/VASPLAB
%
% Inputs:
% filename      String containing the .vasp file name.
% hkl           Triplet with supercell size, integers >= 1.
% elements      Library of elemental parameters
% s             Array of scattering vectors, in Angstroms^-1
% T             Temperature, in Kelvin
%
% Outputs:
% crystal       Structure containing the atom names and Cartesian
%               coordinates (x,y,z), the atomic scattering factors f,
%               and the lattice.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add the VASPLAB filepath.
addpath('vasplab');

% Load the .vasp file.
geometry = import_poscar(filename);

% Generate a supercell of size (h,k,l).
geometryExt = supercell(geometry, hkl);

% Unpack the lattice basis vectors, atom names and counts, and fractional
% coordinates.
lattice = geometryExt.lattice;          % Lattice basis vectors in Cartesian x,y,z coordinates, in Angstroms
atomnames = geometryExt.symbols;

% Generate a list of the atom names.
N = length(geometryExt.coords);
atoms = strings(N,1);
ind = 0;

for m=1:length(geometryExt.atomcount)
    for n=1:geometryExt.atomcount(m)
        atoms(n + ind) = cell2mat(geometryExt.symbols(m));
    end
    ind = geometryExt.atomcount(m);
end

% Calculate the Cartesian coordinates of each atom in the supercell.
r = geometryExt.coords*geometryExt.lattice;

% Offset all the z-coordinates such that z=0 corresponds to the highest
% fractional coordinate.
z0 = max(r(:,3));
r(:,3) = r(:,3) - z0;

% Calculate the total supercell volume.
omega = dot(lattice(1,:), cross(lattice(2,:), lattice(3,:)));

% Calculate atomic scattering factor at each reciprocal space mesh point
% for each element in the lattice.
for m=1:length(atomnames)
    % Calculate scattering factor at each scattering vector.
    fc = GetScatterFact(cell2mat(atomnames(m)), elements, s, T);

    % Pack this array into a structure.
    f.(cell2mat(atomnames(m))) = fc;
end

% Pack up the atom names and coordinates to the output structure.
crystal = struct('atoms', atoms, 'r', r, 'f', f, 'lattice', lattice, 'UClattice', geometry.lattice, 'omega', omega);

end