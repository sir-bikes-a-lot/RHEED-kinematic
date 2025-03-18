function [ f,s ] = GetScatterFact( atom, elements, s, T )
% GetScatterFact
% Calculate complex electronic scattering factor for the specified atom.
%
% INPUTS
% Atom          String containing the atom name.
% elements      Library of elemental parameters
% s             Array of scattering vectors, in Angstroms^-1
% T             Temperature, in Kelvin
%
% Outputs:
% f             Complex electron scattering factor
% s             (Optional) scattering parameter, in Angstroms^-1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fetch the Doyle-Turner and static correlation function coefficients from
% the database of elements.
element = elements.(atom);
DT = element.DoyleTurner;
u11Params = element.u11Params;
u33Params = element.u33Params;

% Calculate the Debye-Waller factor B. Neglect the anisotropy and calculate
% the isotropic factor.
B = DebyeWaller(u11Params, u33Params, element.M, T);

% Calculate the complex electron scattering factor from the Doyle-Turner
% parameters.
f = DoyleTurner(DT, B, s);

end