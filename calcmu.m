function [mu0] = calcmu(crystal, elements, K0, delE, E)
% calcmu
% Calculate the mean absorption coefficient of the crystal using the
% Compton incoherent scattering functions.
%
% Inputs:
% crystal   Structure containing the crystal supercell model parameters
% elements  Library of elemental parameters
% K0        Incident electron beam wavevector, in Angstroms^-1
% delE      Energy loss between incident and scattered electrons, in eV
% E         Incident electron beam energy, in eV
%
% Outputs:
% mu0       Mean absorption coefficient, in Angstroms^-1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the inelastic scattering cross-section for each atom in the crystal. 
% Do this by identifying unique instances of atom names and calculating the
% inelastic scattering cross-section only for those atoms.
a = {};
SigInelTot = 0;

for m=1:length(crystal.atoms)
    % Check the list to see if we have calculated the scattering factors
    % yet.
    if(~ismember(crystal.atoms(m), a))
        % Calculate inelastic scattering cross-section at each scattering 
        % vector.
        SigInelc = calcSigInel(crystal.atoms(m), elements, K0, delE, E);    % Angstroms^2

        % Pack this value into a structure.
        SigInel.(crystal.atoms(m)) = SigInelc;

        % Add the atom name to a list.
        if(isempty(a))
            a = crystal.atoms(m);
        else
            a = [a, crystal.atoms(m)];
        end
    end

    % Sum up the inelastic scattering cross-sections for all the atoms in
    % the crystal.
    SigInelTot = SigInelTot + SigInel.(crystal.atoms(m));   % Angstroms^2
end

% Calculate the mean absorption coefficient by dividing by the supercell
% volume.
mu0 = SigInelTot/crystal.omega;     % Angstroms^-1

end