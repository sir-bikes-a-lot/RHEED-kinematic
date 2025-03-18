addpath('cif');
addpath('gif');
addpath('gif/gif');
addpath('vasplab');
addpath('Crystal models');

psi = 60;           % degrees
% hkl = [5,5,2];
hkl = [5,5,1];
T = 300;            % Kelvin
theta = 3;          % degrees
d = 30;             % cm
radius = 9.3/2;     % cm
E0 = 20e3;          % eV
% filename = 'GaN_0001_2x2_Ga-vacancy.vasp';
filename = 'GaN_0001_2x2_N-H3.vasp';

[r, xd, yd, S, I, crystal, xk1, yk1, xk2, yk2] = CalcRHEED(filename, pi/180*theta, pi/180*psi, T, radius, d, hkl, E0);
[ L0, L ] = PlotRHEED(xd, yd, I, d, theta, 1e10, xk1, yk1, xk2, yk2);
print(gcf, 'GaN_0001_2x2_N-H3_2-1-10', '-dpng','-r600'); 

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('cif');
addpath('gif');
addpath('gif/gif');
addpath('vasplab');
addpath('Crystal models');

% psi = 60;           % degrees
psi = linspace(0, pi/3, 11);
hkl = [5,5,1];
T = 300;            % Kelvin
theta = 3;          % degrees
% mu0 = [1e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1];         % Angstroms^-1
E0 = 20e3;          % eV
d = 40;             % cm
radius = 9.3/2;     % cm
filename = 'GaN_0001_2x2_N-H3.vasp';

for i=1:length(psi)
    [r, xd, yd, S, I, crystal, xk1, yk1, xk2, yk2] = CalcRHEED(filename, theta*pi/180, psi(i), T, radius, d, hkl, E0);
    [ L0, L ] = PlotRHEED(xd, yd, I, d, theta, 1e10, xk1, yk1, xk2, yk2);
    % [r, xd, yd, S, I, crystal] = CalcRHEED(filename, theta*pi/180, psi, T, radius, d, hkl, mu0(i));
    % [ L0, L ] = PlotRHEED(xd, yd, I, d, theta*pi/180, 1e5);
    if(i==1)
        gif('GaN_0001_2x2_N-H3_vary_azimuth.gif', 'DelayTime', 0.5);
    else
        gif('frame',gcf);
    end
end