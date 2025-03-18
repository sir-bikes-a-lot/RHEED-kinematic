function [ L0, L ] = PlotRHEED(xd, yd, I, d, theta, thresh, xk1, yk1, xk2, yk2)
% Development function to plot RHEED patterns.

temp = size(xd);
Nx = temp(1); Ny = temp(2);

% Define a grid of reciprocal space vectors to generate the RHEED pattern.
xlim = 10;
ylim = 10;

Nl = 101;
x = linspace(-xlim, xlim, Nl);             % Ang^-1
y = linspace(-ylim, ylim, Nl);             % Ang^-1

% Specify the broadening of the reciprocal lattice points.
% sigx = 0.1;                             % Ang^-1
% sigy = 0.1;                             % Ang^-1
sigx = 0.001;                             % Ang^-1
sigy = 0.001;                             % Ang^-1

% Create the 2D Gaussian convolution function.
H = zeros(Nl, Nl);

for i=1:Nl
    for j=1:Nl
        H(i,j) = exp(-((x(i))^2)/(2*sigx^2) -((y(j))^2)/(2*sigy^2));
    end
end

% Broaden the reciprocal lattice points using convolution (faster).
L0 = conv2(I,H,'same');

% Apply a max threshold to the intensity map.
L = L0;

for m=1:Nx
    for n=1:Ny
        if (L0(m,n) > thresh)
            L(m,n) = thresh;
        end
    end
end

% Plot the reciprocal lattice diffraction intensity as a slice in the (x,y)
% plane.
% figure;
% contourf(S(:,:,1), S(:,:,2), I, 'LineStyle', 'none');
% % axis([-K,K,-K,K]);
% Ang = char(197);
% xlabel(strcat('S_{x} (', Ang, '^{-1})')); ylabel(strcat('S_{y} (', Ang, '^{-1})'));
% title('Diffraction intensity, I');
% colormap gray;
% savestr = strcat('Ihk_map_psi_', num2str(psi, 3), '.png');
% saveas(gcf, savestr);

% Plot the broadened 2D intensity map of the reciprocal lattice.
% figure;
% contourf(S(:,:,1), S(:,:,2), L, 'LineStyle', 'none');
% % axis([-K,K,-K,K]);
% Ang = char(197);
% xlabel(strcat('k_{x} (', Ang, '^{-1})')); ylabel(strcat('k_{y} (', Ang, '^{-1})'));
% title('Intensity map, L');
% colormap gray;
% savestr = strcat('Lxy_map_psi_', num2str(psi, 3), '.png');
% saveas(gcf, savestr);

% Plot RHEED pattern.
figure;
contourf(xd, yd, L, 'LineStyle', 'none'); hold on;

% Plot the specular reflection.
scatter(0, d*tan(theta), 40, 'MarkerEdgeColor',[0 0.5 0.5],...
              'MarkerFaceColor',[0 0.7 0.7],...
              'LineWidth',1.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEBUG: Turn off Kikuchi line plotting.

% Plot the Kikuchi lines.
sz = size(xk1);
M = sz(1);
plot(xk1(1,:), yk1(1,:), 'w-'); hold on;
plot(xk2(1,:), yk2(1,:), 'w-');
for m=2:M
    plot(xk1(m,:), yk1(m,:), 'w-');
    plot(xk2(m,:), yk2(m,:), 'w-');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hold off;
xlabel('x_{d} (cm)'); ylabel('y_{d} (cm)');
title('RHEED pattern');
colormap gray;
% clim([0, max(max(max(L)), thresh)]);
axis([min(min(xd)), max(max(xd)), min(min(yd)), max(max(yd))]);
axis square
% savestr = strcat('RHEED_pattern_psi_', num2str(psi, 3), '_rad_h', num2str(hkl(1)), '_k', num2str(hkl(2)), '_l', num2str(hkl(3)), '.png');
% saveas(gcf, savestr);

end