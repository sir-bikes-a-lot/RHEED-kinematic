function [W] = GaussFun(x, x0, sig)
% GaussFun
% Set up a Gaussian distribution W normalized such that sum(W) = 1.
% y-offset = 0.
%
% Inputs
% x         1-D vector of input values
% x0        x-intercept
% sig       Standard deviation
%
% Outputs
% W         1-D Gaussian function

W = exp(-((x - x0).^2)/(2*sig^2));
W = W/sum(W);

end