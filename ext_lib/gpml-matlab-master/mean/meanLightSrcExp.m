function [m,dm] = meanLightSrcExp(hyp, x)
% Radiant intensity light source disk model as the mean function. The
% hyperparameters of the mean function are the exponents of the exponential
% function which evalulate to the unknown parameters of the light model.
% This parameterisation ensures that the resulting parameters are always
% positive.
% INPUTS:
%       x - [radius, theta] measurements
%       hyp - [h1, h2, h3] hyperparameters
% OUTPUTS:
%       m - radiant intensity from mean function
%       dm - jacobian of m w.r.t hyp
%
% hyp = [ phi0, mu, r_d ]
%
% phi0 = exp(h1)
% mu = exp(h2)
% r_d = exp(h3)

if nargin<2, m = '3'; return; end             % report number of hyperparameters
if numel(hyp)~=3, error('Exactly three hyperparameter needed.'), end

h1 = hyp(1);
h2 = hyp(2);
h3 = hyp(3);

r = x(:,1);
theta = x(:,2);

%mean
m = (exp(h1).*cos(theta).^exp(h2))./(r.*exp(-h3) + 1).^2;

% directional derivative
dm = @(q) q'*[(exp(h1).*cos(theta).^exp(h2))./(r.*exp((-h3)) + 1).^2, ...
    (log(cos(theta)).*exp((h1 + h2)).*cos(theta).^exp(h2))./(r.*exp((-h3)) + 1).^2, ...
    (2.*r.*exp((-h3)).*exp(h1).*cos(theta).^exp(h2))./(r.*exp((-h3)) + 1).^3];
end