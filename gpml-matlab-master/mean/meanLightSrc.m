function [m,dm] = meanLightSrc(hyp, x)

% disk light source model as the mean function. The output is the radiant
% intensity field.
% Constant mean function. The mean function is parameterized as:

% x = [radius, theta]
%
% m(x) = (phi_0*cos(theta)^mu)/(r/r_d + 1)^2
%
% The hyperparameter is:
%
% hyp = [ phi0, mu, r_d ]


%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2016-04-15.
%
% See also meanFunctions.m.

if nargin<2, m = '3'; return; end             % report number of hyperparameters 
if numel(hyp)~=3, error('Exactly three hyperparameter needed.'), end

phi_0 = hyp(1);
mu = hyp(2);
r_d = hyp(3);

r = x(:,1);
theta = x(:,2);

%mean
% m = (phi_0.*(cos(theta).^mu))./(r/r_d + 1)^2;
m = (phi_0.*(r_d.^2).*(cos(theta).^mu))./((r + r_d).^2);

% q = normalize(ones(length(m),1), 'norm');
% % directional derivative
% dm = [(r_d.^2*cos(theta).^mu)./(r + r_d).^2, ... 
%     (phi_0.*r_d.^2.*log(cos(theta)).*cos(theta).^mu)./(r + r_d).^2, ...
%     (2.*phi_0.*r.*r_d.*cos(theta).^mu)./(r + r_d).^3];
% 
% dm1 = q'*dm;

% directional derivative
dm = @(q) q'*[(r_d.^2*cos(theta).^mu)./(r + r_d).^2, ... 
    (phi_0.*r_d.^2.*log(abs(cos(theta))).*cos(theta).^mu)./(r + r_d).^2, ...
    (2.*phi_0.*r.*r_d.*cos(theta).^mu)./(r + r_d).^3];
