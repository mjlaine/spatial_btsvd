function [cmat,L] = makecovmat2(d,tau2,sig2,phi, covmodel)
% generate covariance matrix given distance matrix and
% covariance function parameters
% [cmat,L] = fun(d,tau2,sig2,phi)

if nargin < 5
  covmodel = 'gaussian';
end

n = size(d,1);

cmat = zeros(n,n);

switch covmodel
  case 'gaussian'
    % Gaussian covariance function
    cfun = @(d) sig2.*exp(-0.5*d.^2/phi.^2) ;
  case 'exponential'
    cfun = @(d) sig2.*exp(-abs(d)/phi) ;
  otherwise
    error('unknown covariance model');
end

%cmat = cfun(d,sig2,phi) + eye(n)*tau2;
cmat = cfun(d) + eye(n)*tau2;

if nargout>1
  L = chol(cmat,'lower');
end
