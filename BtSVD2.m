function [W,alpha,isigma2,SSx]=BtSVD2(U,Y,s)
% [W,alpha,isigma2]=BtSVD_ml(U,y,s)
%
% input:
% U = principal components of data "X" matrix, including the intercept column
% y = observations
% s = singular values of the original X data,
%
% returns:
% W = fitted regression parameters
% alpha = fitted alpha parameters [a0,a]
% isigma2 = fitted 1/sigma2 observation error
%
% Model: y = U*W + eps, eps ~ N(0,I*sigma2), W ~ N(0,1/alpha), alpha = [a0,a,a,...]'
%
% References
% Tipping 2004 page 17, section 4.3.3
% Junttila 2015 

% This version includes intercept and assumes size(y,2) == 1

% add intercept and call it X (need to decide this)
if size(U,2) == length(s)
  X = [ones(size(U,1),1),U];
else
  X = U;
end

[N,M] = size(X);

% initial values
a0      = 0.1;
a       = 0.1;
isigma2 = 1.0;
Imax    = 50; % maximum number of iterations
Imin    = 4;  % minimum number of iterations
Amax    = 1e8; % max alpha
S0      = 1;

a0_prev = a0;
a_prev  = a;
A = [a0;repmat(a,M-1,1)] ;
SSi = [1/S0;1./s(:).^2];

XX = X'*X;
XY = X'*Y;

for ii=1:Imax
  % Calculate inverses by Cholesky decomposition
  SSxich = chol(isigma2*XX+diag(A.*SSi));
  W = SSxich\(SSxich'\XY)*isigma2;             % W = (X'X + A)\X'y
  SSx = SSxich\(SSxich'\eye(size(SSxich)));    % SSx = inv(X'X + A)
  
  gamma = 1-A.*SSi.*diag(SSx);
  a0 = gamma(1)/W(1).^2;
  a = sum(gamma(2:end))./sum(SSi(2:end).*W(2:end).^2);
  a0 = min(Amax,a0);
  a = min(Amax,a);
  A  = [a0;repmat(a,M-1,1)];
  isigma2 = (N-sum(gamma))/sum((Y-X*W).^2);

  if ii>=Imin && (a>=Amax || (abs(a0-a0_prev)+abs(a_prev-a))<1e-4/N)
    break
  else 
    a0_prev = a0;
    a_prev  = a; 
  end
end

alpha = [a0,a];
