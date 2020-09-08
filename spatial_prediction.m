function out=spatial_prediction(results,chain,xnew,xcoordnew,ycoordnew,nsimu)
% prediction from spatial mcmc results

if nargin<6
  nsimu=1;
end

dist = makedistmat([results.xcoord;xcoordnew],[results.ycoord;ycoordnew]);
nobs = size(dist,1);

ychain = zeros(nsimu,1);

bs = @(L,X)linsolve(L,X,struct('LT',true)); % 'backsolve', L\X, L lower diagonal

for isimu=1:nsimu

  if nsimu == 1
    m = mean(chain);
  else
    m = chain(ceil(rand*size(chain,1)),:);
  end

  alpha = m(results.alphaind)';
  beta  = m(results.betaind)';
  theta = m(results.thetaind)';

  if results.opts.useL == 0
    sig = m(end); % obs sig
    ynew = beta(1) + xnew*beta(2:end);
    ynew = ynew + randn(1,1)*sig;
  else

    if isempty(theta)
      tau2 = results.opts.tau2;
      sig2 = results.opts.sig2;
      phi = results.opts.phi;
    else
      tau2 = theta(1);
      sig2 = theta(2);
      phi = theta(3);
    end

    yhat = results.x*beta;
    % for the original data
    cmat = makecovmat2(dist(1:nobs-1,1:nobs-1),tau2,sig2,phi,results.opts.covmodel);
    L = chol(cmat,'lower');
    % tau2 = 0 here!
    cmat2 = makecovmat2(dist,0,sig2,phi,results.opts.covmodel);
    % here not (should we add nugget when in observation location?)
%    cmat2 = makecovmat2(dist,tau2,sig2,phi,results.opts.covmodel);
    cnew = cmat2(1:end-1,end);
    cnew = cnew + (dist(1:end-1,end)==0)*tau2;
    cnew = bs(L,cnew)'; % scale by cmat, cnew = cnew'*cmat^(-1/2)

    ynew = beta(1) + xnew*beta(2:end) + cnew*bs(L,results.y-yhat);
    
    % need to calculate the kriging variance ?
%    xt = bs(L,results.x);
%    dx = [1,xnew]-cnew*xt;
%    tau22 = cnew*cnew' - dx*((xt'*xt)\dx');
%
%    or do we need only the first term, as we already account for beta unc?
    tau22 = cnew*cnew';
    tau20 = sig2 + tau2;

    % more direct formulas to check
%    cnew2 = cmat2(1:end-1,end);
%    dx = [1,xnew] - cnew2'*(cmat\results.x);
%    tau22 = cnew2'*(cmat\cnew2) - dx*((results.x'*(cmat\results.x))\dx');
    
if tau20 - tau22 < -1.0e-13
  warning('something wrong with Kriging variance, please debug')
  keyboard
end

    ynew = ynew + randn(1,1)*sqrt(max(0.0,tau20-tau22));

  end
  ychain(isimu,:) = ynew;
end

out.ychain = ychain;
out.ynew = mean(ychain,1);
out.ystd = std(ychain,1);
