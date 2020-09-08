function [res,chain,sschain] = spatial_mcmcrun(x,y,xcoord,ycoord,uopts)
% mcmc for regression model with hyper and spatial parameters

if size(y,2)>1
  error('Works only for one response, for now');
end

% the default values for the options
opts.nsimu = 10000;
opts.tau2 = 0.1*var(y);
opts.sig2 = var(y);
opts.sig = 1;
opts.phi = [];
opts.diagsig = 1;
opts.useL = 1;
opts.adaptint = 100;
opts.alpha = [1/400.^2 1/400^2]; % initial alpha hyperparameter
opts.FIG = 0;
opts.verbosity = 0;
opts.burnin = 0;
opts.invphi = 1; % sample 1/phi
opts.phiminmax = NaN; % uniform prior for 1/phi
opts.noalpha = 0;
opts.covmodel = 'gaussian';
opts.waitbar = 1;

% new prior parameterization as n0 and s02 (like in mcmcrun)
% s20 is the prior mean for 1/alpha ie for variance.
% n0 is the accuracy of prior mean as virtual observations
opts.alpha_n0  = [2 0];
opts.alpha_s20 = [mean(y) 20].^2;

% prior for diagonal error variance (if diagsig = 1);
opts.N0 = 1;
opts.S20 = -1; %

opts.R = diag(3)*0.1;
opts.pricv =  1;

% check and process options
opts = getoptions(opts,uopts);

% fix
if isnan(opts.phiminmax), opts.phiminmax=[];end
if opts.useL == -1
  opts.useL = 0;
  opts.noalpha = 1;
end

% fix some options
if opts.useL == 1, opts.diagsig = 0; end
if opts.S20  <= 0, opts.S20 = opts.sig.^2; end

n = size(x,1);

if opts.useL

  dist = makedistmat(xcoord,ycoord);
  [cmat,L] = makecovmat2(dist,opts.tau2,opts.sig2,opts.phi,opts.covmodel);

  L0 = chol(cmat,'lower');

  lam = 10;lam = 2;lam = 1;
  L = L0 / lam;
  opts.sig = 1;
  ssold = Inf;
  priold = 0;
  parold = [opts.tau2,opts.sig2,opts.phi];
  par0 = parold;

  if length(opts.pricv)==1
    opts.pricv = [opts.pricv(1),opts.pricv(1),opts.pricv(1)];
  end

  %
  covcholfun=@(par)chol(makecovmat2(dist,par(1),par(2),par(3),opts.covmodel),'lower');
    if opts.invphi == 1
      if isempty(opts.phiminmax)
        priorfun = @(p) (log(p(1)/opts.tau2)/opts.pricv(1)).^2 + ...
            (log(p(2)/opts.sig2)/opts.pricv(2)).^2 + ...
            (log(p(3)*opts.phi)/opts.pricv(3)).^2;
      else
        priorfun = @(p) (log(p(1)/opts.tau2)/opts.pricv(1)).^2 + ...
            (log(p(2)/opts.sig2)/opts.pricv(2)).^2 + ...
            (1/p(3)<opts.phiminmax(1)|1/p(3)>opts.phiminmax(2))*1e99;
      end
    else
      priorfun = @(p) (log(p(1)/opts.tau2)/opts.pricv(1)).^2 + ...
          (log(p(2)/opts.sig2)/opts.pricv(2)).^2 + ...
          (log(p(3)/opts.phi)/opts.pricv(3)).^2;
    end

end

p = size(x,2)-1; % length of beta - 1

if opts.useL==1
  npar = p + 6;
else
  npar = p + 4;
end

% generate names
names{1} = 'intercept';
for i = 2:p+1 %
  names{i} = sprintf('\\beta_{%d}',i-1);
end
names{i+1} = '\sigma_{0}';     % save results as stds
names{i+2} = '\sigma_{\beta}';
if opts.diagsig
  names{i+3} = '\sigma_Y';
else
  names{i+3} = '\tau^2';
  names{i+4} = '\sigma^2';
  names{i+5} = '\phi';
end

chain = zeros(opts.nsimu,npar);
sschain = zeros(opts.nsimu,1);

% initial values from options
alpha = opts.alpha;
R = opts.R;

if opts.waitbar
  hwb = waitbar(0,'doing MCMC, please wait');
end

sig = opts.sig;
nacce = 1;

for isimu = 1:opts.nsimu

  % beta|alpha
  if opts.useL
    if opts.noalpha     % no alpha
      ttt = L\x/sig;
      mu = ttt\[L\y/sig];
    else
      ttt = [L\x/sig;diag([sqrt(alpha(1)),ones(1,p)*sqrt(alpha(2))])];
      mu = ttt\[L\y/sig;zeros(p+1,1)];
    end
  else
    % no L
    if opts.noalpha     % no alpha
      ttt = x/sig;
      mu = ttt\[y/sig];
    else
      ttt = [x/sig;diag([sqrt(alpha(1)),ones(1,p)*sqrt(alpha(2))])];
      mu = ttt\[y/sig;zeros(p+1,1)];
    end
  end

  Sigi = ttt'*ttt;
  beta = mvnorr(1,mu,Sigi,inv(chol(Sigi)));

  chain(isimu,1) = beta(1);
  chain(isimu,2:length(beta)) = (beta(2:end)')';

  % alpha|beta
  % new parameterization
  if not(opts.noalpha)
    alpha(1) = gammar(1,1,(1+opts.alpha_n0(1))/2,...
                      2/(beta(1).^2+opts.alpha_n0(1)*opts.alpha_s20(1)));
    alpha(2) = gammar(1,1,(p+opts.alpha_n0(2))/2,...
                      2/(sum(beta(2:end).^2)+opts.alpha_n0(2)*opts.alpha_s20(2)));

    chain(isimu,p+2) = sqrt(1/alpha(1));
    chain(isimu,p+3) = sqrt(1/alpha(2));
  end

  if opts.diagsig
    % update sig from conjugate gamma distribution
    % n is number of observations
    % ss is the residual sum of squares
    ss  = sum((y-x*beta').^2);
    sig = sqrt(1./gammar(1,1,(opts.N0+n)/2,2./(opts.N0*opts.S20+ss)));
    chain(isimu,end) = sig;

  elseif opts.useL == 1
    % three spatial parameters
    if opts.invphi == 1
      parold(3) = 1/parold(3);
    end
    if isimu == 1
      parnew  = parold; % logN proposal
    else
      parnew  = exp(log(parold) + randn(1,3)*R); % logN proposal
    end
    prinew = priorfun(parnew);
    if opts.invphi == 1
      parnew(3) = 1/parnew(3);
      parold(3) = 1/parold(3);
    end
    Lnew = covcholfun(parnew);
    ssold = norm(L\(y-x*beta')).^2 + 2*sum(log(diag(L)));
    ssnew = norm(Lnew\(y-x*beta')).^2 + 2*sum(log(diag(Lnew)));
    tst = exp(-0.5*(ssnew-ssold + prinew-priold));

    if opts.verbosity>0
      fprintf('%g, %g, %g, %g, %g\n',ssold,ssnew,tst,priold,prinew);
      fprintf('parnew: %g %g %g\n',parnew);
      fprintf('parold: %g %g %g\n',parold);
    end
    if rand<tst % accept
      ssold = ssnew;
      parold = parnew;
      priold = prinew;
      L = Lnew;
      nacce = nacce + 1;
    end
    chain(isimu,end-2) = parold(1);
    chain(isimu,end-1) = parold(2);
    chain(isimu,end)   = parold(3);
    sschain(isimu)     = ssold;

    % adaptation
    if opts.adaptint>0 & fix(isimu/opts.adaptint) == isimu/opts.adaptint
      cnew = cov(log(chain(1:isimu,end-2:end)));
      [Ra,is] = chol(cnew + eye(3)*1e-10);
      if is
        fprintf('Warning cmat singular, not adapting\n');
      else
        R = Ra/2.4*sqrt(3);
      end
    end

  end
  if opts.waitbar & (fix(isimu/200) == isimu/200)
    waitbar(isimu/opts.nsimu,hwb)
  end
end

if opts.burnin > 0
  if opts.burnin > 1
    chain = chain(opts.burnin-1:end,:);
    sschain = sschain(opts.burnin-1:end,:);
  else
    chain = chain(ceil(opts.burnin*opts.nsimu):end,:);
    sschain = sschain(ceil(opts.burnin*opts.nsimu):end,:);
  end
end

if opts.waitbar
  delete(hwb);
end


if opts.FIG

  figure(301); clf
  o = chainstats(chain,names);
  tstat = o(2:end-3-2,1)./o(2:end-3-2,2); % t statistics for beta(1) - beta(p)
  plot(tstat,'-o'); hline(2); hline(-2);
  title('t statistics for regression parameters');
  ylabel('t'); xlabel('\beta index')

end

res.class = 'spmcmc';
res.nsimu = size(chain,1);
res.acce = nacce/opts.nsimu;
res.names = names;
if opts.useL
  res.priorfun = priorfun;
  res.covcholfun = covcholfun;
end
res.alphaind = (p+2):(p+3);
res.betaind = 1:(p+1);
if opts.useL == 1
  res.thetaind = (p+4):(p+6);
else
  res.thetaind = [];
end
res.xcoord = xcoord;
res.ycoord = ycoord;
res.x = x;
res.y = y;
res.R = R;
res.opts = opts;
