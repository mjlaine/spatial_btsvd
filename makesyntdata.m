%% Synthetic data


% no difference in methods if:
% phi too big >0.5  of too small <0.05
% sig2/tau2 small < 2

% test cases to run
% tau2 = 0.1^2
% sig2 = 0.5^2 1.0^2 2.0^2
% phi  = 0.2 0.5
% nobs = 20 50 100

ng = 20; % grid size ng*ng
nall = ng*ng;
ccond = 30; % cond number for X rotation, 30 is OK

if ~exist('INLOOP','var') || INLOOP==0 % these are set in simulation loop, also
  nobs = 30; %  number of observations
  tau2 = 0.1^2; % obs error
%  tau2 = 0.5^2; % obs error
%  sig2 = 0.2^2; % spatial error % sig2/tau2 small ->  no sigf. dep.
  sig2 = 1.0^2; % spatial error
%  sig2 = 2.0^2; % spatial error % sig2/tau2 big -> spatial dependence
%  phi = 0.5; % 0.5 too big (?)
  phi = 0.1;    % spatial correlation length parameter
%  phi = 0.001;    % spatial correlation length parameter
  FIG = 1;
end

nxpar = 10; % number of extra X columns
sigX = 0.01; % std of extra X columns
covmod = 'exponential'; % spatial covariance model 
%covmod = 'gaussian';

% data coordinates on regular grid
[ii,jj]=ndgrid(1:ng);
xy_all = [ii(:),jj(:)]./ng - 0.5;
xx = (1:ng)'./ng-0.5;
yy = (1:ng)'./ng-0.5;

% predictors based on coordinates
X0_all = [xy_all, xy_all(:,1).*xy_all(:,2), xy_all(:,1).^2, xy_all(:,2).^2];
% add extra dummy X columns
X0_all = [X0_all,randn(nall,nxpar)*sigX];

% rotate X to make the columns correlate
C = covcond(ccond,ones(size(X0_all,2),1));
[V,D]=eig(C);
X0_all = X0_all*V';
% add intercept
X_all = [ones(nall,1),X0_all];
npar = size(X_all,2);

beta = zeros(npar,1);
beta(1) = 5; % intercept
beta(2:6) = [1 2 3 3 2]'; % quadratic model
beta(7:end) = 0; % 
% apply X rotation to beta
beta0=beta;
beta(2:end) = V*beta(2:end);

% PCA of full X
meanX = mean(X0_all);
stdX = std(X0_all);
[ux,sx,vx] = svd(scale(X0_all),0);
Tx0_all = ux*sx;
Tx_all = [ones(nall,1) Tx0_all];
Ux_all = [ones(nall,1) ux];

% add spatial correlation to observations
dmat = makedistmat(xy_all(:,1),xy_all(:,2));
[cmat,L] = makecovmat2(dmat,tau2,sig2,phi,covmod);
e =  mvnorr(1,zeros(1,nall),cmat,L')';

% generate observations
y00_all = X_all*beta;
y0_all  = X_all*beta + e;

% select nobs observations at random
inds = randperm(nall,nobs);
X  = X_all(inds,:);
X0 = X0_all(inds,:);
y0 = y0_all(inds,:);
xy = xy_all(inds,:);
T = Tx_all(inds,:);

% virpin alkaa...
Ux = Ux_all(inds,:);
Sx = diag(sx);
% virpin loppuu

% linear fit 
betahat0 = X\y0;
yhat = X*betahat0;
resi = y0-yhat;
betahat = T\y0;
yhat = T*betahat;
resi = y0-yhat;

figure(12);clf
surf(xx,yy,reshape(y0_all,ng,ng));
title('all observations')

figure(1); clf
plot3(xy(:,1),xy(:,2),resi,'o');
grid
zl = zlim;
hold on
plot3(xy(:,1),xy(:,2),ones(nobs,1)*zl(1),'o','markerfacecolor','black');
hold off
title('observations')

% direct fit vs pca fit, should be the same
% [beta,betahat0,[betahat(1);vx*betahat(2:end)]]

% variogram for plot
[vario,vario_raw] = variogram(xy(:,1),xy(:,2),resi, 100);

figure(5); clf
plot(vario(:,1),vario(:,2),'o','markerfacecolor',[0.0 0.5 1.0],'color',[0.0 0.5 1.0])
hline(sig2)
hline(tau2)
hline([],phi)
title('residual variogram')

figure(14); clf
mesh(yy,xx,reshape(e,ng,ng));
title('obs error')

% return

