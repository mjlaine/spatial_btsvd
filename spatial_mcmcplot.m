function out=spatial_mcmcplot(chain,res,type,varargin)
% some plots for spatial_mcmcrun results

if nargin < 3
  type = 'theta';
end

blues = [0.0 0.5 1.0]; %

switch type

  case 'chain'

  case 'beta'

    mcmcplot(chain,res.betaind,res.names,'chainpanel');

  case 'alpha'

%    mcmcplot(chain,res.alphaind,res.names,'chainpanel');
    mcmcplot(chain,res.alphaind,res.names,'denspanel');
    subplot(2,1,1); xl=xlim;xl(1) = max(0,xl(1));xlim(xl);
    hold on
    xx = linspace(xl(1),xl(2));
%    plot(xx,invchi1pf(xx,res.opts.alpha_n0(1),sqrt(res.opts.alpha_s20(1))),'--g');
%res.opts.alpha_n0
    plot(xx,invchipf(xx,res.opts.alpha_n0(1),res.opts.alpha_s20(1)),'--g');
    hold off
    subplot(2,1,2); xl=xlim;xl(1) = max(0,xl(1));xlim(xl)
    hold on
    xx = linspace(xl(1),xl(2));
%    plot(xx,invchi1pf(xx,res.opts.alpha_n0(2),sqrt(res.opts.alpha_s20(2))),'--g');
    plot(xx,invchipf(xx,res.opts.alpha_n0(2),res.opts.alpha_s20(2)),'--g');
    hold off

  case 'theta'

%    mcmcplot(chain,res.thetaind,res.names,'denspanel',1.5);
    mcmcplot(chain,res.thetaind,res.names,'hist', 30);
    subplot(2,2,1); xl=xlim;xl(1) = max(0,xl(1));xlim(xl);
    hold on
    xx = linspace(xl(1),xl(2));
    plot(xx,lognorpf(xx,log(res.opts.tau2),res.opts.pricv(1)),'--b');
    hold off

    subplot(2,2,2); xl=xlim;xl(1) = max(0,xl(1));xlim(xl);
    hold on
    xx = linspace(xl(1),xl(2));
    plot(xx,lognorpf(xx,log(res.opts.sig2),res.opts.pricv(2)),'--b');
    hold off

    subplot(2,2,3); xl=xlim;xl(1) = max(0,xl(1));xlim(xl);
    if isempty(res.opts.phiminmax)
      hold on
      xx = linspace(xl(1),xl(2));
      plot(xx,lognorpf(xx,log(res.opts.phi),res.opts.pricv(3)),'--b');
      hold off
    end

  case 'thetapairs'
    mcmcplot(chain,res.thetaind,res.names,'pairs');
  
  case 'thetachain'
    mcmcplot(chain,res.thetaind,res.names,'chainpanel');

  case 'vario'

    nsimu = 100;
    vario = varargin{1};
    switch res.opts.covmodel
      case 'gaussian'
        vmodel = @(th,x) th(1) + (th(2)-th(1)).*(1-exp(-0.5*(x./th(3)).^2));
      case 'exponential'
        vmodel = @(th,x) th(1) + (th(2)-th(1)).*(1-exp(-(abs(x)./th(3))));
    end
    xx = linspace(0,max(vario(:,1)));
    %  sample variograms from the chain
    for isimu=1:nsimu
      vpars = chain(ceil(rand*size(chain,1)),res.thetaind);
      hold on
      plot(xx,vmodel(vpars,xx),'color',[0.8 0.8 0.8]);
      hold off
    end

    hold on
    plot(vario(:,1),vario(:,2),'o','markerfacecolor',blues,'color',blues)
    hold off


  otherwise
    error('unknown plot type')
end