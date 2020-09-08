%% analyze synthetic data made by makesyntdata.m

opts = struct();
opts.nsimu = 10000;
opts.covmodel = covmod;
opts.tau2 = tau2; % 0.25*var(resi); % tau2, nugget, here its obs unc 
opts.sig2 = sig2; % 0.75*var(resi); % sig2, spatial variance
opts.phi  = phi; % phi, the correlation length
opts.sig  = sqrt(tau2); % not used if useL=1
opts.pricv = [0.1,1,0.5]; % prior CV for tau2 (known better), sig2, phi
% opts.pricv = [1,1,1]; % prior CV for tau2 (known better), sig2, phi
opts.R = diag([0.2,0.2,0.1]);
opts.adaptint = 500;
opts.verbosity = 0;
opts.burnin = 0.5;
%opts.alpha_n0 = [10,0.2];
opts.alpha_n0 = [5,0.2];
%opts.alpha_s20 = [mean(y0)*0.5, 5].^2;
opts.alpha_s20 = [mean(y0), 5].^2;
opts.invphi = 1;
%  opts.phiminmax = [0.1 600];
%opts.phiminmax = [0.00001 5];

opts.noalpha = 0;

nsimu_pred = 100; % how many predictive simulations
                  % takes some time if nsimu large, 
                  % if 1 uses mean(chain) as parameter

if exist('INLOOP','var') && INLOOP==1 % in simulation loop
  opts.waitbar = 0;
end

for useL=0:1
%for useL=1
  opts.useL = useL; % 0 = no L, 1 = L sampling, 2 = fixed L

  resi = residuals(reg(T(:,2:end  ),y0));
  opts.tau2 = 0.25*var(resi); % tau2, nugget, here its obs unc 
  opts.sig2 = 0.75*var(resi); % sig2, spatial variance
  
  % mcmc sampling
  [res,chain,sschain] = spatial_mcmcrun(T,y0,xy(:,1),xy(:,2),opts);

  % from pca to original variables
  %chain(:,res.betaind(2:end)) = chain(:,res.betaind(2:end))*vx';
  
  if FIG 
    % plots
    %  figure(3); clf
    figure(3 + abs(100*useL)); clf
    spatial_mcmcplot(chain,res,'beta')
    %  figure(6); clf
    figure(6 + abs(10*useL)); clf
    spatial_mcmcplot(chain,res,'alpha')
    if opts.useL == 1
      figure(4); clf
      spatial_mcmcplot(chain,res,'theta')
    
      figure(5); clf
      spatial_mcmcplot(chain,res,'vario',vario)
      hline(sig2)
      hline(tau2)
      hline([],phi)
    elseif opts.useL == 0
      figure(4); clf
      mcmcplot(chain,size(chain,2),res.names,'hist', 30);
      xl=xlim;xl(1) = max(0,xl(1));xlim(xl);
      hold on
      xxxx = linspace(xl(1),xl(2));
      plot(xxxx,invchi1pf(xxxx,res.opts.N0,sqrt(res.opts.S20)),'--b');
      hold off
    end
    drawnow
  end
  
  %% predictions
  ypred = zeros(nall,1);
  ypredstd = zeros(nall,1);
  yranks = zeros(nall,1);

  disp('calculating predictions ...')
%  nsimu = 100; % takes some time if nsimu large, if 1 uses mean(chain) as parameters
  for i=1:nall
    xi = xy_all(i,1);
    yi = xy_all(i,2);
    xnew = Tx_all(i,2:end);
%   ypred(i) = getfield(spatial_prediction(res,chain,xnew,xi,yi,nsimu_pred),'ynew');
    ypredi = spatial_prediction(res,chain,xnew,xi,yi,nsimu_pred);
    ypred(i) = ypredi.ynew;
    %ypredstd(i) = ypredi.ystd;
    yranks(i) = findrank(ypredi.ychain,y0_all(i));
    if FIG && i==fix(nall/2)
      figure(221+abs(useL)); clf
      histp(ypredi.ychain); title(sprintf('predictive distribution, %g%%',100*yranks(i))); hline([],y0_all(i));
    end
  end
  disp('done')
 
  if FIG  
    figure(7); clf
    subplot(2,1,1)
    mesh(xx,yy,reshape(y0_all,ng,ng));
    hold on
    plot3(xy(:,1),xy(:,2),y0,'o');
    hold off
    title('true');
  
    subplot(2,1,2)
    mesh(xx,yy,reshape(ypred,ng,ng));
    hold on
    plot3(xy(:,1),xy(:,2),y0,'o');
    hold off
    title('fitted')
  
    figure(13+abs(100*useL));clf
    plot(y0_all,ypred,'ok');
    hold on
    plot([min(y0_all) max(y0_all)],[min(y0_all) max(y0_all)],'color','black');
    plot(y0,ypred(inds),'o','color','black','markerfacecolor','black');
    hold off
    ylabel('predicted');xlabel('observed')
    if useL==0
      title(sprintf('Observation vs. prediction, nobs=%g, nall=%g',nobs,nall))
    else
      title(sprintf('Observation vs. prediction, nobs=%g, nall=%g (spatial)',nobs,nall))
      figure(13);
      hold on
      plot(y0_all,ypred,'o','color','blue');
      plot(y0,ypred(inds),'o','color','blue','markerfacecolor','blue');
      hold off
    end    
    figure(14);clf
    mesh(xx,yy,reshape(ypred-y0_all,ng,ng));
    
    figure(15+100*abs(useL)); clf
    hp=histp(yranks);
    title('prediction rank histogram')
    hp.FaceColor=[0 0.5 1];
    
    drawnow

  end
 
  cbeta = mean(chain(:,res.betaind))';
  cbeta(2:end) = vx*beta(2:end);
%  [beta,betahat0,[betahat(1);vx*betahat(2:end)],cbeta]
  
  if useL==0
    ypred0=ypred;
%   ypredstd0=ypredstd;
  end
end



[ww,aaa,bbb] = BtSVD2(Ux,y0,Sx);
ye = Ux_all*ww;

o = predstats([ypred0,...
               ypred,...
               predict(reg(T(:,2:end),y0),Tx_all(:,2:end)),...
               predict(reg(T(:,2:6  ),y0),Tx_all(:,2:6  )),...
               ye],...
              y0_all, 1, {'MCMC','MCMCsp','PCA','PCA6','BtSVD'});

