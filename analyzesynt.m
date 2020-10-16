%% analyze synthetic data made by makesyntdata.m

opts = struct();
opts.nsimu = 10000;
opts.covmodel = covmod;
opts.pricv = [0.1,1,0.5]; % prior CV for tau2 (known better), sig2, phi
opts.R = diag([0.2,0.2,0.1]);
opts.adaptint = 500;
opts.verbosity = 0;
opts.burnin = 0.5;
opts.alpha_n0 = [5,0.2];
opts.alpha_s20 = [mean(y0), 5].^2;
opts.invphi = 1;
opts.noalpha = 0;

nsimu_pred = 400; % how many predictive simulations
                  % takes some time if nsimu large, 
                  % if 1 uses mean(chain) as parameter

useL = [-1 0 1 2]; % 0 = model with PCA truncation but no spat.corr, 1 = model with spat.corr estimation and PCA truncation, 2 = model with spat.corr but no PCA truncation, -1 = model no spat.corr and no PCA truncation (basic linear regression with MCMC)
useLn = length(useL);

for useLind = 1:useLn
  opts.useL = useL(useLind); % 0 = no L, 1 = L sampling, 2 = fixed L

  resi = residuals(reg(T(:,2:end),y0));
  opts.tau2 = 0.25*var(resi); % tau2, nugget, here its obs unc 
  opts.sig2 = 0.75*var(resi); % sig2, spatial variance
  opts.phi  = sqrt(sum((max(xy_all)-min(xy_all)).^2))*0.01; % phi, the correlation length
  opts.sig  = sqrt(tau2); % not used if useL=1

  % mcmc sampling
  [res,chain,sschain] = spatial_mcmcrun(T,y0,xy(:,1),xy(:,2),opts);

  % to transform from pca to original variables:
  % chain(:,res.betaind(2:end)) = chain(:,res.betaind(2:end))*vx';
  
  if FIG     % plots
    figure(10); 
    spatial_mcmcplot(chain,res,'beta')

    if opts.useL == 0 | opts.useL == 1
        figure(11);
        mcmcplot(chain,res.alphaind);
%        spatial_mcmcplot(chain,res,'alpha')
    end

    if opts.useL == 1 | opts.useL == 2
      figure(12); clf
      mcmcplot(chain,res.thetaind);
      %spatial_mcmcplot(chain,res,'theta')
    
      figure(13); clf
      spatial_mcmcplot(chain,res,'vario',vario)
      hline(sig2)
      hline(tau2)
      hline([],phi)

    elseif opts.useL == 0 | opts.useL == -1 % diagonal spat. covariance matrix
      figure(14);
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
    ypredi = spatial_prediction(res,chain,xnew,xi,yi,nsimu_pred);
    ypred(i) = ypredi.ynew;
    yranks(i) = findrank(ypredi.ychain,y0_all(i));
    if FIG && i==fix(nall/2)
      figure(50); 
      subplot(useLn,1,useLind);
      histp(ypredi.ychain); 
      title(sprintf('predictive distribution for cell id %d, %g%%, L = %d',i,100*yranks(i), opts.useL)); 
      hold on;
      plot(y0_all(i),0,'bo',ypred(i),0,'ro');
      hold off;
      legend('distr','truth','pred');
    end
  end
  if opts.useL == 0 ypred0 = ypred; end 
  if opts.useL == 1 ypred1 = ypred; end 
  if opts.useL == -1 ypred2 = ypred; end 
  if opts.useL == 2 ypred3 = ypred; end 
  disp('done')
 
  if FIG  
    figure(60); 
    subplot(useLn+2,1,1);
    mesh(xx,yy,reshape(y0_all,ng,ng));
    hold on
    plot3(xy(:,1),xy(:,2),y0,'o');
    hold off
    title('true');
 
    subplot(useLn+2,1,useLind+1);
    mesh(xx,yy,reshape(ypred,ng,ng));
    hold on
    plot3(xy(:,1),xy(:,2),y0,'o');
    hold off
    title(sprintf('fitted, L = %d',opts.useL));
  
    figure(70);
    subplot(ceil((useLn+1)/2),2,useLind);
    plot(y0_all,ypred,'ok');
    hold on
    plot([min(y0_all) max(y0_all)],[min(y0_all) max(y0_all)],'color','black');
    plot(y0,ypred(inds),'o','color','black','markerfacecolor','black');
    hold off
    ylabel('predicted');xlabel('observed')
    title(sprintf('y, nobs=%g, nall=%g, L = %d',nobs,nall,opts.useL))
    axis equal;
 
    figure(90);
    subplot(ceil((useLn+1)/2),2,useLind);
    mesh(xx,yy,reshape(ypred-y0_all,ng,ng));
    title(sprintf('residuals, L = %d',opts.useL));
    
    figure(100); 
    subplot(useLn,1,useLind);
    hp=histp(yranks);
    title(sprintf('prediction rank histogram, L = %d',opts.useL));
    hp.FaceColor=[0 0.5 1];
    
    drawnow

  end
 
  cbeta = mean(chain(:,res.betaind))';
  cbeta(2:end) = vx*beta(2:end);
  
  if useL==0
    ypred0=ypred;
  end
end

%% BtSVD2

ww = BtSVD2(Ux,y0,Sx);
ye = Ux_all*ww;
yt = Ux*ww;


figure(60);
subplot(useLn+2,1,useLind+2);
mesh(xx,yy,reshape(ye,ng,ng));
hold on
plot3(xy(:,1),xy(:,2),yt,'o');
hold off
title('fitted, BtSVD2');

figure(70);
subplot(ceil((useLn+1)/2),2,useLind+1);
plot(y0_all,ye,'ok');
hold on
plot([min(y0_all) max(y0_all)],[min(y0_all) max(y0_all)],'color','black');
plot(y0,yt,'o','color','black','markerfacecolor','black');
hold off
ylabel('predicted');xlabel('observed')
title(sprintf('y, nobs=%g, nall=%g, BtSVD2',nobs,nall))
axis equal;

figure(90);
subplot(ceil((useLn+1)/2),2,useLind+1);
mesh(xx,yy,reshape(ye-y0_all,ng,ng));
title('residuals, BtSCVD2');

%% validation set results

fprintf('tau2 = %g, sigma2 = %g, phi = %g, nobs = %d\n',tau2, sig2, phi, nobs)
if useLn > 2
    o = predstats([ypred0(ind_validation),...
               ypred1(ind_validation),...
               ypred2(ind_validation),...
               ypred3(ind_validation),...
               predict(reg(T(:,2:end),y0),Tx_all(ind_validation,2:end)),...
               predict(reg(T(:,2:6),y0),Tx_all(ind_validation,2:6  )),...
               ye(ind_validation)],...
              y0_all(ind_validation), 1, {'MCMCtSVD','MCMCtSVDs','MCMC','MCMCs','PCA','PCA6','BtSVD'});

else
    o = predstats([ypred0(ind_validation),...
               ypred1(ind_validation),...
               predict(reg(T(:,2:end),y0),Tx_all(ind_validation,2:end)),...
               predict(reg(T(:,2:6),y0),Tx_all(ind_validation,2:6  )),...
               ye(ind_validation)],...
              y0_all(ind_validation), 1, {'MCMCtSVD','MCMCtSVDs','PCA','PCA6','BtSVD'});
end