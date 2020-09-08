function [out,out1] = variogram(x,y,ydata, nout)
% build an empirical variogram model
% [out,out1] = variogram(x,y,ydata, nout)
% nout, how to discretize distances
% out discrtized, averaged output
% out1 original variogram data

% pairs
n = length(x);
np = n*(n-1)/2;

out1 = zeros(np,2);

ii=0;
for i=1:n
  for j=i+1:n
    ii=ii+1;
    d = norm([x(i);y(i)]-[x(j);y(j)]);
    v =(ydata(i)-ydata(j)).^2/2;
    out1(ii,:) = [d,v];
  end
end

dd = linspace(0,max(out1(:,1)),nout);
out = zeros(nout-1,3);
for i=1:nout-1
  out(i,1) = dd(i);
  inds  = find(out1(:,1)>=dd(i)&out1(:,1)<dd(i+1));
  out(i,2) = mean(out1(inds,2));
  out(i,3) = length(inds); % weights for fitting
end

