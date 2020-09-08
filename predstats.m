function out=predstats(ypred,ytrue,fid,names)
%PREDSTATS prediction statistics given the true values

if nargin<3, fid=1; end % standard output

nobs = length(ytrue);

o = [];

for i=1:size(ypred,2)

  resid = ypred(:,i)-ytrue;
  bias  = mean(resid);
  rmse  = sqrt(sum(resid.^2)/nobs);
  ym = mean(ytrue);
  ys = std(ytrue);
  %loorel = resid./abs(ytrue);
  %bias_rel = mean(loorel)*100;
  %rmse_rel = sqrt(sum(loorel.^2)/nobs)*100;
  bias_rel = bias./ym*100;
  rmse_rel = rmse./ym*100;
  Q2 = (1- sum(resid.^2)./sum((ytrue-mean(ytrue)).^2))*100;

  o = [o;[bias,bias_rel,rmse,rmse_rel,Q2]];
end

if nargout>0
  out = o;
end

fprintf(fid,'mean(y): %10.5g\n',ym);
fprintf(fid,'std(y):  %10.5g\n\n',ys);

if nargin>3
  fprintf(fid,'           ');
  for i=1:size(ypred,2)
    fprintf(fid,'%10s',names{i});
  end
  fprintf(fid,'\n');
end

p(fid,'Bias',o(:,1));
p(fid,'Bias%',o(:,2));
p(fid,'RMSE',o(:,3));
p(fid,'RMSE%',o(:,4));
p(fid,'Q2',o(:,5));


function p(fid,s,x);

fprintf(fid,'%10s:',s);
for i=1:length(x)
  fprintf(fid,'% 10.5g',x(i));
end
fprintf(fid,'\n');
