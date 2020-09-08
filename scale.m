function y = scale(x,m,s)
% scale by mean and std
if nargin<2, m = mean(x); end
if nargin<3, s = std(x); end
y =  bsxfun(@rdivide,center(x,m),s);
