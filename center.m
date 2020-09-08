function y = center(x,m)
% center by removing the mean from each column
if nargin<2, m = mean(x); end
y =  bsxfun(@minus,x,m);
