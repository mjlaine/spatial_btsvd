function y=findrank(x,x0)
% y=findrank(x,x0) find the rank of x0 on vector x as pct

n = length(x);
i = find(x0<=sort(x),1)-1;
if isempty(i)
  y = 1;
else
  y = i/n;
end
