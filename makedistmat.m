function d= makedistmat(x,y)
% generate distance matrix given data coordinates

n = length(x);
d = zeros(n,n);
for i=1:n
  d(i,:) = sqrt(sum(bsxfun(@minus,[x(i),y(i)],[x,y]).^2,2));
end

