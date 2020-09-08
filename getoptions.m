function opts=getoptions(defopts,useropts)
% process option structure and fill it with defaults

if nargin<2 | isempty(useropts), useropts = struct; end

opts = defopts;
ufields = fieldnames(useropts);
deffields = fieldnames(defopts);

for i=1:length(ufields)
 
  if isfield(defopts,ufields{i})
    opts = setfield(opts,ufields{i},getfield(useropts,ufields{i}));
  else
    fprintf('Unknown option %s\n',ufields{i});
    fprintf('Available options are:\n');
    for j=1:length(deffields)
      fprintf('\t%s\n',deffields{j});
    end
    error('Please check your options.');
  end

end

% chech for empty options that need to be set
found = 0;
for i=1:length(deffields)
  if isempty(getfield(opts,deffields{i}))
    fprintf('No default value for option ''%s'' defined.\n',deffields{i});
    found = 1;
  end
end
if found
  error('Please set all needed options.');
end

