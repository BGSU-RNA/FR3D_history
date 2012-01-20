
load PDBInfo

[a,b] = size(t);

Order = 1:a;

[y,i] = sort(-n(:,2));                    % decreasing order of size
t = t(i,:);
n = n(i,:);

for i = 1:a,
  if strcmp(t{i,1},t{i,10}),              % this represents others
    j = find(ismember(t(:,10),t{i,1}));   % find the others
    for jj = 1:length(j),
      Order(j(jj)) = i + 1/n(j(jj),2);    % largest next
    end
    Order(i) = i;                         % put this one first
  end
end

[y,p] = sort(Order);

t = t(p,:);
n = n(p,:);
Order = Order(p);

for i = 1:a,
  fprintf('%4s %4s %4d %10.6f\n', t{i,1}, t{i,10}, n(i,2), Order(i));
end
