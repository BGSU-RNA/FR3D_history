
function [void] = zShowTable(rowlabels, collabels, table)

ColLabel = [];
for i=1:length(collabels),
  if ~isempty(collabels{i}),
    ColLabel{i} = collabels{i};
  end
end

[y,i] = sort(ColLabel);

fprintf('\n');
fprintf('Code      ');
fprintf(' %4s ',ColLabel{i});
fprintf('  Total');
fprintf('\n');

T = sum(sum(table));                    % total count

for j=1:length(table(:,1)),
  fprintf('%5s     ', rowlabels{j});
  for k=1:length(table(1,:)),
    fprintf(' %4d ', table(j,i(k)));
  end
  fprintf('  %5d', sum(table(j,:)));
  fprintf('  %5.1f%%\n', 100*sum(table(j,:))/T);
end

fprintf('Total     ');
for k=1:length(table(1,:)),
  fprintf('%5.0f ', sum(table(:,i(k))));
end
fprintf('  %4d total count\n',T);

fprintf('Percent   ');
for k=1:length(table(1,:)),
  fprintf('%5.0f%%', 100*sum(table(:,i(k)))/T);
end
fprintf('  %4d total count\n',T);
