% zShowCounts(RowLabel,ColLabel,Count,c) displays a table of counts
% Row sums use only columns specified in c
% Labels should all have the same width

function [void] = zShowCounts(RowLabel,ColLabel,Count,c,W)

T = sum(sum(Count(:,c)));                    % total count

for t = 1:length(W),
  if W(t) == 1,
    fprintf('      ');
    for k=1:length(Count(1,:)),
      fprintf('%8s ', ColLabel{k});
    end
    fprintf('RowTotal\n');
    
    for j=1:length(Count(:,1)),
      fprintf('%3s   ', RowLabel{j});
      fprintf('%8d ', Count(j,:));
      fprintf('%8d\n', sum(Count(j,c)));
    end
    
    fprintf('Total ');
    fprintf('%8d ', sum(Count));
    fprintf('%8d\n',T);
  end
  
  if W(t) == 2,
    fprintf('      ');
    for k=1:length(Count(1,:)),
      fprintf('%8s ', ColLabel{k});
    end
    fprintf('Row%%ofTotal\n');
    
    for j=1:length(Count(:,1)),
      fprintf('%3s  ', RowLabel{j});
      fprintf(' %7.1f%%', 100*Count(j,:)/sum(Count(j,c)));
      fprintf('%8.2f%%\n',100*sum(Count(j,c))/T);
    end
    
    fprintf('Percent ');
    fprintf('%5.1f%%   ', 100*sum(Count)/T);
    fprintf('%6d total count\n',T);
  end

  fprintf('\n');
end
