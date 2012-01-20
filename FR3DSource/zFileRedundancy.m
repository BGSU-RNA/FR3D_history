% zFileRedundancy(File,Detail) explores possible redundancy between PDB files

% File = zAddNTData('Nonredundant_3_list',2);
% zFileRedundancy(File,0);              % to list just possible redundancies
% zFileRedundancy(File,3);              % to list all files

% edit the list, then, to remove the ones you've dropped:
% [File, FIndex] = zAddNTData('Nonredundant_3_list',2,File);
% File = File(FIndex);

function [void] = zFileRedundancy(File,Detail)

a = zeros(length(File),7);

for f=1:length(File),
  a(f,1) = length(File(f).NT);
  if ~isempty(File(f).Info.Resolution),
    a(f,2) = File(f).Info.Resolution;
  else
    a(f,2) = 99999;
  end
  E = abs(triu(File(f).Edge));
  a(f,3) = sum(sum((E > 0) .* (E < 16)));
end

[A,i] = sortrows(a,[-1 2 -3]);

File = File(i);                       % sort by decreasing NumNT, resolution

B = [];

j = 1;

while j <= length(i),

 sl = find(A(:,1) == A(j,1));         % files with same number of nucleotides

 if (Detail > 2) || (length(sl) > 1) || (A(j,1) == 0),
   for x = 1:length(sl),
     fprintf('%4s has %4d nucleotides, %4d pairs, ', File(sl(x)).Filename, A(sl(x),1), A(sl(x),3));
     if isempty(File(sl(x)).Info.Resolution),
       fprintf('resolution  ---- ');
     else
       fprintf('resolution %5.2f ', File(sl(x)).Info.Resolution);
     end

     Info = File(sl(x)).Info;

     fprintf(' %s %s %s %s\n', Info.ReleaseDate, Info.Source, Info.Descriptor, Info.Author);
   end
 end

 if (A(j,1) > 0) && (length(sl) > 1),

  fprintf('Sequence agreement measure\n');
  fprintf('           ');
  for x = 1:length(sl),
   fprintf(' %4s ', File(sl(x)).Filename);
  end
  fprintf('\n');
  for x = 1:length(sl),
   xB = [];
   for k = 1:A(sl(x),1),
    xB = [xB File(sl(x)).NT(k).Code];           % store bases for this file
   end

   fprintf('%4s', File(sl(x)).Filename);
   if isempty(File(sl(x)).Info.Resolution),
     fprintf(' ---- ');
   else
     fprintf('%5.2f ', File(sl(x)).Info.Resolution);
   end

   for y = 1:length(sl),
    yB = [];
    for k = 1:A(sl(y),1),
     yB = [yB File(sl(y)).NT(k).Code];           % store bases for this file
    end

    if x == y,
     fprintf('      ');
    else
     fprintf(' %5.0f', 100*sum(xB == yB)/length(xB));
    end
   end

   fprintf(' %s\n', File(sl(x)).Info.Source);
  end
         
  fprintf('Geometric discrepancy\n');
  fprintf('             ');
  for x = 1:length(sl),
   fprintf(' %4s   ', File(sl(x)).Filename);
  end
  fprintf('\n');
  for x = 1:length(sl),
    fprintf('%4s', File(sl(x)).Filename);
    if isempty(File(sl(x)).Info.Resolution),
      fprintf(' ---- ');
    else
      fprintf('%5.2f ', File(sl(x)).Info.Resolution);
    end
    for y = 1:length(sl),
      d = xDiscrepancy(File(sl(x)),1:A(sl(x),1),File(sl(y)),1:A(sl(y),1));
      if x==y,
       fprintf('        ');
      else
       fprintf(' %7.4f',d);
      end
    end
    fprintf(' %s\n', File(sl(x)).Info.Source);
  end
  fprintf('\n');

 end

 j = j + length(sl);

end


