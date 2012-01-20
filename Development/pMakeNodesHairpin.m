% Add a hairpin

        n = n + 1;
        Node(n).type = 'Hairpin';
        Node(n).MiddleIndex = [a:B];
        Node(n).LeftIndex   = a;                 % added 12/7/08.  OK?
        Node(n).RightIndex  = B;                 % added 12/7/08.  OK?
        Node(n).LeftLetter        = cat(2,File.NT(a:B).Base);
        Node(n).RightLetter       = '';
        if a > B,
          Node(n).MiddleIndex = [a];
        end
        Node(n).P       = ones(17,1);
        Node(n).PIns    = 1;
        Node(n).subtype = 'XXXX';                  % revise this later!

        LI = a;                      % use the first
        RI = B;                      % use the last

        Node(n).Comment = [ ' // Hairpin node ' File.NT(LI).Base File.NT(LI).Number ':' File.NT(RI).Base File.NT(RI).Number];

        % if we were adding interactions in hairpins automatically, ...
        % Node(n).InteractionComment{i} = 

        % if we were adding insertions in hairpins automatically, ...
        % Node(n).InsertionComment{i} = 

        if Verbose > 0,
          fprintf('%3d Hairpin   %s:%s %s\n', n, File.NT(a).Number, File.NT(B).Number, cat(2,File.NT(a:B).Base));
        end
        EndLoop = 1;

LR = File.Edge .* (File.Crossing > 2);
[i,j,e] = find(LR(a:B,:));
i = i + a - 1;
[yy,ii] = sort(abs(e));
i = i(ii);
j = j(ii);
e = e(ii);

if Verbose > 1,
  if length(i) > 0,
    fprintf('Hairpin has these long-range interactions: ===================\n');
    for z = 1:length(i),
      NT1 = File.NT(i(z));
      NT2 = File.NT(j(z));
      fprintf('Pair %s %s%5s_%s - %s%5s_%s %s\n', File.Filename, NT1.Base,NT1.Number,NT1.Chain,NT2.Base,NT2.Number,NT2.Chain, zEdgeText(File.Edge(i(z),j(z))));
    end
  else
    fprintf('Hairpin has no long-range interactions ***********************\n');
  end

  c = cat(1,File.NT(1:File.NumNT).Center); % nucleotide centers
  cent = mean(c);

  fprintf('Distance from center of molecule: %8.4f\n', norm(cent-File.NT(a).Center));

  if Verbose > 2,
    fprintf('Press any key to continue\n');
    pause
  end
end
