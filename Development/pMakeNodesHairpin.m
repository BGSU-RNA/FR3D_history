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
          fprintf('%3d Hairpin  %s:%s\n', n, File.NT(a).Number, File.NT(B).Number);
        end
        EndLoop = 1;
