% probe for insertions

        aa = a;                            % initial values of these
        BB = B;

        LeftIns = 0;
 
        while sum(G(a,(a+1):B)) == 0,      % no interaction w/in this loop
          if Verbose > 0,
            fprintf('    Insertion %3s\n', File.NT(a).Number);
          end
          LeftIns = LeftIns + 1;           % increase mean number of insertions
          a = a + 1;                               % next base on left
        end

        RightIns = 0;

        while sum(G(a:(B-1),B)) == 0,      % no interaction w/in this loop
          if Verbose > 0,
            fprintf('    Insertion     %4s\n', File.NT(B).Number);
          end
          RightIns = RightIns + 1;         % increase mean number of insertions
          B = B - 1;                               % next base on right
        end
 
        if (LeftIns > 0) || (RightIns > 0),
          if strcmp(Node(n).type,'Basepair'),   % add insertions to basepair
            Node(n).lpar = (0.01+LeftIns)  * [ones(16,1); 0];
            Node(n).rpar = (0.01+RightIns) * [ones(16,1); 0];
            Node(n).leftLengthDist  = subPoisson(0.01+LeftIns);
            Node(n).leftLetterDist  = [0.25 0.25 0.25 0.25];
            Node(n).rightLengthDist = subPoisson(0.01+RightIns); 
            Node(n).rightLetterDist = [0.25 0.25 0.25 0.25];
            Node(n).LeftLetter  = [Node(n).LeftLetter cat(2,File.NT((Node(n).LeftIndex+1):(a-1)).Base)];
            Node(n).RightLetter = [cat(2,File.NT((B+1):(Node(n).RightIndex-1)).Base) Node(n).RightLetter];

          elseif strcmp(Node(n).type,'Cluster'),
            n = n + 1;                          % initial node after cluster
            Node(n).type = 'Initial';
            Node(n).nextnode  = n+1;            % index of next node in tree
            Node(n).lpar = LeftIns;
            Node(n).rpar = RightIns;
            Node(n).LeftIndex  = aa;
            Node(n).RightIndex = BB;
            Node(n).LeftLetter  = cat(2,File.NT(aa:(aa+LeftIns-1)).Base);
            Node(n).RightLetter = cat(2,File.NT((BB-RightIns+1):BB).Base);
            Node(n).leftLengthDist  = subPoisson(Node(n).lpar);
            Node(n).leftLetterDist  = [0.25 0.25 0.25 0.25];
            Node(n).rightLengthDist = subPoisson(Node(n).rpar); 
            Node(n).rightLetterDist = [0.25 0.25 0.25 0.25];

            Node(n).Comment = [' // Initial node ' File.NT(aa).Base File.NT(aa).Number ' - ' File.NT(BB).Base File.NT(BB).Number ' from junction ' Node(n).id];

            if Verbose > 0,
              fprintf('%3d Initial %s:%s and %s:%s\n', n, File.NT(aa).Number, File.NT(aa+LeftIns).Number, File.NT(BB).Number, File.NT(BB+RightIns).Number);
            end
          end
        end
