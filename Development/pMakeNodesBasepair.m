% set up basepair

        n = n+1;  
        Node(n).type        = 'Basepair';        % node type
        Node(n).nextnode    = n+1;               % index of next node in tree
        Node(n).LeftLetter  = File.NT(a).Base;
        Node(n).RightLetter = File.NT(b).Base;
        Node(n).Edge        = Interact{a}.Categ(1);
        Node(n).Delete      = 0.01 + (b-a)/1000000;        % deletion prob
        Node(n).lpar        = [.01*ones(16,1); 0]; % left insertion param
        Node(n).rpar        = [.01*ones(16,1); 0]; % right insertion param
        Node(n).LeftIndex   = a;
        Node(n).RightIndex  = b;

        Node(n).SubsProb = pIsoScore(File.Edge(a,b),Node(n).LeftLetter, ...
                                 Node(n).RightLetter,method,ExemplarIDI,ExemplarFreq);

        Node(1).Edge(a,b) = File.Edge(a,b);

        L = Node(n).lpar(1,1);
        R = Node(n).rpar(1,1);
        X = 0:10;     
        Node(n).Z = sum(L.^X*exp(-L)./factorial(X)) ...
                  * sum(R.^X*exp(-R)./factorial(X));
        Node(n).leftLengthDist  = subPoisson(Node(n).lpar(1,1));
        Node(n).leftLetterDist  = [0.25 0.25 0.25 0.25];
        Node(n).rightLengthDist = subPoisson(Node(n).rpar(1,1)); 
        Node(n).rightLetterDist = [0.25 0.25 0.25 0.25];

        Node(n).Comment = [' // Basepair node ' File.NT(a).Base File.NT(a).Number ' - ' File.NT(b).Base File.NT(b).Number ' ' zEdgeText(File.Edge(a,b))];

        if Verbose > 0,
          fprintf('%3d Basepair %4s %4s %s%s %s\n',n, File.NT(a).Number, File.NT(b).Number,File.NT(a).Base,File.NT(b).Base,zEdgeText(File.Edge(a,b)));
        end

        a = a + 1;                                % current base on left
        B = b - 1;

