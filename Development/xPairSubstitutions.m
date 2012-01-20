
%Filenames = {'2avy','2aw4'};
%File = zAddNTData(Filenames);

%f = 1;

%i = zIndexLookup(File(f),'55');
%j = zIndexLookup(File(f),'56');

function [D,D1] = xPairSubstitutions(File,f,i,j)

for a = 1:length(File),
  Filenames{a} = File(a).Filename;
end

Verbose = 0;

Query.Name           = 'Pairwise interaction search';
Query.Description    = 'Pairwise interaction search';
Query.Filename       = File(f).Filename;
Query.NTList         = {File(f).NT(i).Number File(f).NT(j).Number};
Query.ChainList      = {File(f).NT(i).Chain File(f).NT(j).Chain}; 
Query.Edges{1,2}     = [zEdgeText(File(f).Edge(i,j)) ' n' zEdgeText(File(f).Edge(i,j))];
Query.DiscCutoff     = 0.9;      

Query = xConstructQuery(Query,File(f));

UsingLibrary = 1;                  

xFR3DSearch

L = length(Candidates(:,1));        % number of candidates found

D = zeros(L,3);

for c = 1:L,
  D(c,1) = File(Candidates(c,3)).NT(Candidates(c,1)).Code;
  D(c,2) = File(Candidates(c,3)).NT(Candidates(c,2)).Code;
  D(c,3) = Discrepancy(c);
end

D1 = D;

%Search
%Search.Query

[Discrepancy, Candidates] = xRankCandidatesIDI(File,Search.Query,Candidates,Verbose);

D = zeros(L,3);

for c = 1:L,
  D(c,1) = File(Candidates(c,3)).NT(Candidates(c,1)).Code;
  D(c,2) = File(Candidates(c,3)).NT(Candidates(c,2)).Code;
  D(c,3) = Discrepancy(c);
end

% xListCandidates(Search,20);
