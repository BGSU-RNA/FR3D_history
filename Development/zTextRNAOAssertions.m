% Produce a text output in N3 format of the FR3D-classified interactions in the given RNA molecule

function [T] = zTextRNAOAssertions(File)

r = 1;

RNAO_Base_Code{1} = 'RNAO_000101';
RNAO_Base_Code{2} = 'RNAO_000102';
RNAO_Base_Code{3} = 'RNAO_000103';
RNAO_Base_Code{4} = 'RNAO_000104';

% implement the following as a function, so it can accept negative arguments


RNAO_Interaction_Name{1} = 'pairs_with_cWW';


T{r} = '@prefix rna: <http://purl.obolibrary.org/obo/>.';
r = r + 1;
T{r} = '@prefix ro: <http://purl.obolibrary.org/obo/>.';
r = r + 1;
T{r} = '@prefix : <http://example.org>.';
:m1 rdfs:label "".
:m1 :pdb_code "1A50"
:m1 rdf:type rnao:RNAO_0000168.            # m1 is a molecule or something

for n = 1:length(File.NT),
  r = r + 1;
  T{r} = ['m1 rnao:has_proper_part :nt' num2str(n)];
  r = r + 1;
  T{r} = ['nt' num2str(n) ' rdf:type rnao:' RNAO_Base_Code{File.NT(n).Code} '.'];
end

for n = 1:(length(File.NT)-1),
:nt1 rnao:three_prime_to_five_prime_to :nt2

[i,j] = find((abs(File.Edge) > 0) .* (abs(File.Edge) < 25) .* (triu(File.Edge)));

for k = 1:length(i),
  r = r + 1;
  T{r} = ['nt' num2str(i(k)) ' rnao:' InteractionName{File.Edge(i(k),j(k))} ':nt' num2str(j(k))];
end

