% Produce a text output in N3 format of the FR3D-classified interactions in the given RNA molecule

% File = zAddNTData('1J1U');

% function [T] = zTextRNAOAssertions(File)

clear T

r = 1;

RNAO_Base_Code{1} = 'RNAO_0000104';
RNAO_Base_Code{2} = 'RNAO_0000102';
RNAO_Base_Code{3} = 'RNAO_0000105';
RNAO_Base_Code{4} = 'RNAO_0000103';

% implement the following as a function, so it can accept negative arguments

T{r} = '@prefix rnao: <http://purl.obolibrary.org/obo/>.';
r = r + 1;
T{r} = '@prefix rdfs: <http://www.w3.org/2000/01/rdf-schema#>.';
r = r + 1;
T{r} = '@prefix ro: <http://purl.obolibrary.org/obo/>.';
r = r + 1;
T{r} = '@prefix : <http://example.org>.';
r = r + 1;
T{r} = ':m1 rdfs:label "".';
r = r + 1;
T{r} = [':m1 :pdb_code "' File.Filename '".'];
%r = r + 1;
%T{r} = ':m1 rdf:type rnao:RNAO_0000168.            # m1 is a molecule';

% ---------------------------------------- introduce the nucleotides

for n = 1:length(File.NT),
  r = r + 1;
  T{r} = [':m1 ro:has_proper_part :nt' num2str(n) '.'];
end

% ---------------------------------------- label for each nucleotide

for n = 1:length(File.NT),
  r = r + 1;
%  T{r} = [':nt' num2str(n) ' rdfs:label "nucleotide ' num2str(n) ' in chain ' File.NT(n).Chain ' is ' File.NT(n).Base File.NT(n).Number '".'];
  T{r} = [':nt' num2str(n) ' rdfs:label "' File.NT(n).Base File.NT(n).Number '_' File.NT(n).Chain '".'];
end

% ---------------------------------------- base identity of each nucleotide

for n = 1:length(File.NT),
%  r = r + 1;
%  T{r} = [':nt' num2str(n) ' rdf:type rnao:' RNAO_Base_Code{File.NT(n).Code} '.      # RNA base ' File.NT(n).Base];
end

% ---------------------------------------- nucleotide number from PDB file

for n = 1:length(File.NT),
%  r = r + 1;
%  T{r} = [':nt' num2str(n) ' rnao:nucleotide_number "' File.NT(n).Number '".'];
end

% ---------------------------------------- chain from PDB file

for n = 1:length(File.NT),
%  r = r + 1;
%  T{r} = [':nt' num2str(n) ' rnao:chain "' File.NT(n).Chain '".'];
end

% ---------------------------------------- covalent connectivity - pretend!

for n = 1:(length(File.NT)-1),
%  r = r + 1;
%  T{r} = [':nt ' num2str(n) 'rnao:three_prime_to_five_prime_to :nt' num2str(n+1) '.'];
end

% ---------------------------------------- pairwise interactions

E = abs(File.Edge);

[i,j] = find((E > 0) .* (E < 13));

[y,k] = sort(i);
i = i(k);
j = j(k);

for k = 1:length(i),
  r = r + 1;
  T{r} = [':nt' num2str(i(k)) ' rnao:' zRNAOPairwiseInteractions(File.Edge(i(k),j(k))) ' :nt' num2str(j(k)) '.'];
end

% ---------------------------------------- stacking interactions

[i,j] = find((E > 20) .* (E < 24));

[y,k] = sort(i);
i = i(k);
j = j(k);

for k = 1:length(i),
  r = r + 1;
  T{r} = [':nt' num2str(i(k)) ' rnao:' zRNAOPairwiseInteractions(File.Edge(i(k),j(k))) ' :nt' num2str(j(k)) '.'];
end

% ---------------------------------------- print the output

for r = 1:length(T),
  fprintf('%s\n', T{r});
end

% ---------------------------------------- write the output to a file

fid = fopen([File.Filename '_FR3D.n3'],'w');
for r = 1:length(T),
  fprintf(fid,'%s\n', T{r});
end
fclose(fid);

% zFTPToRNA([File.Filename '_FR3D.n3'],'/home/FR3D/AnalyzedStructures/1J1U');


