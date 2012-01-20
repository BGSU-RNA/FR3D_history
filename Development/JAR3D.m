
function [void] = JAR3D(FastaFile,ModelFile,NumSequences,DNA,Range)

if nargin < 3,
  NumSequences = 10;
end

if nargin < 4,
  DNA = 0;
end

if nargin < 5,
  Range = 20;
end

JAR3D_path

C = ['java -cp "' JAR3Dpath filesep 'bin" JAR3DMoleculeAligner "' pwd '" ' FastaFile ' ' ModelFile ' ' num2str(NumSequences) ' ' num2str(DNA) ' ' num2str(Range)]

system(C);







%['java -cp "' JAR3Dpath filesep 'bin" JAR3DMoleculeAligner "' pwd '" ' FastaFile ' ' ModelFile]);



% system(['java -cp "' JAR3Dpath filesep 'bin" JAR3DMoleculeAligner "' pwd '" 16S_sequences_from_1j5e_2AVY.fasta 16S_from_2AVY.txt']);

% system(['java -cp "' JAR3Dpath filesep 'bin" JAR3DMoleculeAligner "' pwd '" 16S_sequence_from_2avy.fasta 16S_from_2AVY.txt']);


