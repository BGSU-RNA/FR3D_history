
function [void] = JAR3D(FastaFile,ModelFile)

JAR3D_path

system(['java -cp "' JAR3Dpath filesep 'bin" JAR3DMoleculeAligner "' pwd '" ' FastaFile ' ' ModelFile]);



% system(['java -cp "' JAR3Dpath filesep 'bin" JAR3DMoleculeAligner "' pwd '" 16S_sequences_from_1j5e_2AVY.fasta 16S_from_2AVY.txt']);

% system(['java -cp "' JAR3Dpath filesep 'bin" JAR3DMoleculeAligner "' pwd '" 16S_sequence_from_2avy.fasta 16S_from_2AVY.txt']);


