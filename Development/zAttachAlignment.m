% zAttachAlignment reads the spreadsheet StructureToAlignmentMap to see if there is an entry for File, and if so, stores the variable FastaCol in File

function [File] = zAttachAlignment(File,Verbose)

[n,t,r] = xlsread('Alignments\StructureToAlignmentMap.xls');

for f = 1:length(File),
  p = find(ismember(upper(t(:,1)),upper(File(f).Filename)));

  if length(p) > 0,
    for a = 1:length(p),
      if r{p(a),4} >= 0,
        FASTA = zReadFASTA(['Alignments' filesep r{p(a),3}]);
        if Verbose > 0,
          fprintf('Read alignment %s for %s\n', r{p(a),3}, File(f).Filename);
        end
        File(f) = zAlignToFASTA(File(f),r{p(a),2},FASTA,r{p(a),4},Verbose);

        if Verbose > 0,
          fprintf('Incorporated sequence data with %s\n', File(f).Filename);
        end
      end
    end
  else
    if Verbose > 0,
      fprintf('zAttachAlignment: No alignment found for %s.\n', File(f).Filename);
    end
  end
end