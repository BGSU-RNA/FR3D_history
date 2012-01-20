% zAttachAlignment reads the spreadsheet StructureToAlignmentMap to see if there is an entry for File, and if so, stores the variable FastaCol in File

function [File] = zAttachAlignment(File,Verbose)

[n,t,r] = xlsread(['Alignments' filesep 'StructureToAlignmentMap.xls']);

for f = 1:length(File),
  p = find(ismember(upper(t(:,1)),upper(File(f).Filename)));  % find line(s)
  if length(p) > 0,                                 % if one or more,
    for a = 1:length(p),                            % loop through them
      if r{p(a),4} >= 0,
        FASTA = zReadFASTA(['Alignments' filesep r{p(a),3}]);
        if Verbose > 0,
          fprintf('Read alignment %s for %s\n', r{p(a),3}, File(f).Filename);
        end

        Chain = r{p(a),2};
        Entry = r{p(a),4};

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
