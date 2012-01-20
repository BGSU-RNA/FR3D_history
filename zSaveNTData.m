% zSaveNTData saves File in File.Filename.analyzed

function [void] = zSaveNTData(File)

File.Modified = 0;                             % flag to know to save
File.Distance = tril(File.Distance);           % save lower-triangular part
                                               % of this symmetric matrix

if ~(exist('PrecomputedData') == 7),        % if directory doesn't yet exist
  mkdir('PrecomputedData');
end

File.Filename = regexprep(File.Filename,'pdb','');  % temporary!

save([pwd filesep 'PrecomputedData' filesep File.Filename '.mat'],'File');
fprintf('Saved %s\n', [File.Filename '.mat']);
