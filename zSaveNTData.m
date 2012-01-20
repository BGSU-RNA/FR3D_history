% zSaveNTData saves File in File.Filename.analyzed

function [void] = zSaveNTData(File)

File.Distance = tril(File.Distance);           % save lower-triangular part
                                               % of this symmetric matrix

if ~(exist('PrecomputedData') == 7),        % if directory doesn't yet exist
  mkdir('PrecomputedData');
end

save([pwd filesep 'PrecomputedData' filesep File.Filename '.mat'],'File');
fprintf('Saved %s\n', [File.Filename '.mat']);
