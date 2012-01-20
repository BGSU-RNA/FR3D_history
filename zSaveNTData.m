% zSaveNTData saves File in File.Filename.analyzed

function [void] = zSaveNTData(File)

File.Modified = 0;                             % flag to know to save
File.Distance = [];                            % clear; recompute on load

if ~(exist('PrecomputedData') == 7),        % if directory doesn't yet exist
  mkdir('PrecomputedData');
end

save([pwd filesep 'PrecomputedData' filesep File.Filename '.mat'],'File');
fprintf('Saved %s\n', [File.Filename '.mat']);
