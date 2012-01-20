% xWriteCandidates writes candidate list to a file

function [] = xWriteCandidates(File,Search)

Model       = Search.Query;
Candidates  = Search.Candidates;
Discrepancy = Search.Discrepancy;
[s,t]       = size(Candidates);
N           = Model.NumNT;

% ---------------------------------- Write results to a file

if ~(exist('SearchSaveFiles') == 7),        % if directory doesn't yet exist
  mkdir('SearchSaveFiles');
end

OUT    = strcat(['SearchSaveFiles' filesep Search.SaveName '_Results.txt']);
fidOUT = fopen(OUT,'w+');

fprintf(fidOUT,'Query information\n');
if isfield(Model,'Filename'),
  fprintf(fidOUT,'  Filename %s\n', Model.Filename);
  fprintf(fidOUT,'  Nucleotides');
  for j=1:Model.NumNT,
    fprintf(fidOUT,'%3s',Model.NT(j).Base);    
    fprintf(fidOUT,'%4s',Model.NT(j).Number);    
  end
  fprintf(fidOUT,'\n');
  if isfield(Model,'ChainList'),
    fprintf(fidOUT,'  Chain');
    for j=1:Model.NumNT,
      fprintf(fidOUT,'%3s',Model.ChainList{j});    
    end
    fprintf(fidOUT,'\n');
  end
end

if isfield(Model,'Mask'),
  fprintf(fidOUT,'  Mask %s\n', Model.Mask);
end

if isfield(Model,'Diagonal'),
  for j=1:Model.NumNT,
    fprintf(fidOUT,'Diagonal{%d} = %s\n', j, Model.Diagonal{j});
  end
end

if isfield(Model,'MaxDiff'),
  fprintf(fidOUT,'  Maximum differences');
  for j=1:(Model.NumNT-1),
    fprintf(fidOUT,'%4d',Model.MaxDiff(j));    
  end
  fprintf(fidOUT,'\n');
end

if isfield(Model,'MinDiff'),
  fprintf(fidOUT,'  Minimum differences');
  for j=1:(Model.NumNT-1),
    fprintf(fidOUT,'%4d',Model.MaxDiff(j));    
  end
  fprintf(fidOUT,'\n');
end

if isfield(Model,'ReqInter'),
  fprintf(fidOUT,'  Interactions\n');
  for i=1:Model.NumNT,
   for j=(i+1):Model.NumNT,
     if ~isempty(Model.ReqInter{i,j}),
       fprintf(fidOUT,'   Nucleotides %3d and %3d interact as',i,j);
       for k=1:length(Model.ReqInter{i,j}),
         fprintf(fidOUT,' %4.1f',Model.ReqInter{i,j}(k));
       end
       fprintf(fidOUT,'\n');
     end
   end
  end
  fprintf(fidOUT,'\n');
end

if isfield(Model,'Edges'),
  fprintf(fidOUT,'  Edges\n');
  for i=1:Model.NumNT,
   for j=(i+1):Model.NumNT,
     if ~isempty(Model.Edges{i,j}),
       fprintf(fidOUT,'   Edges{%d,%d} = %s\n',i,j,Model.Edges{i,j});
     end
   end
  end
  fprintf(fidOUT,'\n');
end

if Model.Geometric > 0,
  fprintf(fidOUT,'  AngleWeight   '); 
  fprintf(fidOUT,' %4.2f', Model.AngleWeight); 
  fprintf(fidOUT,'\n');
  fprintf(fidOUT,'  LocationWeight'); 
  fprintf(fidOUT,' %4.2f', Model.LocWeight); 
  fprintf(fidOUT,'\n');
  fprintf(fidOUT,'  Guaranteed Cutoff %7.4f\n', Model.DiscCutoff);
  fprintf(fidOUT,'  Relaxed    Cutoff %7.4f\n', Model.RelCutoff);
end
fprintf(fidOUT,'\n');

if (Model.Geometric > 0) & (Model.RelCutoff > Model.DiscCutoff),
  fprintf(fidOUT,'Some motifs with discrepancy between %7.4f and %7.4f might not appear below\n\n', Model.DiscCutoff, Model.RelCutoff);
end


% ----------------------------------- Write candidate list to file

for i=1:s,                                    % loop through candidates
  if Discrepancy(i) >= 0,
    f = Candidates(i,N+1);                    % file number
    fprintf(fidOUT,'%15s', File(f).Filename); % print filename
    if Model.Geometric > 0,
      fprintf(fidOUT,'%9.4f',Discrepancy(i)); % print discrepancy
    end
    for j=1:N
      fprintf(fidOUT,'%3s',File(f).NT(Candidates(i,j)).Base);    
      fprintf(fidOUT,'%4s',File(f).NT(Candidates(i,j)).Number);    
    end
    fprintf(fidOUT,'\n');
  end
end

fclose(fidOUT);
