% xListCandidates(File,Candidates,NumToOutput,ModelNumber)

function [] = xListCandidates(File,Search,SortOrder)

Model       = Search.Query;
Candidates  = Search.Candidates;
Discrepancy = Search.Discrepancy;

NumToOutput = 30;

[s,t] = size(Candidates);

N = Model.NumNT;

if SortOrder == 2,
  Candidates = sortrows(Candidates,[N+1 1 2 N+2]);
end

Notify = 1;

switch SortOrder

case 1,

for i=1:s,
  if (Discrepancy(i) >= 0) & (i <= NumToOutput),
    if (Model.Geometric > 0),
      if (Discrepancy(i) > Model.DiscCutoff) & (Notify == 1),
        fprintf('\nSome candidates with discrepancy below %7.4f might not appear below\n\n', Model.RelCutoff);
        Notify = 0;
      end
    end
    f = double(Candidates(i,N+1));
    fprintf('%15s', File(f).Filename);
    if Model.Geometric > 0,
      fprintf('%11.4f',Discrepancy(i));
    end
    for j=1:N
      fprintf('%3s',File(f).NT(Candidates(i,j)).Base);    
      fprintf('%4s',File(f).NT(Candidates(i,j)).Number);    
    end
    if N == 2,
      fprintf('   C1*-C1*: %8.4f', norm(File(f).NT(Candidates(i,1)).Sugar(1,:) - ...
                            File(f).NT(Candidates(i,2)).Sugar(1,:)));
    end
    fprintf('\n');
  end
  LastOK = i;

end

case 2,

for i=1:s,
  if (Discrepancy(i) >= 0) & (i <= 50),
    f = Candidates(i,N+1);
    fprintf('%15s', File(f).Filename);
    if Model.Geometric > 0,
      fprintf('%11.4f',Discrepancy(i));
    end
    for j=1:N
      fprintf('%3s',File(f).NT(Candidates(i,j)).Base);    
      fprintf('%4s',File(f).NT(Candidates(i,j)).Number);    
    end
    fprintf('\n');
  end
  LastOK = i;
end

end



OUT   = strcat(['Query_' Model.Name '_Results.txt']);
fidOUT      = fopen(OUT,'w+');

fprintf(fidOUT,'%s\n','FR3D last results');

fprintf(fidOUT,'Model information\n');
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

fprintf(fidOUT,'  Mask %s\n', Model.Mask);

if isfield(Model,'MaxDiff'),
  fprintf(fidOUT,'  Differences');
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

Notify = 1;
for i=1:s,
  if Discrepancy(i) >= 0,

    if Model.Geometric > 0,
      if (Discrepancy(i) > Model.DiscCutoff) & (Notify == 1),
        fprintf(fidOUT,'\nSome candidates with discrepancy below %7.4f might not appear below\n\n', Model.RelCutoff);
        Notify = 0;
      end
    end
    f = Candidates(i,N+1);
    fprintf(fidOUT,'%15s', File(f).Filename);
    if Model.Geometric > 0,
      fprintf(fidOUT,'%9.4f',Discrepancy(i));
    end
    for j=1:N
      fprintf(fidOUT,'%3s',File(f).NT(Candidates(i,j)).Base);    
      fprintf(fidOUT,'%4s',File(f).NT(Candidates(i,j)).Number);    
    end
    fprintf(fidOUT,'\n');
  end
end

fclose(fidOUT);

if s > 0,
  if LastOK > NumToOutput,
    fprintf('This list is incomplete; see %s for the complete list\n', OUT);
  else
    fprintf('This list is also printed in %s\n',OUT);
  end
end