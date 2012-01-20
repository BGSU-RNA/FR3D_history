% zWriteHTMLFileList produces large tables of all PDB files with RNA and with only the non-redundant dataset

% Format = 1 writes the page all RNA-containing PDB files
% Format = 2 writes the page of non-redundant PDB files

NRDate = '2010-05-19';

load PDBInfo

for r = 1:length(t(:,1)),
  if isempty(t{r,10}),
    t{r,10} = t{r,1};
  end
end

Reports = [[1 Inf]; [2 20]; [2 4]; [2 3.5]; [2 3]; [2 2.5]; [2 2]; [2 1.5]];

for r = 1:length(Reports(:,1)),

Format = Reports(r,1);
MaxRes = Reports(r,2);

switch Format
case 1, 
  i = find(n(:,2)>0);
  [tt,nn] = zEquivalents(t(i,:),n(i,:),5,0);        % 
  fid = fopen('Web\AnalyzedStructures\index.html','w');
  fprintf(fid,'<html>');
  fprintf(fid,'<head>\n');
  fprintf(fid,'<title>List of all PDB structures containing RNA nucleotides</title>\n');
  fprintf(fid,'</head>\n');
  fprintf(fid,'<body>\n');
  fprintf(fid,'<h3>List of all RNA-containing structures</h3>\n');
  fprintf(fid,'Equivalent structures are grouped together.<br>\n');
  fprintf(fid,'Click here for the non-redundant list of high-resolution x-ray structures up to: ');
  for a = 2:length(Reports(:,1)),
    MR = num2str(Reports(a,2));
    MRText = strrep(num2str(MR),'.',',');
    fprintf(fid,'%s',['<a href="Nonredundant_' MRText 'A.html">' MR 'A</a> ']);
  end
  fprintf(fid,'\n');

case 2,
  i = find((n(:,1) > 0) .* (n(:,1) <= MaxRes) .* (n(:,3) > 0));   
                                % Not NMR, res to MaxRes, at least one pair
  [tt,nn] = zEquivalents(t(i,:),n(i,:),5,0);      % get equivalents

  Keep = [];
  Keep(1) = 1;
  for j = 2:length(tt(:,1)),
    if ~strcmp(tt{j,10},tt{j-1,10}),              % first instance from group
      Keep(j) = 1;
    end
  end

  i = find(Keep);
  tt = tt(i,:);
  nn = nn(i,:);

  MRText = strrep(num2str(MaxRes),'.',',');
  fid = fopen(['Web\AnalyzedStructures\Nonredundant_' MRText 'A.html'],'w');
  fprintf(fid,'<html>');
  fprintf(fid,'<head>\n');
  fprintf(fid,'<title>Non-redundant list of PDB structures from X-ray crystallography with resolution $le; %7.1fA</title>\n', MaxRes);
  fprintf(fid,'</head>\n');
  fprintf(fid,'<body>\n');
  fprintf(fid,'<h3>Non-redundant list of %d PDB structures from X-ray crystallography with resolution &le; %7.1fA and containing at least one RNA basepair</h3>\n', length(i), MaxRes);
  fprintf(fid,'<a href="index.html">Click here for the full list</a>\n');
end

fprintf(fid,'<table border=1>');

fprintf(fid,'<tr><th>PDBID</th><th>FR3D analysis</th><th>NTs</th><th>Pairs</th><th>Res</th><th>Date</th><th>Organism</th><th>Description provided by author</th><th>Author</th>');

switch Format,
case 1,
  fprintf(fid,'<th>Represented by</th></tr>');
case 2,
  fprintf(fid,'<th>Also represents</th></tr>');
end

for r = 1:length(tt(:,1)),
  fprintf(fid,'<tr>');
  fprintf(fid,'<td><a href="http://www.rcsb.org/pdb/cgi/explore.cgi?pdbId=%s">%s</a></td>', upper(tt{r,1}), upper(tt{r,1}));
  fprintf(fid,'<td><a href="%s/index.html">FR3D</a></td>', upper(tt{r,1}));

  if nn(r,1) == 0 || isnan(nn(r,1)),
    res = tt{r,3};
  else
    res = sprintf('%7.1f',nn(r,1));
  end
%  res = strrep(res,'ELECTRON MICROSCOPY','EM');
%  res = strrep(res,'SOLUTION NMR','NMR');

  fprintf(fid,'<td align="center">%d</td><td align="center">%d</td><td align="center">%s</td>', nn(r,2), nn(r,3), res);

  org = tt{r,8};
  if isempty(org),
    org = '&nbsp;';
  end

  fprintf(fid,'<td>%s</td><td>%s</td><td>%s</td><td>%s</td>',tt{r,4},org,tt{r,2},tt{r,5});

  switch Format,
  case 1,
    fprintf(fid,'<td>');
    if strcmpi(tt{r,10},tt{r,1}),
      fprintf(fid,'%s',tt{r,1});
    elseif ~isempty(tt{r,10}),
      fprintf(fid,'%s',tt{r,10});
    else
      fprintf(fid,'&nbsp;');
    end
    fprintf(fid,'</td>');
  case 2,
    g = find(ismember(t(:,10),tt{r,10}));         % find others in this group
    fprintf(fid,'<td width="150">');
    if length(g) > 1,
     for h = 1:length(g),
      if ~strcmpi(t{g(h),1},tt{r,1}),
        fprintf(fid,'%s ',t{g(h),1});
      end
     end
    else
     fprintf(fid,'&nbsp;');
    end
    fprintf(fid,'</td>');
  end

  fprintf(fid,'</tr>\n');
end

fprintf(fid,'</table>\n');
fprintf(fid,'</html>');


switch Format
case 1,
  fid = fopen(['PDBFiles' filesep 'Nonredundant_' NRDate '_All_list.pdb'],'w');
  for i = 1:length(tt(:,1)),
    if strcmp(tt{i,1},tt{i,10}),               % file represents itself
      fprintf(fid,'%s\n', tt{i,1});                % write to the list
    end
  end
  fclose(fid);

case 2,


end



end
