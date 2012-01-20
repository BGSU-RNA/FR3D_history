
load PDBInfo

for r = 1:length(t(:,1)),
  if isempty(t{r,10}),
    t{r,10} = '&nbsp;';
  end
end

Format = 2;

switch Format
case 1, 
  i = 1:length(t(:,1));
  fid = fopen('Web\AnalyzedStructures\index.html','w');
  [y,j] = sort(-n(i,2));                  % sort by number of nucleotides
  i = i(j);


case 2,
  OK = zeros(1,length(t(:,1)));
  for r = 1:length(t(:,1)),
    if strcmpi(t{r,1},t{r,10}),
      OK(1,r) = 1;
    end
  end

  i = find(OK);
  [y,j] = sort(-n(i,2));                  % sort by number of nucleotides
  i = i(j);

  fid = fopen('Web\AnalyzedStructures\nonredlist.html','w');
end


tt = t(i,:);
nn = n(i,:);


fprintf(fid,'<html>');

switch Format
case 1,
  fprintf(fid,'<head>\n');
  fprintf(fid,'<title>List of all RNA-containing structures</title>\n');
  fprintf(fid,'</head>\n');
  fprintf(fid,'<body>\n');
  fprintf(fid,'<a href="nonredlist.html">Click here for the reduced-redundancy list</a>\n');
case 2,
  fprintf(fid,'<head>\n');
  fprintf(fid,'<title>Reduced-redundancy list of RNA-containing structures from X-ray crystallography with resolution <= %7.1f</title>\n', max(nn(:,1)));
  fprintf(fid,'</head>\n');
  fprintf(fid,'<body>\n');
  fprintf(fid,'<a href="index.html">Click here for the full list</a>\n');
end

fprintf(fid,'<table border=1>');

fprintf(fid,'<tr><th>PDBID</th><th>FR3D analysis</th><th>NTs</th><th>Pairs</th><th>Res</th><th>Date</th><th>Organism</th><th>Description provided by author</th><th>Author</th>');

switch Format,
case 1,
  fprintf(fid,'<th>Represented by</th></tr>');
case 2,
  fprintf(fid,'<th>Represents</th></tr>');
end

for r = 1:length(tt(:,1)),
  fprintf(fid,'<tr>');
  fprintf(fid,'<td><a href="http://www.rcsb.org/pdb/cgi/explore.cgi?pdbId=%s">%s</a></td>', upper(tt{r,1}), upper(tt{r,1}));
  fprintf(fid,'<td><a href="%s/index.html">FR3D</a></td>', upper(tt{r,1}));

  res = nn(r,1);
  if res == 0,
    res = tt{r,3};
  else
    res = sprintf('%7.1f',res);
  end

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
      fprintf(fid,'&nbsp;');
    elseif ~isempty(tt{r,10}),
      fprintf(fid,'%s ',tt{r,10});
    else
      fprintf(fid,'&nbsp;');
    end
    fprintf(fid,'</td>');
  case 2,
    g = find(ismember(t(:,10),tt{r,1}));
    fprintf(fid,'<td width="150">');
    if length(g) > 1,
     for h = 1:length(g),
      if ~strcmpi(t{g(h),1},tt{r,1}),
        fprintf(fid,'%s ',t{g(h),1});
      end
     end
    else
     fprintf(fid,'Unique');
    end
    fprintf(fid,'</td>');
  end

  fprintf(fid,'</tr>\n');
end

fprintf(fid,'</table>\n');
fprintf(fid,'</html>');