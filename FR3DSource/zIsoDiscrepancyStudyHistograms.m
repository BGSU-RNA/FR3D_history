% zIsoDiscrepancyStudyWithin finds and histograms the isodiscrepancy between
% instances of each basepair family

if ~exist('File'),                           % if no molecule data is loaded,
  [File,SIndex] = zAddNTData('NonRedundant_2008_02_21_list',2);   % load PDB data
  File = File(SIndex);
else
  [File,SIndex] = zAddNTData('NonRedundant_2008_02_21_list',2,File); % add PDB data if needed
  File = File(SIndex);
end                       

% ------------------------------ Only consider high resolution files

Res = zeros(1,length(File));
for f = 1:length(File),
  if ~isempty(File(f).Info.Resolution),
    Res(f) = File(f).Info.Resolution;
  else
    Res(f) = 10;
  end
end
HFile = File(find((Res > 0) .* (Res < 3)));


StudyNum = -4;
Verbose = 2;
Decimal = 1;               % 1 means pay attention to subcategory

N = 200;                   % number of instances of each pair to keep

Param = [];

% 1-AA  2-CA  3-GA  4-UA  5-AC  6-CC  7-GC  8-UC 
% 9-AG 10-CG 11-GG 12-UG 13-AU 14-CU 15-GU 16-UU

Pairs = {'AA', 'CA', 'GA', 'UA', 'AC', 'CC', 'GC', 'UC', 'AG', 'CG', 'GG', 'UG', 'AU', 'CU', 'GU', 'UU'};

switch StudyNum,
case -4,
  CL = zClassLimits;
  pcodes = [1 5 6 7 9 11 13 14 15 16];    % pair codes to work on
  Param = [];
  for j = 1:length(pcodes),             % run through all pair codes specified
    pc = pcodes(j);
    CLE = CL(:,1,pc);                   % basepair category numbers
    CLE = CLE(find(CLE));               % leave out empty entries
    CLE = CLE(find(abs(CLE) < 13));     % basepairing only
    CLE = CLE(find(CLE == fix(CLE)));
    for row = 1:length(CLE),
      Param = [Param; [CLE(row) pc]];
    end
  end
end

List = [];
Study = '';

Start = [];
Start(1) = 1;                     % index of the first instance for Param 1

Text = ['Comparing ' Pairs{Param(1,2)} ' ' zEdgeText(Param(1,1)) ' with '];
Text = [Text Pairs{Param(2,2)} ' ' zEdgeText(Param(2,1))];
for i = 3:length(Param(:,1)),
  Text = [Text ', ' Pairs{Param(i,2)} ' ' zEdgeText(Param(i,1))];
end
fprintf('%s\n', Text);

Discreps = [];

for i = 1:length(Param(:,1)),
  List  = zFindPairs(HFile,Param(i,:),1);
  fprintf('Found %5d instances of %2s %4s class %5.1f\n', s, Pairs{Param(i,2)}, zEdgeText(Param(i,1),Decimal,Param(i,2)), Param(i,1));

  [s,t] = size(List);
  [y,j] = sort(rand(s,1));
  List = List(j,:);                     % order randomly
  List = List(1:min(s,N),:);            % take the first N instances

  [L,M] = size(List);

  ID = zeros(L,L);
  for k = 1:L,                   % Very slow nested loop
    f1    = List(k,3);
    Model = List(k,[1 2]);
    NT1 = HFile(f1).NT(Model(1));
    NT2 = HFile(f1).NT(Model(2));
    for m = (k+1):L,
      f2    = List(m,3);
      Cand  = List(m,[1 2]);
      NT3 = HFile(f2).NT(Cand(1));
      NT4 = HFile(f2).NT(Cand(2));

      [d,ang,t1,t2,cp,ga1,ga2,od] = zIsoDiscrepancy(NT1,NT2,NT3,NT4);

      ID(k,m)  = d;
    end
  end

  Discreps = [Discreps; nonzeros(ID)];        % accumulate non-zero values

  figure(1)
  subplot(1,2,1)
  hist(nonzeros(ID),30);                      % new pairs
  ax = axis;
  ax(2) = 5;
  axis(ax);
  title([Pairs{Param(i,2)} zEdgeText(Param(i,1),Decimal,Param(i,2))]);

  subplot(1,2,2)
  hist(Discreps,50);                          % all pairs
  ax = axis;
  ax(2) = 5;
  axis(ax);

  if max(nonzeros(ID)) > 5,
    drawnow
    z = sort(nonzeros(ID));
    fprintf('%s\n', [Pairs{Param(i,2)} zEdgeText(Param(i,1),Decimal,Param(i,2))]);
    z((end-10):end)
  end
end

figure(2)
clf
colormap('default')
hist(Discreps,75);                          % all pairs
ax = axis;
ax(1) = 0;
ax(2) = 5;
axis(ax);
title('IsoDiscrepancy within basepairs');
saveas(gcf,[strrep('IsoDiscrepancy within basepairs',' ','_') '.png'],'png');


q = [0.9 0.95 0.98 0.99 0.999 0.9999];
for i = 1:length(q),
  fprintf('%7.2f%% quantile is at %7.2f\n', 100*q(i), quantile(Discreps,q(i)));
end
