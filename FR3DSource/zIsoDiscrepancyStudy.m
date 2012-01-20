% zIsoDiscrepancyStudy compares instances of two or more basepair families, making a variety of graphs 

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

StudyNumbers = [-4 1 5];
Verbose = 0;
Decimal = 1;              % do not include subcategories

N = 200;                   % number of instances of each pair to keep

Param = [];

% 1-AA  2-CA  3-GA  4-UA  5-AC  6-CC  7-GC  8-UC 
% 9-AG 10-CG 11-GG 12-UG 13-AU 14-CU 15-GU 16-UU

Pairs = {'AA', 'CA', 'GA', 'UA', 'AC', 'CC', 'GC', 'UC', 'AG', 'CG', 'GG', 'UG', 'AU', 'CU', 'GU', 'UU'};

for sn = 1:length(StudyNumbers),

StudyNum = StudyNumbers(sn);

switch StudyNum,
case -8,
  Param = [4 4];              % tWH UA versus tWH GG
  Param = [Param; [4 11]];
case -7,
  Param = [-4 13];              % tWH UA versus tWH CA
  Param = [Param; [-4 5]];
case -6,
  Param = [1 6];              % cWW CC versus cWW UU
  Param = [Param; [1 16]];
case -5,
  Param = [1 14];              % cWW CU versus cWW UU
  Param = [Param; [1 16]];
case -4,
  Param = [1 7];               % cWW GC versus cWW UA
  Param = [Param; [-1 4]];
case -3,
  Param = [2 13];
  Param = [Param; [-2 4]];
case -2,
  Param = [1 7];
  Param = [Param; [2 7]];
case -1,
  Param = [1 1];
  Param = [Param; [1 6]];
case 0,
  Param = [1 7];
  Param = [Param; [3 11]];
case 1,
  Param = [1 7];               % cWW GC versus cWW GU
  Param = [Param; [1 15]];
case 2,
  Param = [1 13];              % cWW AU versus cWS AU
  Param = [Param; [5 13]];
case 3,
  Param = [1 13];              % cWW AC versus AU
  Param = [Param; [1 5]];
case 4,
  Param = [1 9];               % cWW AG versus cWW AU
  Param = [Param; [1 13]];
case 5,
  Param = [1 15];              % cWW GU versus cWW UG
  Param = [Param; [-1 12]];
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

NText = ['IsoDiscrepancy between ' Pairs{Param(1,2)} ' ' zEdgeText(Param(1,1)) ' and ' Pairs{Param(2,2)} ' ' zEdgeText(Param(2,1))];

for i = 1:length(Param(:,1)),
  Lis  = zFindPairs(HFile,Param(i,:),1);
  [s,t] = size(Lis);
  [y,j] = sort(rand(s,1));
  Lis = Lis(j,:);                     % order randomly

  Start(i+1) = Start(i)+min(s,N);
  fprintf('Found %5d instances of %2s %4s class %5.1f\n', s, Pairs{Param(i,2)}, zEdgeText(Param(i,1),Decimal,Param(i,2)), Param(i,1));
  List = [List; Lis(1:min(s,N),:)];
  Study = [Study ' ' Pairs{Param(i,2)} zEdgeText(Param(i,1),Decimal,Param(i,2))];
end

[L,M] = size(List);

% xDisplayCandidates(HFile,List);

    PD = zeros(L,L);               % Initialize the pair discrepancy
    ID = zeros(L,L);
    ND = zeros(L,L);
    Ang = zeros(L,L);
    GAng = zeros(L,L);
    T1  = zeros(L,L);
    T2  = zeros(L,L);
    CP  = zeros(L,L);
    GA1  = zeros(L,L);
    GA2  = zeros(L,L);
    for k = 1:L,                   % Very slow nested loop
      f1    = List(k,3);
      Model = List(k,[1 2]);
      NT1 = HFile(f1).NT(Model(1));
      NT2 = HFile(f1).NT(Model(2));
      for m = (k+1):L,
        f2    = List(m,3);
        Cand  = List(m,[1 2]);
        PD(k,m) = xDiscrepancy(HFile(f1),Model,HFile(f2),Cand);

        NT3 = HFile(f2).NT(Cand(1));
        NT4 = HFile(f2).NT(Cand(2));

        [d,ang,t1,t2,cp,ga1,ga2,od] = zIsoDiscrepancy(NT1,NT2,NT3,NT4);
        ND(k,m)  = d;
        Ang(k,m) = ang;
        T1(k,m)  = sqrt(t1(1)^2 + t1(2)^2);
        T2(k,m)  = sqrt(t2(1)^2 + t2(2)^2);
        CP(k,m)  = cp;
        GA1(k,m) = ga1;         % difference between glyc bond angles
        GA2(k,m) = ga2;         % difference between glyc bond angles
        GAng(k,m) = sqrt(ga1^2 + ga2^2)/2;
        ID(k,m)  = od;
      end
    end

%    PD = sqrt(PD)/2;            % finish discrepancy calculation

    PD = PD + PD';
    ID = ID + ID';
    ND = ND + ND';

    if (length(List(:,1)) > 2) && (Verbose > 1),
      figure(1)
      clf
      q = zClusterGraph(ND,{},1,1:L);
      colormap('default');
      map = colormap;
      map = map((end-8):-1:8,:);
      colormap(map);
      caxis([0 8]);
      colorbar('location','eastoutside');

      title(['IsoDiscrepancy between instances of ' Study]);
      FN = ['Instance_NewIsoDiscrepancies_' Study];
      saveas(gcf,['Exemplars' filesep FN '.png'],'png');

      drawnow

      figure(2)
      clf
      q = zClusterGraph(ID,{},1,1:L);
      colormap('default');
      map = colormap;
      map = map((end-8):-1:8,:);
      colormap(map);
      caxis([0 16]);
      colorbar('location','eastoutside');

      title(['Old IsoDiscrepancy between instances of ' Study]);
      FN = ['Instance_IsoDiscrepancies_' Study];
      saveas(gcf,['Exemplars' filesep FN '.png'],'png');

      drawnow
    end

    if (length(List(:,1)) > 2) && (Verbose > 0),
      figure(3)
      clf
      subplot(3,3,4)
      hist(nonzeros(PD(1:Start(2)-1,Start(2):Start(3)-1)),30);
      title('Discrepancy between groups 1 and 2');
      ax = axis;
      ax(1) = 0;
      axis(ax);

      subplot(3,3,1)
      hist(nonzeros(PD(1:Start(2)-1,1:Start(2)-1)),30);
      title('Discrepancy within group 1');
      axis(ax);

      subplot(3,3,7)
      hist(nonzeros(PD(Start(2):Start(3)-1,Start(2):Start(3)-1)),30);
      title('Discrepancy within group 2');
      axis(ax);

      subplot(3,3,5)
      hist(nonzeros(ID(1:Start(2)-1,Start(2):Start(3)-1)),30);
      title('Old IsoDiscrepancy between groups 1 and 2');
      ax = axis;
      ax(1) = 0;
      axis(ax);

      subplot(3,3,2)
      hist(nonzeros(ID(1:Start(2)-1,1:Start(2)-1)),30);
      title('Old IsoDiscrepancy within group 1');
      axis(ax);

      subplot(3,3,8)
      hist(nonzeros(ID(Start(2):Start(3)-1,Start(2):Start(3)-1)),30);
      title('Old IsoDiscrepancy within group 2');
      axis(ax);

      subplot(3,3,6)
      hist(nonzeros(ND(1:Start(2)-1,Start(2):Start(3)-1)),30);
      title('Discrepancy between groups 1 and 2');
      ax = axis;
      ax(1) = 0;
      axis(ax);

      subplot(3,3,3)
      hist(nonzeros(ND(1:Start(2)-1,1:Start(2)-1)),30);
      title('Discrepancy within group 1');
      axis(ax);

      subplot(3,3,9)
      hist(nonzeros(ND(Start(2):Start(3)-1,Start(2):Start(3)-1)),30);
      title('Discrepancy within group 2');
      axis(ax);

      figure(4)
      clf

      g1 = 1:(Start(2)-1);
      g2 = g1;

      s = length(g1)*length(g2);
      a = reshape(Ang(g1,g2),1,s);
      b = reshape(sqrt(T1(g1,g2).^2+T2(g1,g2).^2)/2,1,s);
      c = reshape(CP(g1,g2),1,s);
      scatter3(a,b,c,4,'r','filled');

      hold on

      g1 = 1:(Start(2)-1);
      g2 = Start(2):(Start(3)-1);

      s = length(g1)*length(g2);
      a = reshape(Ang(g1,g2),1,s);
      b = reshape(sqrt(T1(g1,g2).^2+T2(g1,g2).^2)/2,1,s);
      c = reshape(CP(g1,g2),1,s);
      scatter3(a,b,c,4,'b','filled');

      g2 = Start(2):(Start(3)-1);
      g1 = g2;

      s = length(g1)*length(g2);
      a = reshape(Ang(g1,g2),1,s);
      b = reshape(sqrt(T1(g1,g2).^2+T2(g1,g2).^2)/2,1,s);
      c = reshape(CP(g1,g2),1,s);
      scatter3(a,b,c,4,'g','filled');
      ax = axis;
      ax(2) = pi;
      axis(ax);

      title(Text);
      xlabel('Angle of rotation');
      ylabel('Length of shift');
      zlabel('Difference in C1 distance');

      rotate3d on
      view(2)

      figure(5)
      clf

      g1 = 1:(Start(2)-1);
      g2 = g1;

      s = length(g1)*length(g2);
      a = reshape(GAng(g1,g2),1,s);
      b = reshape(sqrt(T1(g1,g2).^2+T2(g1,g2).^2)/2,1,s);
      c = reshape(CP(g1,g2),1,s);
      scatter3(a,b,c,4,'r','filled');

      hold on

      g1 = 1:(Start(2)-1);
      g2 = Start(2):(Start(3)-1);

      s = length(g1)*length(g2);
      a = reshape(GAng(g1,g2),1,s);
      b = reshape(sqrt(T1(g1,g2).^2+T2(g1,g2).^2)/2,1,s);
      c = reshape(CP(g1,g2),1,s);
      scatter3(a,b,c,4,'b','filled');

      g2 = Start(2):(Start(3)-1);
      g1 = g2;

      s = length(g1)*length(g2);
      a = reshape(GAng(g1,g2),1,s);
      b = reshape(sqrt(T1(g1,g2).^2+T2(g1,g2).^2)/2,1,s);
      c = reshape(CP(g1,g2),1,s);
      scatter3(a,b,c,4,'g','filled');
      ax = axis;
      ax(2) = pi;
      axis(ax);

      title(Text);
      xlabel('Glycosidic angle of rotation');
      ylabel('Length of shift');
      zlabel('Difference in C1 distance');

      rotate3d on
      view(2)

    end

      figure(6)
subplot(2,2,sn+1);
      NBin = hist(nonzeros(ND(1:Start(2)-1,Start(2):Start(3)-1)),50);
      hist(nonzeros(ND(1:Start(2)-1,Start(2):Start(3)-1)),50);
      title(NText);
      axis([0 6 0 1.1*max(NBin)]);
      saveas(gcf,[strrep(NText,' ','_') '.png'],'png');

q = [0.9 0.95 0.96 0.97 0.98 0.99 0.999 0.9999];
for i = 1:length(q),
  fprintf('%7.2f%% quantile is at %7.2f\n', 100*q(i), quantile(nonzeros(ND(1:Start(2)-1,Start(2):Start(3)-1)),q(i)));
end


if sn < length(StudyNumbers),
  disp('Press any key');
%  pause
end

end
 