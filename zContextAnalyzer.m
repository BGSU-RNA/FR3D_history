% zContextAnalyzer(File,SP,Param,ViewParam) displays information about
% bases and interactions in the neighborhood of a given type of pair, as
% specified using zContextExplorer

function [void] = zContextAnalyzer(File,SP,Param,ViewParam)

    clear Configs
    m = 1;
%    for k = 1:length(SP),                 % do all pairs
    for k = 1:min(200,length(SP)),         % only do the first 200
      a = SP(k).B1Index;
      b = SP(k).B2Index;
      switch Param.Context,
        case 1, C.Indices = [a b a+1 b-1];
        case 2, C.Indices = [a b a-1 b+1];
        case 3, C.Indices = [a b a-1 b+1 a+1 b-1];
      end
      C.Filenum = SP(k).Filenum;
      if (min(C.Indices) > 0) & (max(C.Indices) <= length(File(C.Filenum).NT)),
        C.Filename = File(C.Filenum).Filename;
        for i=1:length(C.Indices),
          C.NT(i) = File(C.Filenum).NT(C.Indices(i));
        end
        C.NumNT=length(C.Indices);
        Configs(m) = C;
        m = m + 1;
      end
    end

  if length(Configs) == 0,
    fprintf('No such configurations found\n');
  else

    for k = 1:length(Configs),
      fprintf('%10s:',Configs(k).Filename);
      for j = 1:length(Configs(k).Indices)/2,
        i1 = Configs(k).Indices(2*j-1);
        i2 = Configs(k).Indices(2*j);
        f  = Configs(k).Filenum;
        N1 = Configs(k).NT(2*j-1);
        N2 = Configs(k).NT(2*j);
        fprintf(' %1s%4s-%5.1f-%1s%4s  ', N1.Base, N1.Number, File(f).Inter(i1,i2), N2.Base, N2.Number);
      end
      fprintf('\n');
    end

    for k = 1:length(Configs),
      fprintf('%10s:',Configs(k).Filename);
      f  = Configs(k).Filenum;
      i  = Configs(k).Indices;
      L  = File(f).Inter(i,i);
      for j = 1:length(Configs(k).Indices),
        i1 = Configs(k).Indices(2*j-1);
        i2 = Configs(k).Indices(2*j);
        N1 = Configs(k).NT(2*j-1);
        N2 = Configs(k).NT(2*j);
        fprintf(' %1s%4s-%5.1f-%1s%4s  ', N1.Base, N1.Number, File(f).Inter(i1,i2), N2.Base, N2.Number);
      end
      fprintf('\n');
    end
  
    for k = 1:length(Configs),
      z(k) = char(Configs(k).NT(3).Base);
      zz(k) = char(Configs(k).NT(4).Base);
    end
    [Table, Chi2, P, Labels] = crosstab(z',zz');
    fprintf('Frequency count of bases stacked on this pair\n');
    zShowTable(Labels, Table);

    D = zeros(length(Configs));
    for k = 1:length(Configs),
      for j = (k+1):length(Configs),
        D(k,j) = sDiscrepancy(Configs(k).NT,Configs(j).NT);
      end
    end
    D = D + D';
    
    [y,i] = sort(sum(D));

    Configs = Configs(i);

    figure(1)

    ViewParam.LineStyle = '-';
    ViewParam.Sugar = 1;

    for i=1:length(Configs),
      clf
      f = Configs(i).Filenum;
      [R,S] = zPlotSomeNTAtOrigin(File(f),Configs(i).Indices,ViewParam);
      grid on
      view(ViewParam.az,ViewParam.el);

      pause
      [ViewParam.az,ViewParam.el] = view;
    end
  end
