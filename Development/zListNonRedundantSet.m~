% zListNonredundantSet lists PDB ID's and chains that are meant to be non-redundant

Verbose = 0;

NRList = 'Nonredundant_2009-05-14_list';

% Filenames = zReadPDBList(NRList,1);

% File = zAddNTData(NRList,0,[],1);

fid = fopen([pwd filesep 'Web' filesep 'AnalyzedStructures' filesep 'NonredundantList.txt'],'w');

Vers = num2str(File(1).ClassVersion);
fprintf(fid,'# PDB_ID_FR3D_Version_%s\n',Vers);

fprintf(fid,'PDB_ID\tChain(s)\n');

for f = 1:length(File);
  if length(File(f).NT) > 0,
    clear ChainInd ChainRed
    Chain = cat(2,File(f).NT.Chain);                  % chains
    U = unique(Chain);
    for u = 1:length(U),                              % loop through chains
      ChainInd{u} = find(Chain == U(u));              % NTs in this chain
    end
    for u = 1:length(U),
      for v = 1:length(U),
        ChainRed(u,v) = sum(sum(File(f).Redundant(ChainInd{u},ChainInd{v})));
      end
    end

    KeepChains = [];

    if sum(sum(ChainRed)) > 0,                        % redundant chains exist
      clear BPQual IntQual
      E = File(f).Edge;
      Q = (abs(E) > 0) .* (abs(E) < 13);              % indicator of edges
      R = (abs(E) > 0) .* (abs(E) < 25);              % indicator of edges
      for u = 1:length(U),
        for v = 1:length(U),
          BPQual(u,v) = sum(sum(Q(ChainInd{u},ChainInd{v})));
          IntQual(u,v) = sum(sum(R(ChainInd{u},ChainInd{v})));
        end
      end

      QualP = full((BPQual + BPQual^2 + BPQual^3 + BPQual^4 + BPQual^5) > 0);

      if Verbose > 1,
        clf
        U
        zCircularDiagram(File(f))
        drawnow
        full(ChainRed)
        full(IntQual)
        QualP
      end

      if all(QualP > 0),
        KeepChains = 1:length(U);
        if Verbose > 1,
          fprintf('Keeping all chains\n');
        end
      else
        G  = 1-(sum(QualP) > 0);             % chains used, or no inter
        m  = 0;                              % maximum quality so far
        n  = 0;
        while any(G == 0),                   % unused chains
          i = min(find(G == 0));              % first unused chain
          KC = find(QualP(i,:) > 0);          % chains with 1st chain
          G(KC) = ones(1,length(KC));          % these have been considered
          Q  = sum(sum(BPQual(KC,KC)));       % quality of these chains
          QInt = sum(sum(IntQual(KC,KC)));
          if Q > m,
            KeepChains = KC;
            m = Q;
            n = QInt;
          elseif Q == m,
            if QInt > n,
              KeepChains = KC;
              n = QInt;
            end
          end
        end
        if Verbose > 1,
          fprintf('Keeping chains ');
          fprintf('%d ', KeepChains);
          fprintf('%s ', U(KeepChains));
          fprintf('\n');
        end
      end

      if Verbose > 1,
        pause
      end
    else
      KeepChains = 1:length(U);                 % keep them all
    end

    fprintf(fid,'%s\t',File(f).Filename);
    fprintf(fid,'%s',U(KeepChains));
    fprintf(fid,'\n');

  end
end

fclose(fid);



















