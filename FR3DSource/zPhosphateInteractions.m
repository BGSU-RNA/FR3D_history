% zPhosphateInteractions checks all nearby pairs of bases for base-phosphate
% interactions, and stores them in a sparse matrix field BasePhosphate


function [File,D] = zPhosphateInteractions(File,Verbose)

if nargin == 1,
  Verbose = 0;
end

D = [];                   % where to save data if Verbose

zStandardBases
Sugar = {'C1*','C2*','O2*','C3*','O3*','C4*','O4*','C5*','O5*','P','O1P','O2P','O3 of next'};

t = cputime;

for f = 1:length(File),

  if isempty(File(f).Distance),
    c = cat(1,File(f).NT(1:File(f).NumNT).Center); % nucleotide centers
    File(f).Distance = zMutualDistance(c,16); % compute distances < 16 Angs
  end

  File(f).BasePhosphate = sparse(zeros(File(f).NumNT));

  % -------- First screening of base pairs ----------------------------------- 

  DistCutoff = 16;                              % max distance for interaction
  [i,j] = find((File(f).Distance < DistCutoff).*(File(f).Distance > 0)); 
                                                % screen by C-C distance

  i = [i; (1:length(File(f).NT))'];             % allow self interactions
  j = [j; (1:length(File(f).NT))'];             % allow self interactions

  % -------- Screen and analyze pairs ----------------------------------------

  p   = [9 11 12 13];                           % rows of the phosphate oxygens
  pn  = {'O5*','O1P','O2P','O3*'};              % names of phosphate oxygens

  for k = 1:length(i),                          % loop through possible pairs

    N1 = File(f).NT(i(k));                      % nucleotide i information
    N2 = File(f).NT(j(k));                      % nucleotide j information

    switch N1.Code
      case 1,                         % Base A
              h   = [11 12 14 15];    % rows of the base hydrogens
              hn  = {'H2','H8','1H6','2H6'}; % names of the base hydrogens
              m   = [ 9  7  6  6];    % rows of the corresponding massive atoms
              e   = [ 1  4  2  3];    % code for location of the interaction
      case 2,                         % Base C
              h   = [10 11 12 13];
              hn  = {'H6','H5','1H4','2H4'}; % names of the base hydrogens
              m   = [ 7  8  6  6];
              e   = [ 9  8  6  5];
      case 3,                         % Base G
              h   = [12 13 15 16];
              hn  = {'H1','H8','1H2','2H2'}; % names of the base hydrogens
              m   = [ 4  7 11 11];
              e   = [13 14 10 11];
      case 4,                         % Base U
              h   = [ 9 11 12];
              hn  = {'H5','H3','H6'}; % names of the base hydrogens
              m   = [ 8  4  7];
              e   = [16 15 17];
    end

    %  Meanings of the codes:
    %  1  A C2-H2  interacts with oxygen of phosphate, called 2BP
    %  2  A N6-1H6 interacts with oxygen of phosphate, called 6BP
    %  3  A N6-2H6 interacts with oxygen of phosphate, called 7BP
    %  4  A C8-H8  interacts with oxygen of phosphate, called 0BP
    %  5  C N4-2H4 interacts with oxygen of phosphate, called 6BP
    %  6  C N4-1H4 interacts with oxygen of phosphate, called 7BP
    %  7  C N4-1H4 and C5-H5 interact with oxygen(s) of phosphate, called 8BP
    %  8  C C5-H5  interacts with oxygen of phosphate, called 9BP
    %  9  C C6-H6  interacts with oxygen of phosphate, called 0BP
    % 10  G N2-1H2 interacts with oxygen of phosphate, called 1BP
    % 11  G N2-2H2 interacts with oxygen of phosphate, called 3BP
    % 12  G N2-2H2 and N1-H1 interacts with oxygen(s) of phosphate, called 4BP
    % 13  G N1-H1  interacts with oxygen of phosphate, called 5BP
    % 14  G C8-H8  interacts with oxygen of phosphate, called 0BP
    % 15  U N3-H3  interacts with oxygen of phosphate, called 5BP
    % 16  U C5-H5  interacts with oxygen of phosphate, called 9BP
    % 17  U C6-H6  interacts with oxygen of phosphate, called 0BP

    dis = zDistance(N1.Fit(m,:), N2.Sugar(p,:)); % distances between mass & O's
    dis = dis .* (dis < 4.5);          % massive-oxygen pairs close enough

    g = [];                           % internal classification number

    for mm = 1:length(m),             % massive atom to consider
     pp = find(dis(mm,:));            % oxygens close enough to consider
     for n = 1:length(pp),            % loop through potential oxygens
      Angle(n)=zAngle(N1.Fit(m(mm),:),N1.Fit(h(mm),:),N2.Sugar(p(pp(n)),:));
                                      % base massive - hydrogen - oxygen angle
      Dist(n) = dis(mm,pp(n));        % distance
     end
     [u,v] = min(-Angle+60*Dist);     % order by quality of potential bond

     for n = 1:length(pp),            % loop through potential oxygens
      if Angle(n) > 100,              % good enough to be "near" base-phosph

        if ((Angle(n) > 150) && (Dist(n) < 3.6)) % looks good
          g = [g e(mm)];              % assign a non-near class.
        else
          g = e(mm) + 100;            % > 100 means "near"
        end

        if Verbose > 1,
          % store information for later display
          ox = (N2.Sugar(p(pp(n)),:)-N1.Fit(1,:)) * N1.Rot; % oxygen displ
          ph2= (N2.Sugar(10,:)-N1.Fit(1,:)) * N1.Rot; % phosphorus displacement

          a = [f i(k) j(k) N1.Code g(end) mm pp(n) Angle(n) Dist(n) ox ph2 File(f).Distance(i(k),j(k)) (v(1)==n)];

% [ g(end) v(1) == n]

          % Columns:
          % 1  file number
          % 2  index of base
          % 3  index of nucleotide using phosphate
          % 4  code of base
          % 5  classification number for this massive-oxygen pair
          % 6  which massive atom is interacting
          % 7  which oxygen is interacting
          % 8  angle of interaction, in degrees
          % 9  distance from massive atom to oxygen, in Angstroms
          %10  displacement of oxygen atom relative to C1' of base
          %13  displacement of phophorus atom relative to C1' of base
          %16  distance between centers of the two bases
          %17  1 if this is the best oxygen for this hydrogen, 0 otherwise

          D = [D; a];                  % append data to data matrix

          if Verbose > 3,
            fprintf('%6s base %s%5s %3s %3d phosphate %s%5s %13s length %6.2f angle %6.2f interaction %s\n', File(f).Filename, N1.Base, N1.Number, AtomNames{h(mm),N1.Code}, g(end), N2.Base, N2.Number, Sugar{p(pp(n))}, dis(mm,pp(n)), Angle, zEdgeText(File(f).Edge(i(k),j(k))));

g

          end
        end

      end

     end  % loop over potential oxygens

     if length(g) > 0,
       if (min(g) < 100) && (max(g) > 100),
         g = g(find(g < 100));               % remove near classifications
       end
       if length(g) > 1,                     % multiple bonds
         g = sort(g);
         if (g(1) == 6) && (g(end) == 8),
           g = 7;                            % two bonds formed
         elseif (g(1) == 11) && (g(end) == 13),
           g = 12;                           % two bonds formed
         else
           fprintf('Another case to consider\n');
           g
           [f i(k) j(k)]
         end
       end

       File(f).BasePhosphate(i(k),j(k)) =   g(1);  % record classification
     end

    end   % loop over massive atoms
  end     % loop over nucleotide pairs
end       % loop over files


if Verbose > 1,

  fprintf('Classifying base-phosphate interactions took %8.2f minutes\n', (cputime-t)/60);

  zPhosDisplay                     % temporary arrangement

end

return

% File = zAddNTData({'1s72','1j5e','2avy','2aw4','2j01'});
% zPhosphateInteractions(File,3);

for f = 1:length(File);
  if isempty(File(f).Info.Resolution),
    Res(f) = 10;
  else
    Res(f) = File(f).Info.Resolution;
  end
end
F = File(find(Res < 4.0));

