% the function xDisplayCandidates(File,Model,Candidates) is used to display
% the Model and the possible matches to the model

function xDisplayCandidates(File,Search)

Model       = Search.Query;
Candidates  = Search.Candidates;
Discrepancy = Search.Discrepancy;

N = Model.NumNT;

warning off

% if there is no model, use the first candidate

if ~isfield(Model,'Filename'),
  f             = Candidates(1,N+1);
  Model.Indices = double(Candidates(1,1:N));
  Model.NT      = File(f).NT(Model.Indices);
  Model.Centers = cat(1,Model.NT(:).Center);
end

if Model.Geometric == 0,
  Model.WeightedCenter = mean(Model.Centers);
  Model.WeightedCenteredCenters = Model.Centers-ones(N,1)*mean(Model.Centers);
end

if isempty(Candidates)
    disp('There are no candidates to display')
    return
end

if Model.Sequential == 1,
  NowSequential = 1;
else
  NowSequential = 0;
end

L = length(Candidates(:,1));
Sequential=0;

Written = zeros(1,L);         % keep track of what was written already
Align   = zeros(1,L);

for i=1:2
  Display(i).n = 1;
  Display(i).sugar = 1;
  Display(i).neighborhood = 0;
  Display(i).superimpose  = 0;
  Display(i).supersugar   = 0;
end

PlotMotif(File,Model,Candidates,Discrepancy,Display,1,Written);  % sub function
PlotMotif(File,Model,Candidates,Discrepancy,Display,2,Written);

Numplots = 2;
stop     = 1;

while stop>0,                              % display Matches to the Model

    k=menu('Display options','Next plot','Previous Plot','Align Plots',...
           'Add plot','Toggle Sequential','Larger Neighborhood', ...
           'Restore original','Toggle sugar','Toggle superimpose', ...
           'Write Current','Quit');
    i=gcf;
    figure(1);
    [az,el]=view;
    x=XLim;
    y=YLim;
    z=ZLim;
    switch k
        case 1                                      % next plot
            figure(i);
            [az,el]=view;
            if i>1
                Display(i).n=min(L,Display(i).n+1);
                PlotMotif(File,Model,Candidates,Discrepancy,Display,i,Written);
            end
            view(az,el);

        case 2                                      % Previous Plot
            figure(i);
            [az,el]=view;
            if i>1
                Display(i).n=max(1,Display(i).n-1);
                PlotMotif(File,Model,Candidates,Discrepancy,Display,i,Written);
            end
            view(az,el);
        case 3                                      % Align Plot
            figure(1)
            % now find the camera angles
            [az,el]=view;
            x=XLim;
            y=YLim;
            z=ZLim;
            for i=2:Numplots
                %if sum(ismember(get(i,'Visible'),'n'))   %check to see if plot is there
                    figure(i)
                    % now change figure i's camera angles
                    view([az,el]);
                    set(gca,'XLim',x);
                    set(gca,'YLim',y);
                    set(gca,'ZLim',z);
                    %end
            end
        case 4                                      % Addplot
            figure(1)
            [az,el]=view;
            Numplots=Numplots+1;
            Display(Numplots).neighborhood=0;
            Display(Numplots).n=1;
            Display(Numplots).sugar=1;
            Display(Numplots).superimpose = 0;
            
            PlotMotif(File,Model,Candidates,Discrepancy,Display,Numplots,Written);
            view([az,el])
         case 5      % Toggle Sequential
            if Model.Sequential == 1,
              disp('Only sequential motifs have been retained in the search');
            elseif NowSequential == 0,
              OrigCandidates = Candidates;
              Candidates = xSequential(Candidates,N,3);
              if length(Candidates(:,1)) > 0,
                NowSequential = 1;
                figure(i)
                Display(i).n = 1;
                PlotMotif(File,Model,Candidates,Discrepancy,Display,i,Written);
              else
                disp('No sequential motifs were found');
                Candidates = OrigCandidates;
              end
            else
              Candidates = OrigCandidates;
              NowSequential = 0;
              figure(i)
              Display(i).n = 1;
              PlotMotif(File,Model,Candidates,Discrepancy,Display,i,Written);
            end

        case 6                                      % toggle Neighborhood
            dn = [4 4 4 4 6 6 8 8 0];
            Display(i).neighborhood = dn(1+Display(i).neighborhood);
            PlotMotif(File,Model,Candidates,Discrepancy,Display,i,Written);
            view(az,el);
        case 7                                      % Restore original
            Display(i).neighborhood=0;
            Display(i).sugar=1;
            PlotMotif(File,Model,Candidates,Discrepancy,Display,1,Written);
        case 8                                      % toggle sugar
            figure(i)
            [az,el]=view;

            if Display(i).superimpose == 0,
              Display(i).sugar = 1 - Display(i).sugar;
            elseif (Display(i).sugar == 0) & (Display(i).supersugar == 0),
              Display(i).sugar = 1;
            elseif (Display(i).sugar == 1) & (Display(i).supersugar == 0),
              Display(i).supersugar = 1;
            elseif (Display(i).sugar == 1) & (Display(i).supersugar == 1),
              Display(i).sugar = 0;
            elseif (Display(i).sugar == 0) & (Display(i).supersugar == 1),
              Display(i).supersugar = 0;
            end

            PlotMotif(File,Model,Candidates,Discrepancy,Display,i,Written);
            view(az,el);
        case 9                                      % toggle superimpose
            Display(i).superimpose = 1-Display(i).superimpose;
            PlotMotif(File,Model,Candidates,Discrepancy,Display,i,Written);
            view(az,el);
        case 10                            % Print Current - Added by Ali
            if ~exist('fidOUT','var')
              OUTfile = ['Model_' Model.Name 'Selected.txt'];
              fidOUT = fopen(OUTfile,'w+');
                                                % write header
              fprintf(fidOUT,'%s\t%s','PDB File','Discrepancy');
              for b=1:Model.NumNT,
                fprintf(fidOUT,'\t%s',strcat('NT',num2str(b)));
              end
              fprintf(fidOUT,'\n');
            end

            if Written(Display(i).n) < 1,
              ccc = double(Candidates(Display(i).n,:));
              f   = ccc(Model.NumNT+1);
              fprintf(fidOUT,'%s\t%8.4f',File(f).Filename,Discrepancy(Display(i).n));
              for b=1:Model.NumNT,
                fprintf(fidOUT,'\t%s%4s',File(f).NT(ccc(b)).Base,File(f).NT(ccc(b)).Number);
              end
              fprintf(fidOUT,'\n');
              Written(Display(i).n) = 1;
            end
            for j=2:Numplots
              PlotMotif(File,Model,Candidates,Discrepancy,Display,j,Written);
            end
            figure(i);
            [az,el]=view;
            if i>1
              Display(i).n=min(L,Display(i).n+1);
              PlotMotif(File,Model,Candidates,Discrepancy,Display,i,Written);
            end
            view(az,el);
        case 11                                      % quit Display
            if exist('fidOUT','var')
               fclose(fidOUT);
            end
            stop=0;
            for i=1:Numplots
%                set(figure(i),'Visible','off')
            end
    end
    set(gca,'XLim',x);
    set(gca,'YLim',y);
    set(gca,'ZLim',z);
end  % end while

%-------------------------------------------------------------------------

function  PlotMotif(File,Model,Candidates,Discrepancy,Display,i,Written)

    N = Model.NumNT;

    figure(i)
    clf

    if (i > 1) & (Display(i).superimpose == 1),
      Indices = Model.Indices;
      if (Model.NumNT > 2),
        R = eye(3);
        S = mean(cat(1,Model.NT.Center));
      else
        R = Model.NT(1).Rot;
        S = mean(cat(1,Model.NT.Center));
      end
      MVP.Rotation = R;
      MVP.Shift    = S;
      MVP.LineStyle = '-.';
      MVP.LineThickness = '1';
      MVP.Sugar    = Display(i).supersugar;
      MVP.ConnectSugar = 0;
      zDisplayNT(Model,1:N,MVP);
    end

    if (i == 1),
      Indices = 1:N;
      zShowInteractionTable(Model,1:N);
    else
      [s,t] = size(Candidates);
      n       = Display(i).n;
      f       = Candidates(n,N+1);
      Disc    = Discrepancy(n);
      Indices = double(Candidates(n,1:N));
      zShowInteractionTable(File(f),double(Indices));
    end

    [az,el]=view;
    VP.Sugar    = Display(i).sugar;

    if (i==1) & (Model.NumNT > 2),
      R = eye(3);
      S = Model.WeightedCenter;
    elseif (i==1) & (Model.NumNT == 2),
      R = Model.NT(1).Rot;
      S = mean(cat(1,Model.NT(Indices).Center));
    elseif Model.NumNT > 2,
      MC = Model.WeightedCenteredCenters;          % align to the model
      CandiCenters = cat(1,File(f).NT(Indices).Center);
      CC = CandiCenters - ones(N,1)*mean(CandiCenters);

      R = zBestRotation(MC, CC);
      S = mean(CandiCenters);
    else
      R = File(f).NT(Indices(1)).Rot;
      S = mean(cat(1,File(f).NT(Indices).Center));
    end

    VP.Rotation = R;
    VP.Shift    = S;

    if Display(i).neighborhood > 0,
      a = zeros(1,File(f).NumNT);
      v = Display(i).neighborhood;
      for j=1:length(Indices),
        a = a + (File(f).Distance(Indices(j),:) < v) .* ...
                (File(f).Distance(Indices(j),:) > 0);
      end
      a(Indices) = zeros(1,length(Indices));  % take out ones in Indices
      B = find(a);
      Indices = [Indices B];
    end

    if i==1,
      VP.ConnectSugar = 0;
      zDisplayNT(Model,Indices,VP);
    else
      zDisplayNT(File(f),Indices,VP);
    end

    if i==1
      if Model.Geometric == 1,
        xlabel('Query motif');
      else
        xlabel('First candidate');
      end
    else
      xlabel(['Plot ',int2str(n),' of ',int2str(s),'   Disc: ',...
             num2str(Disc)]);
      if Written(n) == 1;
        ylabel('Written to file');
      else
        ylabel('');
      end
    end
    grid on
    set(gcf,'Renderer','OpenGL')
%    set(gcf,'Renderer','zbuffer')


