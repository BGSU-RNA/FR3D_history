% zPairViewer loads one or more pdb files, allows the user to specify
% selection criteria, then displays pair information in various formats

clear ViewParam

if ~exist('File'),                  % If files haven't already been loaded
  zSpecifyPDBFiles;        % The master list
  File = zGetNTData(Filenames,0);
end

ViewParam.FigNum = 2;

while ViewParam.FigNum > 0,
  if exist('Param'),
    Param = zEnterPairSelection(Param);
  else
    Param = zEnterPairSelection([]);
  end
  SP    = zSelectPairs(File,Param);

  if length(SP) > 0,
    ViewParam.Mode   = 1;
    while ViewParam.Mode(1) > 0,
      ViewParam = zEnterViewMode(Param,ViewParam);

      if any(ViewParam.Mode == 1),
        SP = zColorPairs(File,SP,Param,ViewParam);
      end

      if ViewParam.Sort == 1,
        SP = zSortPairs(File,SP,ViewParam);
      else
        for i=1:length(ViewParam.Mode),
          switch ViewParam.Mode(i),
            case 1, FigsDone = zScatterPairs(File,SP,Param,ViewParam);
                    ViewParam.FigNum = ViewParam.FigNum + FigsDone;
            case 2, zListPairs(File,SP,2);
            case 3, zListPairs(File,SP,3);
            case 4, [File,SP,ViewParam]= zHandClassifyPairs(File,SP,ViewParam);
                    for f=1:length(File),
                      if File(f).Modified == 1,
                        zWriteHandFile(File(f));
                        File(f).Modified = 0;
                        zSaveNTData(File(f));
                      end
                    end
            case 5, zContextViewer(File,SP,Param,ViewParam);
          end
        end
      end
    end
  end
end  
