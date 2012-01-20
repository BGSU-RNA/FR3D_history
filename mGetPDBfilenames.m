%mGetPDBfilenames %By Ali Mokdad - April 14 2006
%This looks in the folder <PDBFiles> for files with extensions .pdb-1repetitions

% Extended to search Matlab's path for PDB files and to use the native
% file separator rather than \  - CLZ 2006-07-18

clear s temp

a = dir(['PDBFiles' filesep '*.pdb']);
a = [a; dir(['PrecomputedData' filesep '*.mat'])];

p = path;
c = [0 strfind(p,pathsep) length(p)+1];
for i=1:length(c)-1,
  a = [a; dir([p(c(i)+1:c(i+1)-1) filesep '*.pdb'])];
end

if ~isempty(a)
    for i=1:length(a)
      temp{i} =regexprep(a(i).name,'.pdb|.mat| ','');   % strip extensions and spaces
    end
    
    temp = sort(temp);

    s{1,1}=temp{1};
    count=2;
    for i=2:length(temp)
        if ~strcmp(temp{i},s{count-1,1})
            s{count,1}=temp{i};
            count=count+1;
        end
    end
else
    set(handles.Status,'String','ALERT: there are no PDB files or saved PDB data in the folders "PDBFiles" and "PrecomputedData". Please put some PDB files in these locations and try again!');
    s{1,1} =' ';%to prevent error
    set(handles.ReadQuery,'Visible','Off'); %just to prevent the user from clicking it anyway!
end

