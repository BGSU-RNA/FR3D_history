%mGetPDBfilenames %By Ali Mokdad - April 14 2006
%This looks in the folder <PDBFiles> for files with extensions .pdb
%and also looks in the folder <PrecomputedData> for files with extensions .mat
%Then it removes all extensions and makes a unified sorted list with no repetitions


x=ls('PDBFiles/*.pdb*'); %Beta, this includes *.pdbs
y=ls('PrecomputedData/*.mat');

z=strvcat(x,y);
z=sortrows(z);

if ~isempty(z)
    for i=1:length(z(:,1))
        temp{i}=regexprep(z(i,:),'.pdb*|.mat| ','');
    end

    s{1}=temp{1};
    count=2;
    for i=2:length(temp)
        if ~strcmp(temp{i},s{count-1})
            s{count}=temp{i};
            count=count+1;
        end
    end
else
    set(handles.Status,'String','ALERT: there are no PDB files or saved PDB data in the folders "PDBFiles" and "PrecomputedData". Please put some PDB files in these locations and try again!');
    s=' ';%to prevent error
    set(handles.ReadQuery,'Visible','Off'); %just to prevent the user from clicking it anyway!
end

