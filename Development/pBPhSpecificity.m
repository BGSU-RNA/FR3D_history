
function [Q] = pBPhSpecificity(bph,bp)

Code = [1 1 1 1 2 2 2 2 2 3 3 3 3 3 4 4 4 2 3]; % which base makes which BPh

Q = 0.10*ones(4,1);
Q(Code(bph),1) = 0.7;                   % 70% conservation of the base
                                        % is about the lowest observed in any
                                        % combination of circumstances

if bp == 1,                             % the base is part of a pair

  % more specific calculations can go here

else


end
