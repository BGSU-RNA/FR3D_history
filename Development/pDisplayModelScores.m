
clear S

pModelScores;                            % read model scores copied from JAR3D

[s,t] = size(S);

for i = 1:s,
  S(i,:) = S(i,:) - max(S(i,:));         % subtract the maximum; make it zero
end

pcolor(S);
shading flat
caxis([-5 1]);
