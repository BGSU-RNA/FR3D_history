
% xCandidateMap

function [void] = xCandidateMap(File,Search)

C = Search.Candidates;
N = length(C(1,:));
f = find(C(:,N)==1);       % find the candidates from file 1

for i=1:length(f)
hold on
        V=cat(1,File(1).NT(C(i,1:N-1)).Center);
        scatter3(V(:,1),V(:,2),V(:,3),28,[0 0 1],'filled');      %plots the points in 3 space
        for j=1:N-2
            for k=j+1:N-1
              plot3([V(j,1),V(k,1)],[V(j,2),V(k,2)],[V(j,3),V(k,3)],'b');    %draws the lines between dots
            end
        end
end
hold off