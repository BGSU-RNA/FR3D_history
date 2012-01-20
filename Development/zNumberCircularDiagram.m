
function [A] = zNumberCircularDiagram(File,View,Thickness,r)

N = length(File.NT);

if r >= 1,
  ra = 1.02 * r;
  rb = 1.04 * r;
  rc = 1.11 * r;
else
  ra = r/1.02;
  rb = r/1.04;
  rc = r/1.2;
end

% ------------------------------------- Determine where to plot each NT

A = zeros(1,N);           % angle around circle for nucleotides
A(1) = 1;

if View(8) == 1,
  spaces = [1 4 8 12];
else
  spaces = [1 0 0 12];
end

nc = [1];                                % new chain starting
for t = 1:(N-1),
  A(t+1) = A(t) + spaces(1);
  if File.Covalent(t,t+1) == 0,
    A(t+1) = A(t+1) + spaces(2);
  end
  if (File.NT(t).Chain ~= File.NT(t+1).Chain),
    A(t+1) = A(t+1) + spaces(3);
    nc = [nc t+1];
  end
end
A(end+1) = A(end) + spaces(4);
mA = max(A);
d = find(diff(A) > 1);                  % locations of jumps in A
d = [0 d];                              % prepend for the first nucleotide

A = A * (-2*pi) / mA;                   % convert A to radians

% ---------------------------------- Draw the outside circle

for m = 1:(length(d)-1),
  theta = [A(d(m)+1):-0.01:A(d(m+1)) A(d(m+1))];
  plot(r*cos(theta),r*sin(theta),'k');
  hold on
end

for n = 1:length(nc),
  theta = A(nc(n));

  if cos(theta) > 0,
    angle = 180*theta/pi;
    if r >= 1,
      ha    = 'left';
    else
      ha    = 'right';
    end
    va    = 'top';
  else
    angle = 180*(theta - pi)/pi;
    if r >= 1,
      ha    = 'right';
    else
      ha    = 'left';
    end
    va    = 'bottom';
  end

  text(rc*cos(theta), rc*sin(theta), ['Chain ' File.NT(nc(n)).Chain], 'Rotation', angle, 'HorizontalAlignment', ha, 'VerticalAlignment', va, 'FontSize', 6);

end

if N > 1000,
  step  = 50;
  sstep = 10;
elseif N > 500,
  step = 20;
  sstep = 5;
elseif N > 300,
  step = 10;
  sstep = 2;
elseif N > 100,
  step = 5;
  sstep = 1;
else
  step = 1;
  sstep = 1;
end

kk = 1;

flag = 0;
for k = 1:N,
  kkk = str2num(File.NT(k).Number);
  if ~isempty(kkk),                               % it's really a number
    kk = kkk;
    flag = 0;                                     % OK to use next w/ letter
  else                                            % it has a letter in it
    kkk = str2num(File.NT(k).Number(1:(end-1)));  % omit last character
    if ~isempty(kkk),
      kk = kkk;                                   % use this for display
      flag = 1;                                   % but only once
    end
  end
  theta = A(k);
  if cos(theta) > 0,
    angle = 180*theta/pi;
    if r >= 1,
      ha    = 'left';
    else
      ha    = 'right';
    end
  else
    angle = 180*(theta - pi)/pi;
    if r >= 1,
      ha    = 'right';
    else
      ha    = 'left';
    end
  end

  if (mod(kk,sstep) == 0) && (mod(kk,step) ~= 0) && (Thickness < 1) && (flag == 0),
    plot([r*cos(theta) ra*cos(theta)], [r*sin(theta), ra*sin(theta)],'k');
    text(rb*cos(theta), rb*sin(theta), File.NT(k).Number,'FontSize',4, 'Rotation', angle, 'HorizontalAlignment', ha, 'VerticalAlignment', 'middle');
  end
  if mod(kk,step) == 0 && (flag == 0),
    plot([r*cos(theta) rb*cos(theta)], [r*sin(theta), rb*sin(theta)],'k');
    text(rc*cos(theta), rc*sin(theta), File.NT(k).Number, 'Rotation', angle, 'HorizontalAlignment', ha, 'VerticalAlignment', 'middle');
  end
end