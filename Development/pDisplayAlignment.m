
  fprintf('JAR3D %d\n',jj);
  fprintf('%s\n', Alig.get(0));
  for i = 1:NumSeq,
    if 3*i <= size(Alig),
      fprintf('%s %30s %s\n', Alig.get(3*i-2), Alig.get(3*i-1), Alig.get(3*i));
    end
  end
