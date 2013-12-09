nb_iter = 500;
nb_outlaws = 0;

load powerSwitchUpData;
for k_iter = 1:nb_iter
  clear X y;


  rp = randperm(size(XRaw, 1));
  XRaw = XRaw(rp, :);

  m = round(size(XRaw(:, 1)) * 0.7);
  XRawTest = XRaw(m:end, :);

  y(:, 1) = (XRaw(:, 1) - XRaw(:, 3));
  y(:, 2) = (XRaw(:, 2) - XRaw(:, 4));

  ones = XRaw(:, 1);
  ones(:, 1) = 1;
  length1 = sqrt((XRaw(:, 5) - XRaw(:, 3)) .* (XRaw(:, 5) - XRaw(:, 3)) + (XRaw(:, 6) - XRaw(:, 4)) .* (XRaw(:, 6) - XRaw(:, 4)));
  length2 = sqrt((XRaw(:, 7) - XRaw(:, 3)) .* (XRaw(:, 7) - XRaw(:, 3)) + (XRaw(:, 8) - XRaw(:, 4)) .* (XRaw(:, 8) - XRaw(:, 4)));
  angle1 = atan2((XRaw(:, 6) - XRaw(:, 4)), (XRaw(:, 5) - XRaw(:, 3)));
  angle2 = atan2((XRaw(:, 8) - XRaw(:, 4)), (XRaw(:, 7) - XRaw(:, 3)));

  firstValues = [ones, length1, length2];
  secondValues = [ones, sin(angle1), cos(angle1), cos(angle2), sin(angle2)];

  firstCValues = cellstr(["1"; "length1"; "length2"]);
  secondCValues = cellstr(["1"; "sin(angle1)"; "cos(angle1)"; "cos(angle2)"; "sin(angle2)"]);

  num = 0;
  for i = firstValues
    for j = secondValues
      num += 1;
      X(:, num) = i .* j;
    end;
  end;

  XTest = X(m:end, :);
  yTest = y(m:end, :);

  X = X(1:m, :);
  y = y(1:m, :);

  theta = normalEqn(X, y);

  ypred = XTest * theta;
  error = calculateError(XTest, theta, y);
  distances = sqrt((ypred(:, 1) + XRawTest(:, 3) - XRawTest(:, 1)) .* (ypred(:, 1) + XRawTest(:, 3) - XRawTest(:, 1)) + (ypred(:, 2) + XRawTest(:, 4) - XRawTest(:, 2)) .* (ypred(:, 2) + XRawTest(:, 4) - XRawTest(:, 2)));


  nb_outlaws += sum(distances > 15) / size(XRawTest, 1);
end;
nb_outlaws = nb_outlaws / nb_iter;
nb_outlaws

%if (false)


file_id = fopen(strcat('c/code.c'), 'w');
fprintf(file_id, "float v1[%d];\n", length(firstCValues));
fprintf(file_id, "float v2[%d];\n", length(secondCValues));

numI = 0;
for i = 1:length(firstCValues)
  fprintf(file_id, "v1[%d] = %s;\n", numI, firstCValues{i});
  numI += 1;
end;

numJ = 0;
for j = 1:length(secondCValues)
  fprintf(file_id, "v2[%d] = %s;\n", numJ, secondCValues{j});
  numJ += 1;
end;

fprintf(file_id, "Point p;\n");
fprintf(file_id, "p.x = ");
num = 0;
for i = 1:length(firstCValues)
  for j = 1:length(secondCValues)
    num += 1;
    fprintf(file_id, "%f * v1[%d] * v2[%d]", theta(num, 1), i - 1, j - 1);
    if (i != length(firstCValues) || j != length(secondCValues))
      fprintf(file_id, " + ");
    end
  end;
end;
fprintf(file_id, ";\n");

fprintf(file_id, "p.y = ");
num = 0;
for i = 1:length(firstCValues)
  for j = 1:length(secondCValues)
    num += 1;
    fprintf(file_id, "%f * v1[%d] * v2[%d]", theta(num, 2), i - 1, j - 1);
    if (i != length(firstCValues) || j != length(secondCValues))
      fprintf(file_id, " + ");
    end
  end;
end;
fprintf(file_id, ";\n");
fprintf(file_id, "return p;");


fclose(file_id);
%end;
