clear ; close all; clc

bestNb = 0;
bestLambda = 0;
bestF1 = 0;
load "bundlePCA"
for nb_hidden_layer = 70:20:170
  for lambda = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30]
    fileName = ["trained/", num2str(nb_hidden_layer), "-", num2str(lambda), "-4.mat"];
    [reussite, F1] = getScore(fileName);
    if (F1 > bestF1)
      bestF1 = F1;
      bestNb = nb_hidden_layer;
      bestLambda = lambda;
    end;
  end
end;

bestF1
bestNb
bestLambda
