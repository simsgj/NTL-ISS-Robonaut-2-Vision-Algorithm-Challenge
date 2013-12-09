function [reussite, F1, prec, rec] = getScore(fileName)

load "bundlePCA"
str = sprintf("load -binary %s",fileName);
eval(str);
fraction = 0.7;
X = X((round(size(X, 1) * 0.6) + 1):(round(size(X, 1) * 0.8)), :);
y = y((round(length(y) * 0.6) + 1):(round(length(y) * 0.8)));
m = size(X, 1);
[vals, ids, idx] = unique(y);

rp = randperm(m);

score = 0;

predictions = [];
for i = 1:m
  
  
  pred = predict(Theta1, Theta2, X(rp(i),:));
  score += vals(pred) == y(rp(i));
  predictions(rp(i)) = vals(pred);

end

predictions = predictions';

reussite = score / m;


    % Number of true positives
    tp = sum(predictions != 0 & y != 0);

    % Number of false positives
    fp = sum(predictions != 0 & y == 0);
    
    % Number of false negatives
    fn = sum(predictions == 0 & y != 0);
    
    % Precision
    prec = tp / (tp + fp);
    
    % Recall
    rec = tp / (tp + fn);
    
    % F1 score
    F1 = (2 * prec * rec) / (prec + rec);

end;
