clear ; close all; clc
fraction = 0.6;
load "bundle"
XDisplay = X;
yDisplay = y;
XDisplay = X((round(size(X, 1) * fraction) + 1):end, :);
yDisplay = y((round(length(y) * fraction) + 1):end);
load "bundlePCA"
load "best/80-0.01-4.mat"

[vals, ids, idx] = unique(y);
X = X((round(size(X, 1) * fraction) + 1):end, :);
y = y((round(length(y) * fraction) + 1):end);
m = size(X, 1);


rp = randperm(m);

score = 0;

predictions = [];
for i = 1:m
  
  
  pred = predict(Theta1, Theta2, X(rp(i),:));
  score += vals(pred) == y(rp(i));
  predictions(rp(i)) = vals(pred);
  %if (false)%(y(rp(i)) == 3) %(false) %(vals(pred) != y(rp(i)))
  displayData(XDisplay(rp(i), :));
  fprintf('\nNeural Network Prediction: expected: %d result: %d (correct: %d, %d)\n', y(rp(i)), vals(pred), vals(pred) == y(rp(i)), yDisplay(rp(i)));
  
  fprintf('Program paused. Press enter to continue.\n');
  
  pause
  %end
  %pause
end

predictions = predictions';

reussite = score / m;


    % Number of true positives
    tp = sum(predictions == 1 & y == 1);

    % Number of false positives
    fp = sum(predictions == 1 & y == 0);
    
    % Number of false negatives
    fn = sum(predictions == 0 & y == 1);
    
    % Precision
    prec = tp / (tp + fp);
    
    % Recall
    rec = tp / (tp + fn);
    
    % F1 score
    F1 = (2 * prec * rec) / (prec + rec);
    
fprintf("Taux de reussite : %f\n", reussite * 100);
fprintf("Precision : %f\n", prec);
fprintf("Recall : %f\n", rec);
fprintf("F1 : %f\n", F1);

