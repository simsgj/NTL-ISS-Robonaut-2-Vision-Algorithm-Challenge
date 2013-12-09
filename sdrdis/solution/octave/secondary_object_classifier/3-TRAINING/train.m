function train(hidden_layer_size, lambda, nb_iterations, fraction)

% Load Training Data
fprintf('Loading and Visualizing Data ...\n')

load('bundlePCA');

[vals, ids, idx] = unique(y);
num_labels = size(vals, 1);
y = idx;


X = X(1:round(size(X, 1) * fraction), :);
y = y(1:round(length(y) * fraction));
m = size(X, 1);
input_layer_size = size(X, 2);

% Randomly select 100 data points to display
%if (false) {
%  fprintf('Loading and Visualizing Data ...\n');
%  sel = randperm(size(X, 1));
%  sel = sel(1:100);
%
%  displayData(X(sel, :));
%}

fprintf('\nInitializing Neural Network Parameters ...\n')

Theta1 = randInitializeWeights(input_layer_size, hidden_layer_size);
Theta2 = randInitializeWeights(hidden_layer_size, num_labels);

% Unroll parameters
initial_nn_params = [Theta1(:) ; Theta2(:)];

fprintf('\nTraining Neural Network... \n')

nb = 0;
step = 25;
while (nb * step <= nb_iterations)

  nb+= 1;
  
  fileName = ["trained/", num2str(hidden_layer_size), "-", num2str(lambda), '-', num2str(nb), ".mat"];

  if (exist(fileName) == 2)
    str = sprintf("load -binary %s",fileName);
    eval(str);
    %fprintf ("STEP: %d completed, cost: %f...\n", nb, cost(end));
    initial_nn_params = [Theta1(:) ; Theta2(:)];
  else
  

    options = optimset('MaxIter', step);

    costFunction = @(p) nnCostFunction(p, ...
                                       input_layer_size, ...
                                       hidden_layer_size, ...
                                       num_labels, X, y, lambda);

    [nn_params, cost] = fmincg(costFunction, initial_nn_params, options);
    
    initial_nn_params = nn_params;

    Theta1 = reshape(nn_params(1:hidden_layer_size * (input_layer_size + 1)), ...
                     hidden_layer_size, (input_layer_size + 1));

    Theta2 = reshape(nn_params((1 + (hidden_layer_size * (input_layer_size + 1))):end), ...
                     num_labels, (hidden_layer_size + 1));
                     
    
    str = sprintf("save -binary %s hidden_layer_size lambda Theta1 Theta2 cost",fileName);
    eval(str);
    
    fprintf ("STEP: %d completed, cost: %f...\n", nb, cost(end));
  end;
end
