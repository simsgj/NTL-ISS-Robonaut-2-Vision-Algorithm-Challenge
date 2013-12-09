clear ; close all; clc


load bundle



mu = mean(X, 2);
sigma = std(X, 1, 2);
sigma = sigma + (sigma == 0);
X = (X - repmat(mu, 1, size(X, 2))) ./ repmat(sigma, 1, size(X, 2));

[U, S] = pca(X);

%  Visualize the top 36 eigenvectors found
displayData(U(:, 1:49)');

fprintf('Program paused. Press enter to continue.\n');
pause;

K = 100;
Z = projectData(X, U, K);

X_rec  = recoverData(Z, U, K);

% Display normalized data
subplot(1, 2, 1);
displayOCT(X(1:49,:), 30);
title('Original');
axis square;

% Display reconstructed data from only k eigenfaces
subplot(1, 2, 2);
displayOCT(X_rec(1:49,:), 30);
title('Recovered');
axis square;
