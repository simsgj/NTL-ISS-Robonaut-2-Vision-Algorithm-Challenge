function displayOCT(X, imgSize)
  X = (X + 2) / 4;
  X = min(X, 1);
  X = max(X, 0);
  #sample = reshape(X, 30, 30);
  displayData(X, imgSize);
  #imshow(sample);
end
