function [res] = calculateError(data_input, theta, data_output)
  res = 0;
  prediction = data_input * theta;
  for i = 1:size(data_input, 1)
    res += sqrt((prediction(i, 1) - data_output(i, 1)) .* (prediction(i, 1) - data_output(i, 1)) + (prediction(i, 2) - data_output(i, 2)) .* (prediction(i, 2) - data_output(i, 2)));
  end
  res /= (size(data_input, 1) * 1.0);
end
