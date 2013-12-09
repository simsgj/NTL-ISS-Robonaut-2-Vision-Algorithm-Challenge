function training(from_index)

for nb_hidden_layer = (from_index * 40 + 70):20:((from_index + 1) * 40 + 69)
  for lambda = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30]
    train(nb_hidden_layer, lambda, 75, 0.7);
  end
end;

end
