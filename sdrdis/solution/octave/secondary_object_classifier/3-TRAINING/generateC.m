function generateC(className, nb_hidden_layer, lambda, nb_trained)
  fileName = ["trained/", num2str(nb_hidden_layer), "-", num2str(lambda), "-", num2str(nb_trained), ".mat"];
  str = sprintf("load -binary %s",fileName);
  eval(str);
  
  Theta1 = Theta1';
  Theta2 = Theta2';
  
  # H file

  file_id = fopen(strcat('c/', tolower(className), '.h'), 'w');
  
  defineClassName = strcat(toupper(className), '_H');
  
  fprintf(file_id, "#ifndef %s\n", defineClassName);
  fprintf(file_id, "#define %s\n\n", defineClassName);
  
  fprintf(file_id, "#include \"neuralnetwork.h\"\n\n");
  
  fprintf(file_id, "class %s : public NeuralNetwork\n{\npublic:\n    %s();\n};\n\n", className, className);
  
  fprintf(file_id, "#endif // %s", defineClassName);
  
  fclose(file_id);
  
  # CPP file
  file_id = fopen(strcat('c/', tolower(className), '.cpp'), 'w');
  
  fprintf(file_id, "#include \"%s\"\n\n", strcat(tolower(className), '.h'));
  
  layersConstantName = strcat(toupper(className), "_LAYERS");

  fprintf(file_id, "static size_t %s[3] = {%d, %d, %d};\n", layersConstantName, size(Theta1, 1) - 1, nb_hidden_layer, size(Theta2, 2));
  
  weightsConstantNamePrefix = strcat(toupper(className), "_WEIGHTS_");
  
  str = "";
  for i = 1:size(Theta1, 1)
    for j = 1:size(Theta1, 2)
      if (j != 1 || i != 1)
        str = [str, ", "];
      end;
      str = [str, num2str(Theta1(i, j))];
    end;
  end;
  
  fprintf(file_id, "static float %s1[%d] = {%s};\n", weightsConstantNamePrefix, size(Theta1, 1) * size(Theta1, 2), str);
  
  str = "";
  for i = 1:size(Theta2, 1)
    for j = 1:size(Theta2, 2)
      if (j != 1 || i != 1)
        str = [str, ", "];
      end;
      str = [str, num2str(Theta2(i, j))];
    end;
  end;
  
  fprintf(file_id, "static float %s2[%d] = {%s};\n\n", weightsConstantNamePrefix, size(Theta2, 1) * size(Theta2, 2), str);

  
  
  
  
  
  
  load "bundlePCA"
  [vals, ids, idx] = unique(y);
  X = X((round(size(X, 1) * 0.8) + 1):end, :);
  y = y((round(length(y) * 0.8) + 1):end);
  m = size(X, 1);

  Theta1 = Theta1';
  Theta2 = Theta2';

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
  

  fprintf(file_id, "%s::%s() {\n", className, className);
  fprintf(file_id, "    layers = new size_t[3];\n");
  fprintf(file_id, "    layers = %s;\n", layersConstantName);
  fprintf(file_id, "    nbLayers = 3;\n\n");
  
  fprintf(file_id, "    weights = new float * [2];\n");
  fprintf(file_id, "    weights[0] = new float [%d];\n", size(Theta1, 1) * size(Theta1, 2));
  fprintf(file_id, "    weights[0] = %s1;\n", weightsConstantNamePrefix);
  fprintf(file_id, "    weights[1] = new float [%d];\n", size(Theta2, 1) * size(Theta2, 2));
  fprintf(file_id, "    weights[1] = %s2;\n\n", weightsConstantNamePrefix);
  
  fprintf(file_id, "    precision = %s;\n", num2str(prec));
  fprintf(file_id, "    recall = %s;\n", num2str(rec));
  fprintf(file_id, "}\n");
  
  fclose(file_id);
end;
