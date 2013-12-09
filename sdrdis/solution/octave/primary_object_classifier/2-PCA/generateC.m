function generateC(className, U, nb)
  # H file

  file_id = fopen(strcat('c/', tolower(className), '.h'), 'w');
  
  defineClassName = strcat(toupper(className), '_H');
  
  fprintf(file_id, "#ifndef %s\n", defineClassName);
  fprintf(file_id, "#define %s\n\n", defineClassName);
  
  fprintf(file_id, "#include \"pca.h\"\n\n");
  
  fprintf(file_id, "class %s : public PCA\n{\npublic:\n    %s();\n};\n\n", className, className);
  
  fprintf(file_id, "#endif // %s", defineClassName);
  
  fclose(file_id);
  
  # CPP file
  
  file_id = fopen(strcat('c/', tolower(className), '.cpp'), 'w');
  
  fprintf(file_id, "#include \"%s\"\n\n", strcat(tolower(className), '.h'));
  
  eigenVectorsConstantName = strcat(toupper(className), "_EIGENVECTORS");
  
  fprintf(file_id, "static float %s[%d] = {\n", eigenVectorsConstantName, nb * size(U, 2));
  
  str = "";
  for i = 1:nb
    #fprintf(file_id, "vector < float > eigen%d;\n", i);
    for j = 1:size(U, 2)
      if (j != 1 || i != 1)
        str = [str, ", "];
      end;
      str = [str, num2str(U(j, i))];
    end
    str = [str, "\n"];
    #fprintf(file_id, "this->eigenvectors.push_back(eigen%d);\n", i);
  end
  
  fprintf(file_id, "%s};\n\n", str);
  
  fprintf(file_id, "%s::%s() {\n    eigenVectors = new float[%d];\n    eigenVectors = %s;\n    nbEigenVectors = %d;\n    nbInputs = %d;\n}\n\n", className, className,  nb * size(U, 2), eigenVectorsConstantName, nb, size(U, 2));
  
  fclose(file_id);
end;
