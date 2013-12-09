clear ; close all; clc

cd "power_switch/UP";

files = ls;

j = 1;
for i=1:270
  filename = strcat(num2str(i), ".txt");
  if (exist(filename))
    eval(strcat("load\t", filename));
    eval(strcat("X(", num2str(j) ,", :) =\tX", num2str(i)), ";");
    j += 1;
  end;
end;

cd ../..
