clear ; close all; clc

cd negatives;

files = ls;

j = 0;
for i=1:length(files)
  j += 1;
  fileName = strtrim(files(i, :));
  [IMG] = imread(fileName);
  IMG = cast(IMG, 'double');
  infos = imfinfo(fileName);
  if infos.TotalColors == 1
    IMG = IMG .* 255;
  end
  X(j, :) = IMG'(:);
  y(j, 1) = 0;
end;

cd ..;
cd screws;
files = ls;

for i=1:length(files)
  j += 1;
  fileName = strtrim(files(i, :));
  [IMG] = imread(fileName);
  IMG = cast(IMG, 'double');
  infos = imfinfo(fileName);
  if infos.TotalColors == 1
    IMG = IMG .* 255;
  end
  X(j, :) = IMG'(:);
  y(j, 1) = 1;
end;

cd ..;
cd leds;
files = ls;

for i=1:length(files)
  j += 1;
  fileName = strtrim(files(i, :));
  [IMG] = imread(fileName);
  IMG = cast(IMG, 'double');
  infos = imfinfo(fileName);
  if infos.TotalColors == 1
    IMG = IMG .* 255;
  end
  X(j, :) = IMG'(:);
  y(j, 1) = 2;
end;

cd ..;
cd numpad;
files = ls;

for i=1:length(files)
  j += 1;
  fileName = strtrim(files(i, :));
  [IMG] = imread(fileName);
  IMG = cast(IMG, 'double');
  infos = imfinfo(fileName);
  if infos.TotalColors == 1
    IMG = IMG .* 255;
  end
  X(j, :) = IMG'(:);
  y(j, 1) = 3;
end;

cd ..;

save -binary bundle X y
