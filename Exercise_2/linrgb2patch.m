function [img] = linrgb2patch(lin_rgb)
%LINRGB2PATCH transforms linear RGB to sRGB and outputs a 100x100 color
%patch

rgb = (lin2rgb(lin_rgb));
img = double(cat(3, ones(100,100).*rgb(1), ones(100,100).*rgb(2), ones(100,100).*rgb(3)));

end

