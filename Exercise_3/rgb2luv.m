function luv= rgb2luv(rgb, wp)
%RGB2LUV Converts linear values to L*u*v* coordinates
%   Code from Stephen Westland "Computational Colour Science using MATLAB
%   2e"

xyz = rgb2xyz(lin2rgb(rgb));
white = whitepoint(wp);

luv = zeros(size(rgb));

% compute u' v' for sample
up = 4*xyz(:,1)./(xyz(:,1) + 15*xyz(:,2) + 3*xyz(:,3));
vp = 9*xyz(:,2)./(xyz(:,1) + 15*xyz(:,2) + 3*xyz(:,3));
% compute u' v' for white
upw = 4*white(1)/(white(1) + 15*white(2) + 3*white(3));
vpw = 9*white(2)/(white(1) + 15*white(2) + 3*white(3));
index = (xyz(:,2)/white(2) > 0.008856);
luv(:,1) = luv(:,1) + index.*(116*(xyz(:,2)/white(2)).^(1/3) - 16);  
luv(:,1) = luv(:,1) + (1-index).*(903.3*(xyz(:,2)/white(2)));  
luv(:,2) = 13*luv(:,1).*(up - upw);
luv(:,3) = 13*luv(:,1).*(vp - vpw);


end

