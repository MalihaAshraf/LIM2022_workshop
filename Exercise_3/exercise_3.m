clear all
close all

load('../M.mat');

%% Compare luminance/lightness values between different colour spaces

figure,
grey_px = 0:0.1:1;
rgb_grey = (cat(1, grey_px, grey_px, grey_px))';
colors = lin2rgb(rgb_grey);

% Luminance from LMS space
LMS_grey = rgb_grey * M_rgb2lms_sdr';
lum_grey = sum(LMS_grey(:,1:2), 2);
subplot(2, 2, 1)
scatter(grey_px, lum_grey, 20, colors, 'filled', 'MarkerEdgeColor', 'k'); hold on
plot([grey_px(1), grey_px(end)], [lum_grey(1), lum_grey(end)], '--k');
xlabel('Pixel value')
ylabel('Luminance (L+M) cd/m^2')
grid on
clear LMS_grey lum_grey

% Luminance from XYZ space
XYZ_grey = rgb_grey * M_rgb2xyz_sdr';
y_grey = XYZ_grey(:,2); % Y is the luminance in XYZ space
subplot(2, 2, 2)
scatter(grey_px, y_grey, 20, colors, 'filled', 'MarkerEdgeColor', 'k'); hold on
plot([grey_px(1), grey_px(end)], [y_grey(1), y_grey(end)], '--k');
xlabel('Pixel value')
ylabel('Luminance (Y) cd/m^2')
grid on
clear XYZ_grey y_grey

% Lightness from L*a*b* space
Lab_grey = rgb2lab(lin2rgb(rgb_grey), 'WhitePoint', 'd65');
L_grey = Lab_grey(:,1);
subplot(2, 2, 3)
scatter(grey_px, L_grey, 20, colors, 'filled', 'MarkerEdgeColor', 'k'); hold on
plot([grey_px(1), grey_px(end)], [L_grey(1), L_grey(end)], '--k');
xlabel('Pixel value')
ylabel('Lightenss (Lab space)')
grid on
clear Lab_grey L_grey


% Lightness from L*u*v* space
Luv_grey = rgb2luv(rgb_grey, 'd65');
L2_grey = Luv_grey(:,1);
subplot(2, 2, 4)
scatter(grey_px, L2_grey, 20, colors, 'filled', 'MarkerEdgeColor', 'k'); hold on
plot([grey_px(1), grey_px(end)], [L2_grey(1), L2_grey(end)], '--k');
xlabel('Pixel value')
ylabel('Lightenss (Luv space)')
grid on
clear Luv_grey L2_grey

% Lightness values are not a linear transform of linear RGB values


%% Compare chromaticity coordinaes between different colorspaces

% Generate uniformly spaced chromaticity coordinates in a*b* plane with 
% L* = 76

L_mid = 76;  % play around with this value to view different gamuts at different lightness levels


[a_mat, b_mat] = meshgrid(-100:10:100, -100:10:100);
ab = [reshape(a_mat, numel(a_mat), 1), reshape(b_mat, numel(b_mat), 1)];
r = sqrt( ab(:,1).^2 + ab(:,2).^2);
ab_chrom = ab(r<= 100, :);
Lab_chrom = [L_mid*ones(length(ab_chrom),1), ab_chrom];
colors = lab2rgb(Lab_chrom, 'WhitePoint', 'd65');
RGB_chrom = rgb2lin(colors); % linear RGB
clear a_mat b_mat r ab

remove_out_of_gamut = true;
if remove_out_of_gamut
    % removes colors which are beyond the display's gamut
    criteria = ~logical(sum(colors < -0, 2));
    Lab_chrom = Lab_chrom(criteria, :);
    ab_chrom = ab_chrom(criteria, :);
    colors = colors(criteria, :);
    RGB_chrom = RGB_chrom(criteria, :);
end

figure, 
% Plot chromaticity coordinates in Lab space
subplot(2, 2, 1)
scatter(ab_chrom(:,1), ab_chrom(:,2), 20, colors, 'filled', 'MarkerEdgeColor', 'k'); hold on
axis square
xlabel('-a* ---------- 0 ---------- +a*')
ylabel('-b* ---------- 0 ---------- +b*')



% Plot chromaticity coordinates in Luv space
subplot(2, 2, 2)
Luv_chrom = rgb2luv(RGB_chrom, 'd65');
uv_chrom = Luv_chrom(:, 2:3); 
scatter(uv_chrom(:,1), uv_chrom(:,2), 20, colors, 'filled', 'MarkerEdgeColor', 'k'); hold on
axis square
xlabel('-u* ---------- 0 ---------- +u*')
ylabel('-v* ---------- 0 ---------- +v*')



% Plot chromaticity coordinates in xy space
subplot(2, 2, 3)
XYZ_chrom = lab2xyz(Lab_chrom, 'WhitePoint', 'd65');
xy_chrom = XYZ2Yxy(XYZ_chrom);
xy_chrom = xy_chrom(:,2:3);
scatter(xy_chrom(:,1), xy_chrom(:,2), 20, colors, 'filled', 'MarkerEdgeColor', 'k'); hold on
axis square
xlabel('x')
ylabel('y')


% Plot chromaticity coordinates in opponent color mechanism
subplot(2, 2, 4)
LMS_chrom = RGB_chrom * M_rgb2lms_sdr';
DKL_chrom = LMS_chrom * M_lms2dkl';
opp_chrom = DKL_chrom(:,2:3);
scatter(opp_chrom(:,1), opp_chrom(:,2), 20, colors, 'filled', 'MarkerEdgeColor', 'k'); hold on
axis square
xlabel('L-M')
ylabel('S-(L+M)')



