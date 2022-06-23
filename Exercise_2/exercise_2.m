clear all
close all

load('../M.mat');


%% LMS <-> DKL matrices

% [0.950470000000000,1,1.08883000000000]
wp_d65 = whitepoint('d65');
lms_gray = wp_d65 * M_xyz2lms';

mc1 = lms_gray(1)/lms_gray(2);
mc2 = (lms_gray(1)+lms_gray(2))/lms_gray(3);

M_lms_dkl = [ 1  1 0;
               1 -mc1 0;
              -1 -1 mc2 ];

M_sheer = [ 1 0 0;
            0 1 0;
            0 0 1 ];

M_lms2dkl = (M_sheer * M_lms_dkl);
M_dkl2lms = inv(M_lms2dkl);

save('../M.mat','M_lms2dkl', 'M_dkl2lms','-append');
clear mc1 mc2 M_sheer M_lms_dkl

%% Isolating luminance mechanism

% For the D65 LMS grey point, DKL coordinates are: [Luminance, 0, 0], since
% it is an achromatic coordinate.

sdr_max_lum_weight = 0.5*max_lum_sdr./(sum(lms_gray(1:2))); % Weight to scale LMS grey coordinates to half intensity of SDR display
DKL_0 = (lms_gray.*sdr_max_lum_weight)*M_lms2dkl';
LMS_0 = DKL_0 * M_dkl2lms';
RGB_0 = LMS_0 * M_lms2rgb_sdr'; % RGB value corresponding to D65 grey point at half intensity of the SDR display


% If we want to increase the luminance without changing the values of
% chromatic mechanisms, we add our desired luminance change to the 1st
% coordinate of DKL_0

del_lum = 20;
DKL_1 = DKL_0 + [del_lum 0 0];
LMS_1 = DKL_1 * M_dkl2lms';
RGB_1 = LMS_1 * M_lms2rgb_sdr'; 
% RGB_0 and RGB_1 are coordinates of stimuli isolating the luminance
% mechansism

figure, imshowpair(linrgb2patch(RGB_0), linrgb2patch(RGB_1), 'montage');
title('Gray point (left); Grey point + Luminance (right)');


%% Isolating RG mechanism

% If we want to increase the response of red-green mechanim without
% changing the response of luminance of achromatic or yellow-violet
% mechanism, we add our desired change to the 2nd coordinate of DKL_0

del_rg = 20; % play around with this value, try negative values as well
DKL_2 = DKL_0 + [0 del_rg 0];
LMS_2 = DKL_2 * M_dkl2lms';
RGB_2 = LMS_2 * M_lms2rgb_sdr'; 
% RGB_0 and RGB_1 are coordinates of stimuli isolating the luminance
% mechansism

figure, imshowpair(linrgb2patch(RGB_0), linrgb2patch(RGB_2), 'montage');
title('Gray point (left); Grey point + RG response (right)');


%% Isolating YV mechanism

% If we want to increase the response of yellow-violet mechanim without
% changing the response of luminance of achromatic or red-green 
% mechanism, we add our desired change to the 3rd coordinate of DKL_0

del_yv = 70; % play around with this value, try negative values as well
DKL_3 = DKL_0 + [0 0 del_yv];
LMS_3 = DKL_3 * M_dkl2lms';
RGB_3 = LMS_3 * M_lms2rgb_sdr'; 
% RGB_0 and RGB_1 are coordinates of stimuli isolating the luminance
% mechansism

figure, imshowpair(linrgb2patch(RGB_0), linrgb2patch(RGB_3), 'montage');
title('Gray point (left); Grey point + YV response (right)');


%% Isolating RG mechaism with respect to saturated red point

RGB_red = [1 0 0];
LMS_red = RGB_red * M_rgb2lms_sdr';
DKL_red = ( LMS_red )*M_lms2dkl';

del_rg = -10; % Only negative values will work here. Increasin rg response further is not possible, because DKL_red is the most saturated red the display is able to produce
DKL_2 = DKL_red + [0 del_rg 0];
LMS_2 = DKL_2 * M_dkl2lms';
RGB_2 = LMS_2 * M_lms2rgb_sdr'; 
% RGB_0 and RGB_1 are coordinates of stimuli isolating the luminance
% mechansism

figure, imshowpair(linrgb2patch(RGB_red), linrgb2patch(RGB_2), 'montage');
title('Peak red (left); Peak red - RG response (right)');



