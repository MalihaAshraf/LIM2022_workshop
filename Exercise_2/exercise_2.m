clear all
close all

load('../M.mat');
M_rgb2lms = M_rgb2lms_sdr;
M_lms2rgb = M_lms2rgb_sdr;

%% LMS <-> DKL matrices

% Based on Brainard's tutorial (see pdf in refs)
% chose background - necessary for DKL space since the transformation
% matrix is defined on incremental cone signals

rgb_grey = [0.5 0.5 0.5]; 
lms_grey = M_rgb2lms * rgb_grey';
bg = lms_grey;

mc1 = bg(1)/bg(2);
mc2 = (bg(1)+bg(2))/bg(3);

M_lms2dkl_raw = [ 1  1 0;
               1 -mc1 0;
              -1 -1 mc2 ];

M_dkl2lms_raw = inv(M_lms2dkl_raw);

% Scale matrix to obtain unit outputs for an input vector of unit length
% with respect to the specified grey background
achr_unit = M_dkl2lms_raw(:,1) / norm(M_dkl2lms_raw(:,1)./bg);
rg_unit = M_dkl2lms_raw(:,2) / norm(M_dkl2lms_raw(:,2)./bg);
yv_unit = M_dkl2lms_raw(:,3) / norm(M_dkl2lms_raw(:,3)./bg);

lum_resp_raw = M_lms2dkl_raw * achr_unit % not normalized
l_minus_m_resp_raw = M_lms2dkl_raw * rg_unit % not normalized
s_minus_lum_resp_raw = M_lms2dkl_raw * yv_unit % not normalized

% Rescale the raw matrix
M_rescale = [1/lum_resp_raw(1) 0 0; ...
             0 1/l_minus_m_resp_raw(2) 0; ...
             0 0 1/s_minus_lum_resp_raw(3)];

M_lms2dkl = M_rescale * M_lms2dkl_raw;

% Test scaling
% All the following calculations should produce unit responses
lum_resp = M_lms2dkl*achr_unit
l_minus_m_resp = M_lms2dkl*rg_unit
s_minus_lum_resp = M_lms2dkl * yv_unit

M_dkl2lms = inv(M_lms2dkl);


% Try different stimuli
M_dkl2lms * [0 0 0]'
M_dkl2lms * [1 0 0]'   %isolates L+M mechanism
M_dkl2lms * [0 1 0]'   %isolates L-M mechanism
M_dkl2lms * [0 0 1]'   %isolates S-(L+M) mechanism


M_lms2dkl * lms_grey   %a purely achromatic modulation results in activation of the L+M channel only

% Check whether this modulation is within gamut
M_lms2rgb*M_dkl2lms * [0 0 0.5]'  % Values outside  of [0-1] range means the color is out of gamut  


save('../M.mat','M_lms2dkl', 'M_dkl2lms','-append');
clear mc1 mc2 

%% Isolating luminance mechanism

% For the LMS grey point, DKL coordinates are: [Luminance, 0, 0], since
% it is an achromatic coordinate.

LMS_0 = lms_grey';
RGB_0 = rgb_grey; % RGB value corresponding to grey point at half intensity of the SDR display


% If we want to increase the luminance without changing the values of
% chromatic mechanisms, we add our desired luminance change to the 1st
% coordinate of the DKL value

del_lum = 1; % normalized dkl units. Play around with this value to change the lightness of the color patch. Test negtive values
DKL_1 =  [del_lum 0 0];
LMS_1 = DKL_1 * M_dkl2lms';
LMS_stim = LMS_0 + LMS_1;
RGB_1 = LMS_stim * M_lms2rgb_sdr'; 
% RGB_0 and RGB_1 are coordinates of stimuli isolating the luminance
% mechansism

figure, imshowpair(linrgb2patch(RGB_0), linrgb2patch(RGB_1), 'montage');
title('Grey point (left); Grey point + Luminance (right)');


%% Isolating RG mechanism

% If we want to increase the response of red-green mechanim without
% changing the response of luminance of achromatic or yellow-violet
% mechanism, we add our desired change to the 2nd coordinate of the DKL
% value

del_rg = 0.1; % play around with this value, try negative values as well
DKL_2 = [0 del_rg 0];
LMS_2 = DKL_2 * M_dkl2lms';
LMS_stim = LMS_0 + LMS_2;
RGB_2 = LMS_stim * M_lms2rgb_sdr'; 
% RGB_0 and RGB_1 are coordinates of stimuli isolating the luminance
% mechansism

figure, imshowpair(linrgb2patch(RGB_0), linrgb2patch(RGB_2), 'montage');
title('grey point (left); Grey point + RG response (right)');


%% Isolating YV mechanism

% If we want to increase the response of yellow-violet mechanim without
% changing the response of luminance of achromatic or red-green 
% mechanism, we add our desired change to the 3rd coordinate of the DKL
% value

del_yv = -0.5; % play around with this value, try negative values as well
DKL_3 = [0 0 del_yv];
LMS_3 = DKL_3 * M_dkl2lms';
LMS_stim = LMS_0 + LMS_3;
RGB_3 = LMS_stim * M_lms2rgb_sdr'; 
% RGB_0 and RGB_1 are coordinates of stimuli isolating the luminance
% mechansism

figure, imshowpair(linrgb2patch(RGB_0), linrgb2patch(RGB_3), 'montage');
title('grey point (left); Grey point + YV response (right)');

