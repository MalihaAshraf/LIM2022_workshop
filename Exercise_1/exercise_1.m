% clear all;
% close all;

%% Display spectra

% Load display spectra
% File structure: Column 1: Wavelength (lambda), Column 2: Spectrum of full
% white, Columns 3-5: Spectra of RGB channels respectively
spec1 = csvread( 'spectrum_sdr_display.csv' ); % Load SDR display
spec2 = csvread( 'spectrum_hdr_display.csv' ); % Load HDR display
lambda_display = spec1(:,1);

% Plot display spectra
figure, plot( lambda_display, spec1(:,2), 'DisplayName', 'SDR display');
hold on,
plot( lambda_display, spec2(:,2), 'DisplayName', 'HDR display');
xlabel('Wavelength (\lambda)');
ylabel('W/sr.m^{-2}');
legend boxoff
title('Spectral power distribution - full white')
axis([300 1000 0 0.024]);


figure,
% Plot RGB spectra for SDR display
plot( lambda_display, spec1(:,3), '-r', 'DisplayName', 'R'); hold on
plot( lambda_display, spec1(:,4), '-g', 'DisplayName', 'G'); hold on
plot( lambda_display, spec1(:,5), '-b', 'DisplayName', 'B'); hold on
xlabel('Wavelength (\lambda)');
ylabel('W/sr.m^{-2}');
legend boxoff
title('Spectral power distribution - RGB channels')
axis([300 1000 0 0.01]);

% Plot sum of spectra
plot( lambda_display, spec1(:,3)+spec1(:,4)+spec1(:,5), '--k', 'DisplayName', 'R+G+B'); hold on

% Plot white
plot( lambda_display, spec1(:,2), '-m','DisplayName', 'White'); hold on,


%% Load Cone fundamentals and luminous efficiency curve

% Load CIE LMS 2006 cone fundamentals. Source: cvrl.org 
cmf_lms = csvread( 'linss2_10e_1.csv' );
% Interpolate LMS response to match wavelength range of the display spectra
cmf_lms = interp1( cmf_lms(:, 1), cmf_lms(:, 2:4), lambda_display, 'linear', 0);
% Load and interpolate luminance response curve
v_lambda = csvread( 'vl_linCIE2008v2e_1.csv' );
v_lambda = interp1( v_lambda(:,1), v_lambda(:,2), lambda_display, 'linear', 0);


%% Calculating cone signals

figure,
% Plot display spectra
subplot(2,2,1), hold on
plot( lambda_display, spec1(:,2), 'DisplayName', 'SDR display');
plot( lambda_display, spec2(:,2), 'DisplayName', 'HDR display');
xlabel('Wavelength (\lambda)');
ylabel('W/sr.m^{-2}');
legend boxoff
title('Spectral power distribution - full white')
axis([300 800 0 0.024]);

% Plot cone abosrption curves
subplot(2,2,2)
plot( lambda_display, cmf_lms(:,1), '-r', 'DisplayName', 'L'); hold on
plot( lambda_display, cmf_lms(:,2), '-g', 'DisplayName', 'M');
plot( lambda_display, cmf_lms(:,3), '-b', 'DisplayName', 'S');
% plot( lambda_display, v_lambda, '-k', 'DisplayName', 'Luminance');
[val ,idx] = max(cmf_lms); 
plot(lambda_display(idx(1)), 0,'r^','DisplayName', num2str(lambda_display(idx(1))));
plot(lambda_display(idx(2)), 0,'g^','DisplayName', num2str(lambda_display(idx(2))));
plot(lambda_display(idx(3)), 0,'b^','DisplayName', num2str(lambda_display(idx(3))));
xlim( [ 300 800 ])
xlabel('Wavelength (\lambda)');
ylabel('Normalized cone respones');
legend boxoff
title('LMS cone sensitivities')

% Calculate colorimetric values
% Non-normalized LMS values for full white
subplot(2,2,3)
k = 683.002; %1 watt at 555nm (= peak of V*(lambda)) = 683 lumen
LMS_sdr = trapz(lambda_display, cmf_lms.*spec1(:,2)).*k; 
LMS_hdr = trapz(lambda_display, cmf_lms.*spec2(:,2)).*k;
% Plot non-normalised LMS for both HDR and SDR display
bar(1,LMS_sdr(1),'r','DisplayName', 'L'); hold on; 
bar(2,LMS_sdr(2),'g','DisplayName', 'M'); 
bar(3,LMS_sdr(3),'b','DisplayName', 'S'); 
xticks([1 2 3]);
xticklabels({'L','M','S'});
axis([0 4 0 800]);
ylabel('Arbitary units');
legend off
title('SDR')


subplot(2,2,4)
bar(1,LMS_hdr(1),'r','DisplayName', 'L'); hold on; 
bar(2,LMS_hdr(2),'g','DisplayName', 'M'); 
bar(3,LMS_hdr(3),'b','DisplayName', 'S'); 
xticks([1 2 3]);
xticklabels({'L','M','S'});
axis([0 4 0 800]);
ylabel('Arbitary units');
legend off
title('HDR')


% Plot luminous efficiency curves
figure,
plot(lambda_display, v_lambda, '-k', 'DisplayName', 'V(\lambda)'); hold on;
[val,idx]=max(v_lambda); 
plot(lambda_display(idx), 0,'k^','DisplayName', num2str(lambda_display(idx))); 
xlim( [ 300 800 ])
xlabel('Wavelength (\lambda)');
ylabel('Luminous efficiency');
legend boxoff
title('(Photopic) Luminous Efficiency');



%% Calculating cone signals - normalized

figure,

% Plot display spectra
subplot(2,2,1), hold on
plot( lambda_display, spec1(:,2), 'DisplayName', 'SDR display');
plot( lambda_display, spec2(:,2), 'DisplayName', 'HDR display');
xlabel('Wavelength (\lambda)');
ylabel('W/sr.m^{-2}');
legend boxoff
title('Spectral power distribution - full white')
axis([300 800 0 0.024]);

% Normalised L M S 
subplot(2,2,2) %normalised L M S 
lms_weights = [0.689903 0.348322 0.0371597];
plot(lambda_display, v_lambda, '-k', 'DisplayName', 'V(\lambda)'); hold on;
plot( lambda_display, lms_weights(1)*cmf_lms(:,1), '-r', 'DisplayName', 'L'); 
plot( lambda_display, lms_weights(2)*cmf_lms(:,2), '-g', 'DisplayName', 'M');
plot( lambda_display, lms_weights(3)*cmf_lms(:,3), '-b', 'DisplayName', 'S');
xlabel('Wavelength (\lambda)');
ylabel('Relative efficiency');
xlim( [ 300 800 ])
legend boxoff
title('LMS cone sensitivities');


% Calculate colorimetric values
% Normalized LMS values for full white
subplot(2,2,3)
k=683; %1 watt at 555nm (= peak of V*(lambda)) = 683 lumen
LMS_sdr = trapz(lambda_display, cmf_lms.*spec1(:,2)).*k; 
LMS_hdr = trapz(lambda_display, cmf_lms.*spec2(:,2)).*k;
%LMS normalization weights
lms_weights = [0.689903 0.348322 0.0371597];
LMS_sdr_norm = lms_weights.*LMS_sdr;   %[159.520593696154,69.9403299121463,4.80661995045993]
% Plot normalised LMS for both HDR and SDR display
bar(1,LMS_sdr_norm(1),'r'); hold on; 
bar(2,LMS_sdr_norm(2),'g'); 
bar(3,LMS_sdr_norm(3),'b'); 
text(1,LMS_sdr_norm(1),num2str(LMS_sdr_norm(1),'   %.1f'),'VerticalAlignment','bottom', 'HorizontalAlignment','center');
text(2,LMS_sdr_norm(2),num2str(LMS_sdr_norm(2),'   %.1f'),'VerticalAlignment','bottom', 'HorizontalAlignment','center');
text(3,LMS_sdr_norm(3),num2str(LMS_sdr_norm(3),'   %.1f'),'VerticalAlignment','bottom', 'HorizontalAlignment','center');
xticks([1 2 3]);
xticklabels({'L','M','S'});
axis([0 4 0 500]);
ylabel('cd/m^2');
legend off
title('SDR')

subplot(2,2,4)
LMS_hdr_norm = lms_weights.*LMS_hdr;
bar(1,LMS_hdr_norm(1),'r','DisplayName', 'L'); hold on; 
bar(2,LMS_hdr_norm(2),'g','DisplayName', 'M'); 
bar(3,LMS_hdr_norm(3),'b','DisplayName', 'S'); 
xticks([1 2 3]);
xticklabels({'L','M','S'});
axis([0 4 0 500]);
ylabel('cd/m^2');
legend off
title('HDR')


%% Luminance derived from L and M signals

% LMS normalization weights
lms_weights = [0.689903 0.348322 0.0371597];
LMS_sdr_norm = lms_weights.*LMS_sdr;   %159.5206   69.9403    4.8066
LMS_hdr_norm = lms_weights.*LMS_hdr;

% Luminance values from normalized LMS values (L+M)
lum_sdr_2 = sum(LMS_sdr_norm(1:2)); % value should be equal to max_lum_sdr   max_lum=229.4609
lum_hdr_2 = sum(LMS_hdr_norm(1:2)); % value should be equal to max_lum_hdr

%% Luminance values from the V(lambda) response function for full white, R,G,B

max_lum_sdr(1) = trapz(lambda_display, v_lambda.*spec1(:,2)).*683.002;  %229.4609
max_lum_sdr(2) = trapz(lambda_display, v_lambda.*spec1(:,3)).*683.002   %
max_lum_sdr(3) = trapz(lambda_display, v_lambda.*spec1(:,4)).*683.002   %  
max_lum_sdr(4) = trapz(lambda_display, v_lambda.*spec1(:,5)).*683.002    %
max_lum_sdr   % W R G B    %229.4615   47.0419  156.1564   21.5444



%% Calculating luminance values for full white for HDR and SDR 

figure,

% Plot display spectra
subplot(2,2,1)
plot( lambda_display, spec1(:,2), 'DisplayName', 'SDR display');
hold on,
plot( lambda_display, spec2(:,2), 'DisplayName', 'HDR display');
xlabel('Wavelength (\lambda)');
ylabel('W/sr.m^{-2}');
legend boxoff
title('Spectral power distribution - full white')
axis([300 800 0 0.024]);

% Plot V (lambda)
subplot(2,2,2) 
plot(lambda_display, v_lambda, '-k', 'DisplayName', 'V(\lambda)'); hold on;
xlabel('Wavelength (\lambda)');
ylabel('Relative efficiency');
xlim( [ 300 800 ])
legend boxoff
title('Luminous Efficiency');


% Plot calculated luminance from integration of white spectra
subplot(2,2,3)
b=bar([max_lum_sdr(1)]);
b.FaceColor = 'flat';
text(1,300,num2str(max_lum_sdr(1)),'HorizontalAlignment','center'); 
xticks([1]);
xticklabels({'SDR'});
axis([0 2 0 500]);
ylabel('cd/m^2');
legend off
title('Luminance')

% Plot calculated luminance from normalized L+M values
subplot(2,2,4)
b=bar([lum_sdr_2]);
b.FaceColor = 'flat';
text(1,400, num2str(LMS_sdr_norm(1)),'HorizontalAlignment','center');
text(1,300, num2str(LMS_sdr_norm(2)),'HorizontalAlignment','center');
xticks([1]);
xticklabels({'SDR'});
axis([0 2 0 500]);
ylabel('cd/m^2');
legend off
title('L+M')


%% Get RGB2LMS matrix from spectra (SDR only)

RGB_sdr_response = spec1(:,3:5);

% Calculate RGB -> LMS matrix from display RGB spectra and LMS responses
cmf_lms_sc = cmf_lms.*lms_weights;
M_rgb2lms_sdr = ((RGB_sdr_response'*cmf_lms_sc))'.*683.002;

% Calculate LMS values from RGB
lms_test1 = [1 1 1]*M_rgb2lms_sdr'  % should be equal to LMS_sdr_norm
lms_R_sdr = [1 0 0]*M_rgb2lms_sdr'
lms_G_sdr = [0 1 0]*M_rgb2lms_sdr'
lms_B_sdr = [0 0 1]*M_rgb2lms_sdr'

% Calculate LMS2RGB matrix
M_lms2rgb_sdr = pinv(M_rgb2lms_sdr);
M_lms2rgb = M_lms2rgb_sdr;
M_rgb2lms = M_rgb2lms_sdr;


%% Cone isolating values

% Test RGB values for specified LMS coordinates (SDR display)
LMS_1 = [100, 50, 2];
RGB_test1 = LMS_1 * M_lms2rgb_sdr'; 

% If we want to change only the value of M while keeping L and S responses
% the same
LMS_2 = [120, 50, 2];
RGB_test2 = LMS_2 * M_lms2rgb_sdr'; 

% RGB_test1 and RGB_test2 are L cone isolating coordinates


%% Calculate XYZ <-> RGB matrices and XYZ <-> LMS using the same method above

cmf_xyz = csvread( 'ciexyz31.csv' );
% Interpolate XYZ response to match wavelength range of the display spectra
cmf_xyz = interp1( cmf_xyz(:, 1), cmf_xyz(:, 2:4), lambda_display, 'linear', 0);


% Calculate RGB <-> XYZ matrix from display RGB spectra and LMS responses
M_rgb2xyz_sdr = ((RGB_sdr_response'*cmf_xyz))'.*683.002;
M_xyz2rgb_sdr = pinv(M_rgb2xyz_sdr);


% Calculate LMS <-> XYZ matrix from XYZ and LMS responses
M_lms2xyz = ((cmf_lms_sc'*cmf_xyz))'.*683.002;
M_xyz2lms = pinv(M_lms2xyz);



%% Save all color transformation matrices

save('../M.mat', 'M_lms2rgb_sdr', 'M_lms2xyz',...
    'M_rgb2lms_sdr', 'M_rgb2xyz_sdr', 'M_xyz2lms',... 
    'M_xyz2rgb_sdr');



%% CALCULATING CONE SIGNALS and Luminance for R,G,B using normalised LMS and SPD - SDR only

figure,  

% Plot display spectra
subplot(2,2,1) 
plot( lambda_display, spec1(:,3), ':','DisplayName', 'R'); hold on;
plot( lambda_display, spec1(:,4),  ':','DisplayName', 'G');
plot( lambda_display, spec1(:,5), ':', 'DisplayName', 'B');
plot( lambda_display, spec1(:,2), ':', 'DisplayName', 'W');
xlabel('Wavelength (\lambda)');
ylabel('W/sr.m^{-2}');
legend boxoff
title('Spectral power distribution - R G B')
axis([300 800 0 0.01]);

% Plot normalized LMS Color matching function
subplot(2,2,2)
lms_weights = [0.689903 0.348322 0.0371597];
plot(lambda_display, v_lambda, '-k', 'DisplayName', 'V(\lambda)'); hold on;
plot( lambda_display, lms_weights(1)*cmf_lms(:,1), '-r', 'DisplayName', 'L'); 
plot( lambda_display, lms_weights(2)*cmf_lms(:,2), '-g', 'DisplayName', 'M');
plot( lambda_display, lms_weights(3)*cmf_lms(:,3), '-b', 'DisplayName', 'S');
xlabel('Wavelength (\lambda)');
ylabel('Relative efficiency');
xlim( [ 300 800 ])
legend boxoff
title('LMS cone sensitivities');

% Plot calculated LMS outputs for RGBW spectra using SPD integration
% melthod
subplot(2,2,3) 
% Calculate colorimetric values
% Normalized LMS values for R G B 
k=683; %1 watt at 555nm (= peak of V*(lambda)) = 683 lumen
LMS_sdr_R = trapz(lambda_display, cmf_lms.*spec1(:,3)).*k; 
LMS_sdr_G = trapz(lambda_display, cmf_lms.*spec1(:,4)).*k; 
LMS_sdr_B = trapz(lambda_display, cmf_lms.*spec1(:,5)).*k; 
LMS_sdr_W = trapz(lambda_display, cmf_lms.*spec1(:,2)).*k; 
%LMS normalization weights
lms_weights = [0.689903 0.348322 0.0371597];
LMS_sdr_R_norm = lms_weights.*LMS_sdr_R;    % 40.5264    6.5153    0.0900
LMS_sdr_G_norm = lms_weights.*LMS_sdr_G;    % 103.2220   52.9341    0.4375
LMS_sdr_B_norm = lms_weights.*LMS_sdr_B;    % 12.3159    9.2284    4.1035
LMS_sdr_W_norm = lms_weights.*LMS_sdr_W;    % 
% Plot LMS for R G B of SDR display
b=bar([1:15], [LMS_sdr_R_norm 0 LMS_sdr_G_norm 0 LMS_sdr_B_norm 0 LMS_sdr_W_norm]); hold on;
text(1,LMS_sdr_R_norm(1),num2str(LMS_sdr_R_norm(1),'%.1f'),'Rotation', 90);
text(2,LMS_sdr_R_norm(2),num2str(LMS_sdr_R_norm(2),'%.1f'),'Rotation', 90);
text(3,LMS_sdr_R_norm(3),num2str(LMS_sdr_R_norm(3),'%.1f'),'Rotation', 90);
b.FaceColor = 'flat';
b.CData(1,:) = [1 0 0];
b.CData(2,:) = [0 1 0];
b.CData(3,:) = [0 0 1];
text(5,LMS_sdr_G_norm(1),num2str(LMS_sdr_G_norm(1),'%.1f'),'Rotation', 90);
text(6,LMS_sdr_G_norm(2),num2str(LMS_sdr_G_norm(2),'%.1f'),'Rotation', 90);
text(7,LMS_sdr_G_norm(3),num2str(LMS_sdr_G_norm(3),'%.1f'),'Rotation', 90);
b.FaceColor = 'flat';
b.CData(5,:) = [1 0 0];
b.CData(6,:) = [0 1 0];
b.CData(7,:) = [0 0 1];
text(9,LMS_sdr_B_norm(1),num2str(LMS_sdr_B_norm(1),'%.1f'),'Rotation', 90);
text(10,LMS_sdr_B_norm(2),num2str(LMS_sdr_B_norm(2),'%.1f'),'Rotation', 90);
text(11,LMS_sdr_B_norm(3),num2str(LMS_sdr_B_norm(3),'%.1f'),'Rotation', 90);
b.FaceColor = 'flat';
b.CData(9,:) = [1 0 0];
b.CData(10,:) = [0 1 0];
b.CData(11,:) = [0 0 1];
text(13,LMS_sdr_W_norm(1),num2str(LMS_sdr_W_norm(1),'%.1f'),'Rotation', 90);
text(14,LMS_sdr_W_norm(2),num2str(LMS_sdr_W_norm(2),'%.1f'),'Rotation', 90);
text(15,LMS_sdr_W_norm(3),num2str(LMS_sdr_W_norm(3),'%.1f'),'Rotation', 90);
b.FaceColor = 'flat';
b.CData(13,:) = [1 0 0];
b.CData(14,:) = [0 1 0];
b.CData(15,:) = [0 0 1];
xticks([2 6 10 14]);
xticklabels({'R','G','B','W'});
axis([0 16 0 250]);
%ylabel('cd/m^2');
legend off
title('LMS for R G B (SPD)')

% Plot calculated L M S output using the matrix RGB2LMS 
subplot(2,2,4)  
LMS_sdr_W_norm_M = [1 1 1]*M_rgb2lms_sdr';  %white
LMS_sdr_R_norm_M = [1 0 0]*M_rgb2lms_sdr';
LMS_sdr_G_norm_M = [0 1 0]*M_rgb2lms_sdr';
LMS_sdr_B_norm_M = [0 0 1]*M_rgb2lms_sdr';
LMS_sdr_Greyty_norm_M = [0.5 0.5 0.5]*M_rgb2lms_sdr'; %grey
% Plot LMS for R G B of SDR display based on matrix
b=bar([1:15], [LMS_sdr_R_norm_M 0 LMS_sdr_G_norm_M 0 LMS_sdr_B_norm_M 0 LMS_sdr_W_norm_M]); hold on;
text(1,LMS_sdr_R_norm_M(1),num2str(LMS_sdr_R_norm_M(1),'%.1f'),'Rotation', 90);
text(2,LMS_sdr_R_norm_M(2),num2str(LMS_sdr_R_norm_M(2),'%.1f'),'Rotation', 90);
text(3,LMS_sdr_R_norm_M(3),num2str(LMS_sdr_R_norm_M(3),'%.1f'),'Rotation', 90);
b.FaceColor = 'flat';
b.CData(1,:) = [1 0 0];
b.CData(2,:) = [0 1 0];
b.CData(3,:) = [0 0 1];
text(5,LMS_sdr_G_norm_M(1),num2str(LMS_sdr_G_norm_M(1),'%.1f'),'Rotation', 90);
text(6,LMS_sdr_G_norm_M(2),num2str(LMS_sdr_G_norm_M(2),'%.1f'),'Rotation', 90);
text(7,LMS_sdr_G_norm_M(3),num2str(LMS_sdr_G_norm_M(3),'%.1f'),'Rotation', 90);
b.FaceColor = 'flat';
b.CData(5,:) = [1 0 0];
b.CData(6,:) = [0 1 0];
b.CData(7,:) = [0 0 1];
text(9,LMS_sdr_B_norm_M(1),num2str(LMS_sdr_B_norm_M(1),'%.1f'),'Rotation', 90);
text(10,LMS_sdr_B_norm_M(2),num2str(LMS_sdr_B_norm_M(2),'%.1f'),'Rotation', 90);
text(11,LMS_sdr_B_norm_M(3),num2str(LMS_sdr_B_norm_M(3),'%.1f'),'Rotation', 90);
b.FaceColor = 'flat';
b.CData(9,:) = [1 0 0];
b.CData(10,:) = [0 1 0];
b.CData(11,:) = [0 0 1];
text(13,LMS_sdr_W_norm_M(1),num2str(LMS_sdr_W_norm_M(1),'%.1f'),'Rotation', 90);
text(14,LMS_sdr_W_norm_M(2),num2str(LMS_sdr_W_norm_M(2),'%.1f'),'Rotation', 90);
text(15,LMS_sdr_W_norm_M(3),num2str(LMS_sdr_W_norm_M(3),'%.1f'),'Rotation', 90);
b.FaceColor = 'flat';
b.CData(13,:) = [1 0 0];
b.CData(14,:) = [0 1 0];
b.CData(15,:) = [0 0 1];
xticks([2 6 10 14]);
xticklabels({'R','G','B','W'});
axis([0 16 0 250]);
%ylabel('cd/m^2');
legend off
title('LMS for R G B (Matrix)')




