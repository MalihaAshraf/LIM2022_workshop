%% Display spectra

% Load display spectra
% File structure: Column1: Wavelength (lambda), Column2: Spectrum of full
% white, Column3-5: Spectra of RGB channels respectively
spec1 = csvread( 'spectrum_sdr_display.csv' ); % Load SDR display
spec2 = csvread( 'spectrum_hdr_display.csv' ); % Load HDR display
lambda_display = spec1(:,1);

% Plot display spectra
figure, plot( lambda_display, spec1(:,2), 'DisplayName', 'SDR display');
hold on,
plot( lambda_display, spec2(:,2), 'DisplayName', 'HDR display');
xlabel('Wavelength (\lambda)');
ylabel('Display emission spectra');
legend boxoff

figure,
% Plot RGB spectra
plot( lambda_display, spec1(:,3), '--r', 'DisplayName', 'R (SDR)'); hold on
plot( lambda_display, spec2(:,3), '-r', 'DisplayName', 'R (HDR)'); hold on
plot( lambda_display, spec1(:,4), '--g', 'DisplayName', 'G (SDR)'); hold on
plot( lambda_display, spec2(:,4), '-g', 'DisplayName', 'G (HDR)'); hold on
plot( lambda_display, spec1(:,5), '--b', 'DisplayName', 'B (SDR)'); hold on
plot( lambda_display, spec2(:,5), '-b', 'DisplayName', 'B (HDR)'); hold on
legend boxoff

%% Load Cone fundamentals

% Load CIE LMS 2006 cone fundamentals. Source: cvrl.org 
cmf_lms = csvread( 'linss2_10e_1.csv' );
% Interpolate LMS response to match wavelength range of the display spectra
cmf_lms = interp1( cmf_lms(:, 1), cmf_lms(:, 2:4), lambda_display, 'linear', 0);
% Load and interpolate luminance response curve
v_lambda = csvread( 'vl_linCIE2008v2e_1.csv' );
v_lambda = interp1( v_lambda(:,1), v_lambda(:,2), lambda_display, 'linear', 0);

figure,
% Plot cone abosrption curves
plot( lambda_display, cmf_lms(:,1), '-r', 'DisplayName', 'L'); hold on
plot( lambda_display, cmf_lms(:,2), '-g', 'DisplayName', 'M');
plot( lambda_display, cmf_lms(:,3), '-b', 'DisplayName', 'S');
plot( lambda_display, v_lambda, '-k', 'DisplayName', 'Luminance');

xlim( [ 300 800 ])
xlabel('Wavelength (\lambda)');
ylabel('Normalized cone respones');
legend boxoff

%% Calculate colorimetric values

% Non-normalized LMS values
LMS_sdr = trapz(lambda_display, cmf_lms.*spec1(:,2)).*683.002; % 683.002 converts the radiometric values to cd/m^2
LMS_hdr = trapz(lambda_display, cmf_lms.*spec2(:,2)).*683.002;

% Luminance values from luminance response function
max_lum_sdr = trapz(lambda_display, v_lambda.*spec1(:,2)).*683.002;
max_lum_hdr = trapz(lambda_display, v_lambda.*spec2(:,2)).*683.002;

% LMS normalization weights
lms_weights = [0.689903 0.348322 0.0371597];
LMS_sdr_norm = lms_weights.*LMS_sdr;
LMS_hdr_norm = lms_weights.*LMS_hdr;

% Luminance values from normalized LMS values (L+M)
lum_sdr_2 = sum(LMS_sdr_norm(1:2)); % value should be equal to max_lum_sdr
lum_hdr_2 = sum(LMS_hdr_norm(1:2)); % value should be equal to max_lum_hdr


%% Get RGB2LMS matrix from spectra

% Max value of RGB [1 1 1] for the displays respectively
rgb_sdr_weights = trapz(lambda_display, spec1(:,3:5).*spec1(:,2));
rgb_hdr_weights = trapz(lambda_display, spec2(:,3:5).*spec2(:,2));

% Normalized RGB response curves of displays
RGB_sdr_response = spec1(:,3:5)./rgb_sdr_weights;
RGB_hdr_response = spec2(:,3:5)./rgb_hdr_weights;

% Calculate RGB -> LMS matrix from display RGB spectra and LMS responses
M_rgb2lms_sdr = (pinv(RGB_sdr_response)*(lms_weights.*cmf_lms.*683.002))';
M_rgb2lms_hdr = (pinv(RGB_hdr_response)*(lms_weights.*cmf_lms.*683.002))';


% Calculate LMS values from RGB

lms_test1 = [1 1 1]*M_rgb2lms_sdr';  % should be equal to LMS_sdr_norm
lms_test2 = [1 1 1]*M_rgb2lms_hdr';  % should be equal to LMS_sdr_norm

% Calculate LMS2RGB matrix
M_lms2rgb_sdr = pinv(M_rgb2lms_sdr);
M_lms2rgb_hdr = pinv(M_rgb2lms_hdr);


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
M_rgb2xyz_sdr = (pinv(RGB_sdr_response)*(cmf_xyz.*683.002))';
M_rgb2xyz_hdr = (pinv(RGB_hdr_response)*(cmf_xyz.*683.002))';
M_xyz2rgb_sdr = pinv(M_rgb2xyz_sdr);
M_xyz2rgb_hdr = pinv(M_rgb2xyz_hdr);


% Calculate LMS <-> XYZ matrix from XYZ and LMS responses
M_lms2xyz = (pinv(lms_weights.*cmf_lms.*683.002)*(cmf_xyz.*683.002))';
M_xyz2lms = pinv(M_lms2xyz);



%% Save all color transformation matrices

save('../M.mat', 'M_lms2rgb_hdr', 'M_lms2rgb_sdr', 'M_lms2xyz', 'M_rgb2lms_hdr',...
    'M_rgb2lms_sdr', 'M_rgb2xyz_hdr', 'M_rgb2xyz_sdr', 'M_xyz2lms', 'M_xyz2rgb_hdr',... 
    'M_xyz2rgb_sdr', 'max_lum_sdr', 'max_lum_hdr');

