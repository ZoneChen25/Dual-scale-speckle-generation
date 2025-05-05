%% Utilizing the combination of Jinc Function to generate dual-scale Turing speckle
% The code is used for the generation of dual-scale Turing speckle with spatial kernel.
% 1. Section "Parameter" shows different parameters for users to modify. It
% should be noted that the the size of the Kernel and the generated image
% should be odd number.
% 2. The recommended number of convolutions is 5.
% 3. Try to change C_extra to change the duty factor.If you find it hard to 
% increase duty factor, try to change Convsize. The recommended duty factor
% is 0.5. 
clear;
close all;
clc;

%% Parameter
Scale_ratio = 5;  % The scale ratio between dual-scale speckles
Convsize = 221;  % The size of the Kernel(odd number)
Convnum = 5;  % Number of convolutions
D1 = 3;  % The target near-filed speckle diameter(pixel) after imaging
D2 = 3;  % The target far-filed speckle diameter(pixel) after imaging
Target_ratio = 1;  % Amplitude ratio between far-field and near-field speckle
C_extra = 80;  % The lower one between the amplitudes of two bands 
Far_field_region = 1200;  % Far_field_region(pixel) = Print_size / dpp; Print_size is the pratical print size(mm); dpp is pixel-to-physical size ratio(mm/pixel)
Rows = 7087;  % The height of the generated speckle image(odd number)
Cols = Rows;  % The width of the generated speckle image(odd number)


%% Generate a Gaussian white noise image
% Create an all-255 image
White_image = ones(Rows, Cols, 'uint8'); 

% Create a Gaussian white noise image
Sigma = 100;
Noise = uint8(randn(Rows, Cols) * Sigma);
White_noise_image = White_image + Noise;

% Transform the image from uint8 to double
Img_o = im2double(White_noise_image);


%% Initialization
% Initialization of Jinc Function
C = 0.001;  %  Parameter c = (a2 - a1)/2 = (a4 - a3)/2; Determine the target band width
Convspan = 4;
u = linspace(-Convspan, Convspan, Convsize);  % The value ranges of Jinc Function
v = u;
[U, V] = meshgrid(u, v);  % Build a mesh between [-Convspan, Convspan]


% Window Function
Window = blackman(length(u));  % Blackman window function
W_2d = Window * Window';  % 2D Window
Rho = sqrt(U.^2 + V.^2);
Rho(Rho == 0) = eps;  % Avoid denominator to be 0


% Initialization of Dual-band Function
C_far = Far_field_region * (Convsize - 1) / 2 / Convspan / Rows;  % The far-field coefficient
C_near = C_far * Scale_ratio;  % The near-field coefficient
A12 = C_far / D2;
A1 = A12 - C;
A2 = A12 + C;
F_uv1 = A1 * besselj(1, A1 * 2 * pi * Rho) ./ Rho .* W_2d;
F_uv2 = A2 * besselj(1, A2 * 2 * pi * Rho) ./ Rho .* W_2d;
A34 = C_near / D1;
A3 = A34 - C;
A4 = A34 + C;
F_uv3 = A3 * besselj(1, A3 * 2 * pi * Rho) ./ Rho .* W_2d;
F_uv4 = A4 * besselj(1, A4 * 2 * pi * Rho) ./ Rho .* W_2d;


%% Build a dual-scale Kernel
% The initial Kernel before self-convolution
Weight_a = Target_ratio; 
Kernel_initial = Weight_a ^ (1 / Convnum) * (F_uv2 - F_uv1) + (F_uv4 - F_uv3);
Span_half = (Rows - 1) / 2;
Freq_axis = -Span_half : Span_half;

% Loop convolution
Kernel_temp = Kernel_initial;
if Convnum > 1    % Choose to self-conv for reducing computing time
    for ii = 1 : (Convnum - 1)
        Kernel_temp = conv2(Kernel_temp, Kernel_initial, 'same');
    end
end

% Adjust the coefficient to ensure the Target_ratio to retain the set value after convolution
Kernel_fft = abs(fftshift(fft2(Kernel_temp, Rows, Rows)));
Kernel_fft_cut = Kernel_fft(Span_half + 1, :);
Peaks = findpeaks(Kernel_fft_cut, 1 : Rows); 
Peaks_sort = sort(Peaks, 'descend');
if Target_ratio < 1
    Low_freq_peak = Peaks_sort(3);
    High_freq_peak = Peaks_sort(1);
else
    Low_freq_peak = Peaks_sort(1);
    High_freq_peak = Peaks_sort(3);
end
Current_ratio = Low_freq_peak / High_freq_peak;
Weight_a = (Target_ratio / Current_ratio) ^ (1 / Convnum);
Kernel_adjust = Target_ratio .^ (1 / Convnum) .* Weight_a .* (F_uv2 - F_uv1) + (F_uv4 - F_uv3);


Kernel_end = Kernel_adjust;
if Convnum > 1
    for ii = 1 : (Convnum - 1)
        Kernel_end = conv2(Kernel_end, Kernel_adjust, 'same');
    end
end

% Ensure the lower one between the amplitudes of two bands to be C_extra
Kernel_end_fft = abs(fftshift(fft2(Kernel_end, Rows, Rows)));
Kernel_fft_cut1 = Kernel_end_fft(Span_half + 1, :);
Peaks1 = findpeaks(Kernel_fft_cut1, 1 : Rows); 
Peaks_sort1 = sort(Peaks1, 'descend');
Kernel_end = C_extra * Kernel_end / Peaks_sort1(3);  % The final Kernel

% DFT of the kernel
Kf_final = abs(fftshift(fft2(Kernel_end, Rows, Rows)));
figure, plot(Freq_axis, Kf_final(Span_half + 1, :), 'LineWidth', 1.5);
title('The spectrum of the Kernel')


% The dual-scale Turing Speckle
DST_speckle = imfilter(Img_o, Kernel_end, 'same', 'replicate', 'conv');  % Apply the kernel
DST_speckle_2 = imbinarize(DST_speckle);  % Transform to binary image
Iuint8 = uint8(DST_speckle_2 * 255);  % Transform to uint8 image

% The duty factor of the generated speckle
Duty_factor = sum(DST_speckle_2(:)) / Rows / Rows;  % Number of 1 / Total number
display(['The duty factor of the generated speckle is ', num2str(Duty_factor)]);

% Display
figure, imshow(DST_speckle_2);
title('Binary Dual-scale Speckle')

% Save the image
% imwrite(Iuint8, 'Dual_scale_speckle');

