%% Utilizing frequency domain filtering method to generate dual-scale speckle
% The code is used for generating dual-scale speckles with a dual-band
% filter. Section "Parameter" shows the variables for users to modify.
clear;
close all;
clc;

%% Parameter
% The size of the generated image
rows = 500;
cols = 500;

% Define the parameters of the filter 
D0_low = 30;  % Center radius of low frequency band
D0_high = 150;  % Center radius of high frequency band
width_low = 2;  % Width of low frequency band
width_high = 2;  % Width of high frequency band
H_low = 1.02;  % Coefficient of low frequency band
H_high = 1;  % Coefficient of high frequency band
num = 5;  % Filter cycle times

%% Generate a Gaussian white noise image
% Create an all-255 image
white_image = ones(rows, cols, 'uint8'); 

% Create a Gaussian white noise image
sigma = 100;
noise = uint8(randn(rows, cols) * sigma);
white_noise_image = white_image + noise;

% Transform the image from uint8 to double
I = im2double(white_noise_image); 

% Convert to the frequency domain
F = fft2(I);

% Move the zero frequency component to the spectrum center
F_shifted = fftshift(F); 

%% Create a dual-scale image
% Build a dual-band filter
[M, N] = size(F);
[u, v] = meshgrid(-N/2:N/2-1, -M/2:M/2-1);
D = sqrt(u.^2 + v.^2);
H = zeros(size(F));
H(D <= D0_low + width_low/2 & D >= D0_low - width_low/2) = H_low;  % The low frequency band
H(D <= D0_high + width_high/2 & D >= D0_high - width_high/2) = H_high;  % The high frequency band

% Apply the filter
G_shifted = H .* F_shifted;

% Loop the filter
for i= 1 : num 
    G_shifted = H .* G_shifted;
end

% Convert to the spatial domain
G = ifftshift(G_shifted);
J = real(ifft2(G));

% Display images before and after filtering
figure
imshow(I);
title('Original Image');

figure
imshow(J, []);
title('Filtered Image');

% Transform to uint8 image
Iuint8 = uint8(mat2gray(J)*255);

% Tranform to binary image if necessary
% I_bi = imbinarize(Iuint8); 

% Save the image
% imwrite(Iuint8, '1.tif');
