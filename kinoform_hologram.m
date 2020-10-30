%By: Christianlly Cena
%..............."Hologram Creation Using a Kinoform Algorithm"..........%

clear all; close all;
%Creating a Sample Image
%I = zeros(250,250);%creates a 250x250 matrix consists of zeroes
%for n=1:10:250
 %I(1:250,n:(n+4))=1;   %loop that will generate strips
%end
I=imread('cameraman.tif', 'tif');
I=double(I);
figure;imshow(mat2gray(I))%converts the matrix to an image and displays it
title('Original Image')

%Creating a Random Phase Mask
pm=rand([256,256]);%creating a 250x250 (size of the image) array of random numbers
I=times(I,exp(-2i*pi*pm));%adding a random phase mask

%Shifting of domain from amplitudes to phases
A=fftshift(fft2(fftshift(I)));%calculates DFT (Discrete Fourier Transform) 
%fftshift places the zero-freq component to the center of the array
figure;imshow(mat2gray(abs(A)));colorbar;
title('Random Phase Mask')

%Phase of the Image Spectrum
B=angle(A);%returns the phase angles of A
figure; imshow(mat2gray(B));colorbar;
title('Phase of the Image Spectrum')

%Reconstruction (FFT)
%From Phase to Complex Amplitude
C=fftshift(ifft2(fftshift(exp(i*B))));%calculates IDFT (Inverse Discrete Fourier Transform)
%exp(i*B)is the Kinoform hologram
figure; 
imshow(mat2gray(abs(C)));
title('Reconstructed image')
