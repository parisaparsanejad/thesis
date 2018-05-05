% CORONAL image Processing
close all;
clear all;
clc
I=imread('CORONAL.jpg');
imshow(I);
I=imcrop(I);
imshow(I)

%% Filtering Steps
%Use double precision and scale image values to the range of [0 1]
I = im2double(I); 
I=I./(max(max(I)));

%2D Gaussian filter
I_Gau=imgaussfilt(I,2);
figure;imshow(I_Gau);title('After Gaussian filtering')

% Adjust image contrast (stretch values)
I_Con= imadjust(I_Gau,[0.4 0.9],[]);
h=figure('Position',[100 100 1000 1000]);hold on; %Normal plot
imshow(I_Con);title('After image adjustment')

%Thresholding (create binary image)
IBi= im2bw(I_Con, 0.25);
figure;imshow(IBi);title('After threshold (binary)')

[x,y] = ginput(1)
figure
imshow(IBi)
hold on
plot(x,y, 'b*')
hold off


A = [x 500]; 
B = [y y]; 

figure
imshow(IBi)
hold on
plot(A,B,'*')
figure
imshow(IBi)
hold on
line(A,B)
hold off



