%Blazed Grating Program
%By Christianlly B. Cena
clear all;
close all;
%%
%User-defined blazed patterns
opengl hardware
disp('Choose Blazed Grating Pattern:');
disp('1. Vertical');
disp('2. Horizontal');
disp('3. Radial');
disp('4. Diagonal');
disp('5. Fork');
blaze = input('Desired Pattern: ');

%Modification and Saving of the pattern
binary=true;%if true,creates a binary image of the selected blazed pattern
ampmod=false;%if true,performs amplitude modulation
collage=false;%if true, creates an image consisting of 4 blazed gratings
FFT=true;%if true, calculates the inverse fourier transform of the desired pattern
saveimage=true;%if true, saves all the pattern in the matlab folder
savefft=false;
%%
%Definitions
A=5;%grating period at large distances away from the fork; for amp mod
fd= input('Desired Number of Fringe: ');%0.5 % This parameter defines the grating period of a blazed grating
lg= input('Desired Azimuthal Number: '); % This is the value of the topological charge of the LG beam
w= input('Input Beam waist: '); %modulation width, for amplitude modulation
imName2='opticaltrap';
imName3='collage';% 0This is the name of the file that has the collage
imName4='opticaltrapcollage';
lg1=num2str(lg);%for file naming
fd1=num2str(fd);%for file naming
A1=num2str(A);%for file naming
w1=num2str(w);%for file naming
Nx = 512; % # of pixels in x-dimension of Meadowlark SLM
Ny = 512;  % # of pixels in y-dimension of Meadowlark SLM
C = ones(Ny,Nx,3); % Initialization of the image matrix

%%
%Phase Calculation
for x = 1:Nx % Double loop to calculate the phase of each pixel
    for y = 1:Ny
        x0 = Nx/2; % locating the center of the image in the x coordinate
        y0 = Ny/2;  %locating the center of the image in the y coordinate
        xr=x-x0;
        yr=y-y0;
        r=sqrt(xr^2+yr^2); % radial coordinate
        phi=atan2(yr,xr);  % angular coordinate
        %Amplitude modulation
        if ampmod % Amplitude modulation
            %am=exp(((-pi)*(r^2))/(w^2));
            am=r^abs(lg)*exp(-r^2/w^2)/((w*sqrt(lg/2))^abs(lg)*exp(-lg/2));
        else
            am=1;
        end
        
        %Blazing
        switch blaze
            case 1
                
                imName='blazed';
                blz=fd*(xr);%phi=1, for a blazed grating
                r1 = mod(blz,2*pi)/(2*pi); %phase modulation
                imName1='vrtcl';
                if binary
                    imName='binary';
                    if r1 >=0.5
                        r1=0;
                    else
                        r1=1;
                    end
                end
                %Flipping the Matrix
                r1=r1*am;  % For Reverse scaling use r1=1-r1*am or r1=r1*am
                %Incorporating each value to the pixel (x,y) of Image C
                for i=1:3
                    C(y,x,i) = r1; %pixel, 3 colors
                end
            case 2
                
                imName='blazed';
                blz=fd*(yr);
                r1 = mod(blz,2*pi)/(2*pi);
                imName1='hzntl'; % This is the name of the file that has the pattern
                if binary
                    imName='binary';
                    if r1 >=0.5
                        r1=0;
                    else
                        r1=1;
                    end
                end
                r1=r1*am; 
                for i=1:3
                    C(y,x,i) = r1;
                end
            case 3
                
                imName='blazed';
                blz=fd*(r); %for zoomed out radial patterns change to r^2
                r1 = mod(blz,2*pi)/(2*pi);
                imName1='rdl';
                if binary
                    imName='binary';
                    if r1 >=0.5
                        r1=0;
                    else
                        r1=1;
                    end
                end
                %Flipping the Matrix
                r1=r1*am;
                for i=1:3
                    C(y,x,i) = r1; %pixel, 3 colors
                end
            case 4
                
                imName='blazed';
                blz=fd*(x+y);
                r1 = mod(blz,2*pi)/(2*pi);
                imName1='dgnl';
                if binary
                    imName='binary';
                    if r1 >=0.5
                        r1=0;
                    else
                        r1=1;
                    end
                end
                
                %Flipping the Matrix
                r1=r1*am; 
                for i=1:3
                    C(y,x,i) = r1; %pixel, 3 colors
                end
                
            case 5
                
                imName='blazed';
                %blz=((lg*phi)-((lg*2*pi)));
                blz=((lg*phi)-((2*pi)/A)*(r*cos(phi)));
                r1 = mod(blz,2*pi)/(2*pi);
                imName1='fork'; %Flipping the Matrix
                if binary
                    imName='binary';
                    if r1 >=0.5
                        r1=0;
                    else
                        r1=1;
                    end
                end
                r1=r1*am; 
                for i=1:3
                    C(y,x,i) = r1; %pixel, 3 colors
                end  
        end
    end
end
%%
%Image File Creation
B=rot90(C,2);
figure;
imshow(B)
title('Original object');
B=double(B);
B=bsxfun(@rdivide,B,max(max(B)));
axis equal;
axis off;
fileType = 'bmp';
imName1=strcat(imName, imName1);
if saveimage
    imfname1=[imName1,'fd',fd1,'lg',lg1,'A',A1,'w',w1,'.bmp']; %saves the pattern in the matlab folder
    imwrite(B, imfname1, fileType);
end
%%Fast Fourier Transform of the pattern
if FFT
    Q=B;
    Q=fftshift(ifft2(fftshift(Q)));
    Q=bsxfun(@rdivide,Q,max(max(Q)));
    figure;
    Q=abs(Q);
    imshow((Q));
    title('Optical Trap');
    if savefft
        imfname2=[imName2,imName1,'fd',fd1,'lg',lg1,'A',A1,'w',w1,'.bmp']; %saves the pattern in the matlab folder
        imwrite((Q), imfname2, fileType);
    end
end
%%
%Combination of Patterns
if collage
    Q1=imresize(C,[x0 y0]); %resizing C to half its normal size
    Q2=permute(Q1,[2 1 3]); %transposing the Q1
    Q3=rot90(Q1,2);
    Q4=rot90(Q2,2);
    Q12=cat(2,Q1,Q2); %combining Q1 and Q2
    Q43=cat(2,Q4,Q3); %combining Q4 and Q3
    Q5=cat(1,Q12,Q43); %combining Q12 and Q43
    %Q5=bsxfun(@rdivide,Q5,max(max(Q5)));
    figure;
    imshow(Q5)
    title('Collage Pattern');
    if saveimage
        imfname3=[imName3,imName1,'fd',fd1,'lg',lg1,'A',A1,'w',w1,'.bmp']; %saves the pattern in the matlab folder
        imwrite(Q5, imfname3, fileType);
    end
    
    if FFT %for fft of the collage pattern
        Q5=fftshift(ifft2(fftshift(Q5)));
        figure;
        Q5=abs(Q5);
        imshow(Q5);
        title('Optical Trap');
        if savefft
            imfname4=[imName4,imName1,'fd',fd1,'lg',lg1,'A',A1,'w',w1,'.bmp']; %saves the pattern in the matlab folder
            imwrite((Q5), imfname4, fileType);
        end
    end
end
