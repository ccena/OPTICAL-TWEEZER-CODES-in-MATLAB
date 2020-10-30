%By: Christianlly B. Cena
%This program utilized the Gerchberg-Saxton algorithm to develop holograms
%Initializing Parameters
Nx=512;%row
Ny=512;%column
L=2; %
k=0.01;
lambda=640;%wavelength
f=5;%focal length
for x = 1:Nx % Double loop to calculate the phase of each pixel
    for y = 1:Ny
        x0 = Nx/2; % locating the center of the image in the x coordinate
        y0 = Ny/2;  %locating the center of the image in the y coordinate
        xr=x-x0;
        yr=y-y0;
        r0=(xr^2+yr^2);
        r=sqrt(xr^2+yr^2); % radial coordinate
        phi=atan2(yr,xr);  % angular coordinate
        
        dpp=(L*angle(x+(i*y)))+(k*x)+((pi/(lambda*f))*r0);
        r1=mod(dpp,2*pi);
        
        r1=1-r1;
        for i=1:3
            C(y,x,i) = r1; %pixel, 3 colors
        end
        end
    end
    figure;
    imshow(C)
