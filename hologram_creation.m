    %Code for generating point traps using various algorithms
%By: Christianlly B. Cena
clc
close all
clear
disp('Choose algorithm:');
disp('1. Random mask (RM)');
disp('2. Superposition of gratings and lenses (S)');
disp('3. Random Superposition (SR)');
disp('4. Gerchberg-Saxton (GS)');
disp('5. Weighted Gerchberg-Saxton (GSW)');
alg = input('Number of algorithm: ');
T=input('Number of traps:');
spacing=input('Trap Spacing:');
saveimage=false;
%%
%Definition of Parameters
f = 2e-3; %focal distance of objective lens
s = 1; %SLM pixel scaling factor
pfoc = 8e-8; %arbitrary pixel width of CCD in focal plane
pslm = 15e-6*s; %arbitrary pixel width of SLM
lambda = 1064e-9; %wavelength of the laser beam
M = 512; %SLM pixel rows
N = 512; %SLM pixel columns
P = 256; %number of available discrete phase values of SLM pixels
ccdnumpixrows = 800; %arbitrary CCD pixel rows
ccdnumpixcols = 600; %arbitrary CCD pixel columns
viewrows = 800; %viewport rows 
viewcols = 600; %viewport columns 
viewrows_start = (ccdnumpixrows-viewrows)/2; %w.r.t. CCD 
viewcols_start = (ccdnumpixcols-viewcols)/2; %w.r.t. CCD
GS_maxit = 30; %maximum number of iterations for GS algorithm
GSW_maxit = 30; %maximum number of iterations for GSW algorithm
%%
%trap definition

traps_viewport = zeros(T,4); %id | row[pixels] | column[pixels] | z[um]
for i = 1:T %Target trapping pattern
row = (mod(i-1,10)*spacing)+spacing; %(mod(i-1,10)*50)+50;
col = (floor((i-1)/5)+1)*spacing;
traps_viewport(i,:) = [i row col 0]; %row and column w.r.t. viewport
end
traps = traps_viewport;
traps(:,2) = viewrows_start + traps_viewport(:,2)-ccdnumpixrows/2;
traps(:,3) = viewcols_start + traps_viewport(:,3)-ccdnumpixcols/2;
%lookup matrices definition
slmx = zeros(M,N); %contains x coordinate [m] for each SLM pixel
slmy = zeros(M,N); %contains y coordinate [m] for each SLM pixel
for i = 1:M
for j = 1:N
slmx(i,j) = (i-M/2)*pslm;
slmy(i,j) = (j-N/2)*pslm;
end
end
ccdx = zeros(ccdnumpixrows,ccdnumpixcols); %contains x co. [m] in foc. plane for each CCD ...
%pixel
ccdy = zeros(ccdnumpixrows,ccdnumpixcols); %contains y co. [m] in foc. plane for each CCD ...
%pixel
for i = 1:ccdnumpixrows
for j = 1:ccdnumpixcols
ccdx(i,j) = (i-ccdnumpixrows/2)*pfoc;
ccdy(i,j) = (j-ccdnumpixcols/2)*pfoc;
end
end

%%
%Algorithms for creating point traps
switch alg
case 1
    %random mask (RM)
    disp('RM algorithm');
    tic;
    randmat = ceil(rand(M,N)*T);%creating a random mask
    phi = 2*pi/lambda/f*(slmx.*reshape(traps(randmat,2),M,N)*pfoc ...
        +slmy.*reshape(traps(randmat,3),M,N)*pfoc) + ...
        pi*reshape(traps(randmat,4),M,N)*10^(-6)/lambda/f^2*(slmx.^2+slmy.^2);
    phi = mod(phi,2*pi);
    phi_slm = round(phi/(2*pi/P))*2*pi/P;
    calctime = toc;
    figure;
    imshow(mat2gray(phi_slm));
    B=double(phi_slm);
    B=bsxfun(@rdivide,B,max(max(B)));
    Q=fftshift(ifft2(fftshift(phi_slm)));
    figure;
    Q=abs(Q);
    imshow((Q));
    imName1='RM'
    fileType = 'bmp';
    if saveimage
    imfname1=[imName1,'.bmp']; %saves the pattern in the matlab folder
    imwrite(B, imfname1, fileType);
    end
%superposition of gratings and lenses (S)
case 2
    disp('S algorithm');
    tic;
    deltatrap = zeros(M,N,T);
    for j = 1:T
        deltatrap(:,:,j) = 2*pi/lambda/f*(slmx.*traps(j,2)*pfoc ...
            +slmy.*traps(j,3)*pfoc) + pi*traps(j,4)*10^(-6)/lambda/f^2*(slmx.^2+slmy.^2);
    end
    phi = mod(angle(sum(exp(1i.*deltatrap),3)),2*pi);
    phi_slm = round(phi/(2*pi/P))*2*pi/P;
    calctime = toc;
    figure;
    imshow(mat2gray(phi_slm));
    B=double(phi_slm);
    B=bsxfun(@rdivide,B,max(max(B)));
    Q=fftshift(ifft2(fftshift(phi_slm)));
    figure;
    Q=abs(Q);
    imshow((Q));
    imName2='S'
    fileType = 'bmp';
    if saveimage
    imfname2=[imName2,'.bmp']; %saves the pattern in the matlab folder
    imwrite(B, imfname2, fileType);
    end
case 3
    %random superposition (SR)
    disp('SR algorithm');
    tic;
    deltatrap = zeros(M,N,T);
    theta = rand(1,T)*2*pi;
    for j = 1:T
        deltatrap(:,:,j) = 2*pi/lambda/f*(slmx.*traps(j,2)*pfoc ...
            +slmy.*traps(j,3)*pfoc) + ...
            pi*traps(j,4)*10^(-6)/lambda/f^2*(slmx.^2+slmy.^2) -ones(M,N)*theta(j);
    end
    phi = mod(angle(sum(exp(1i.*deltatrap),3)),2*pi);
    phi_slm = round(phi/(2*pi/P))*2*pi/P;
    calctime = toc;
    figure;
    imshow(mat2gray(phi_slm));
    B=double(phi_slm);
    B=bsxfun(@rdivide,B,max(max(B)));
    Q=fftshift(ifft2(fftshift(phi_slm)));
    figure;
    Q=abs(Q);
    imshow((Q))
    imName3='SR'
    fileType = 'bmp';
    if saveimage
    imfname3=[imName3,'.bmp']; %saves the pattern in the matlab folder
    imwrite(B, imfname3, fileType);
    end;
case 4
    %Gerchberg?-Saxton (GS)
    disp('GS algorithm');
    tic;
    deltatrap = zeros(M,N,T);
    deltatrap_theta = zeros(M,N,T);
    Vtraps_it = zeros(M,N,T);
    theta = rand(1,T)*2*pi;
    for j = 1:T
        deltatrap(:,:,j) = 2*pi/lambda/f*(slmx.*traps(j,2)*pfoc ...
            +slmy.*traps(j,3)*pfoc) + pi*traps(j,4)*10^(-6)/lambda/f^2*(slmx.^2+slmy.^2);
        deltatrap_theta(:,:,j) = deltatrap(:,:,j)-ones(M,N)*theta(j);
    end
    phi = angle(sum(exp(1i.*deltatrap_theta),3));
    for j = 1:T
        Vtraps_it(:,:,j) = 1/N/M*sum(sum(exp(1i *(phi-deltatrap(:,:,j))),1),2)*ones(M,N);
    end
    disp(' SR guess created');
    for k = 1:GS_maxit
        disp([' iteration ' num2str(k) '/' num2str(GS_maxit)]);
        phi = angle(sum(exp(1i.*deltatrap).*Vtraps_it ...
            ./sqrt(Vtraps_it.*conj(Vtraps_it)),3));
        for j = 1:T
            Vtraps_it(:,:,j) = 1/N/M*sum(sum(exp(1i ...
                *(phi-deltatrap(:,:,j))),1),2)*ones(M,N);
        end
    end
    
    phi = mod(phi,2*pi);
    phi_slm = round(phi/(2*pi/P))*2*pi/P;
    calctime = toc;
    figure;
    imshow(mat2gray(phi_slm));
    B=double(phi_slm);
    B=bsxfun(@rdivide,B,max(max(B)));
    Q=fftshift(ifft2(fftshift(phi_slm)));
    figure;
    Q=abs(Q);
    imshow((Q));
    imName4='GS'
    fileType = 'bmp';
    if saveimage
    imfname4=[imName4,'.bmp']; %saves the pattern in the matlab folder
    imwrite(B, imfname4, fileType);
    end
    
case 5
    disp('GSW algorithm');
    tic;
    deltatrap = zeros(M,N,T);
    deltatrap_theta = zeros(M,N,T);
    Vtraps_it = zeros(M,N,T);
    theta = rand(1,T)*2*pi;
    for j = 1:T
        deltatrap(:,:,j) = 2*pi/lambda/f*(slmx.*traps(j,2)*pfoc ...
            +slmy.*traps(j,3)*pfoc) + pi*traps(j,4)*10^(-6)/lambda/f^2*(slmx.^2+slmy.^2);
        deltatrap_theta(:,:,j) = deltatrap(:,:,j)-ones(M,N)*theta(j);
    end
    phi = angle(sum(exp(1i.*deltatrap_theta),3));
    for j = 1:T
        Vtraps_it(:,:,j) = 1/N/M*sum(sum(exp(1i *(phi-deltatrap(:,:,j))),1),2)*ones(M,N);
    end
    w_it = ones(M,N,T);
    for k = 1:GSW_maxit
        disp([' iteration ' num2str(k) '/' num2str(GSW_maxit)]);
        Vtraps_it_mod = sqrt(Vtraps_it.*conj(Vtraps_it));
        w_it = w_it*mean(Vtraps_it_mod(1,1,:),3)./Vtraps_it_mod;
        phi = angle(sum(exp(1i.*deltatrap).*w_it.*Vtraps_it./Vtraps_it_mod,3));
        for j = 1:T
            Vtraps_it(:,:,j) = 1/N/M*sum(sum(exp(1i ...
                *(phi-deltatrap(:,:,j))),1),2)*ones(M,N);
        end
    end
    phi = mod(phi,2*pi);
    phi_slm = round(phi/(2*pi/P))*2*pi/P;
    calctime = toc;
    figure;
    imshow(mat2gray(phi_slm));
    B=double(phi_slm);
    B=bsxfun(@rdivide,B,max(max(B)));
    Q=fftshift(ifft2(fftshift(phi_slm)));
    figure;
    Q=abs(Q);
    imshow((Q));
    imName5='GSW'
    fileType = 'bmp';
    if saveimage
    imfname5=[imName5,'.bmp']; %saves the pattern in the matlab folder
    imwrite(B, imfname5, fileType);
    end
end

disp('Numerical test of algorithm');
deltatraptest = zeros(M,N,T);
for j = 1:T
deltatraptest(:,:,j) = 2*pi/lambda/f*(slmx.*traps(j,2)*pfoc+ slmy.*traps(j,3)*pfoc) + ...
pi*traps(j,4)*10^(-6)/lambda/f^2*(slmx.^2+slmy.^2);
end
Vtraps = zeros(1,T);
for m = 1:T
Vtraps(m) = 1/N/M*sum(sum(exp(1i*(phi_slm-deltatraptest(:,:,m))),1),2);
end
Itraps = Vtraps.*conj(Vtraps);
e = sum(Itraps)*100;
u = 1-(max(Itraps)-min(Itraps))/(max(Itraps)+min(Itraps));
sigma = sqrt(mean((Itraps-mean(Itraps)).^2))/mean(Itraps);
disp([' calculation time = ' num2str(calctime) 's']);
disp([' efficiency = ' num2str(e)]);
disp([' uniformity = ' num2str(u)]);
disp([' standard deviation = ' num2str(sigma)]);
