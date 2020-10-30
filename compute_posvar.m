function [POSVAR]  = compute_posvar(y)
    %The function computes the variance of trap positions
    [r c] = size(y);
    xeq=mean(y);
    N=ceil(r*1);
    for i=1:N
       
%       xo= y(1:r-i,:);
        xf = y(1+i:r,:);
        %disp([i, size(xo,1), size(xf,1)])
        L =size(xf,1);
        POSVAR(i,:) = sum((xf-xeq).^2)/L;
    end
   
