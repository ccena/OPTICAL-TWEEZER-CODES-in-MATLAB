
function [MSD,T]  = compute_msd(y)
	% Function for computing the MSD
    [r c] = size(y);
    % N = ceil(r*1)
    N=ceil(r*1);
    for i=1:N
       
        xo = y(1:r-i,:);
        xf = y(1+i:r,:);
        % disp([i, size(xo,1), size(xf,1)])
        L =size(xo,1);
        MSD(i,:) = sum ((xf-xo).^2)/L;
    end
    T = 1:N
