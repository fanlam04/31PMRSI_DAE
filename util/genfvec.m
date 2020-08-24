function fvec = genfvec(N,BW)
    % Generates frequency values for the DFT of data with length N and bandwidth
    % BW.
    n1 = -floor(N/2);
    n2 = floor((N+1)/2)-1;
    fvec = (n1:n2)*BW/N;
end

