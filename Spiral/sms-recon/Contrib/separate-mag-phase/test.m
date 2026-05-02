beta11 = 6e4; % regularization paramter for phase (rg2/rg4)
soft = @(t,a) (t - a * sign(t)) .* (abs(t) > a); % soft thresholding function

U = Wavelet3(Nx,Ny,Nz,mask,5);

for iter = 1:niter
    xi = pcg_bls_exp_ep(A, C, yi, mi, xi, beta11, del, ...
        nsubiter(1), 'mask', mask); % with regularizer 4
    
    phxi = exp(1j*xi);
    if firstUpdate
        tmp = img_i + 1/curv * conj(rfph) .* (A' * (yi - A * (rfph.*img_i)));
        img_i = U' * soft(U * tmp, beta22 / curv);
    else
        tmp = U * mi + 1/curv * real(U * (conj(phxi) .*(A'*(yi-A*(phxi.*mi)))));
        mi = U' * soft(tmp, beta22 / (2 * curv));
    end