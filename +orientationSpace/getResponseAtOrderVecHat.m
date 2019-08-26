function responseAtOrder = getResponseAtOrderVecHat(response_hat,Korg,Kg)
%getResponseAtOrderVec Get response at order, vectorized so that each pixel
% has an independent order. Response already has Fourier transform applied
%
% Extracted from orientationSpace.diffusion.refineBifurcation

% Mark Kittisopikul, November 2017
% Northwestern University

    if(isempty(Kg))
%         responseAtOrder = NaN(size(response_hat));
        responseAtOrder = NaN([size(response_hat,1) numel(Kg)]);
        return;
    end
%     x = [0:ceil(R.filter.K)*R.filter.sampleFactor -ceil(R.filter.K)*R.filter.sampleFactor:-1];
%     x = [0:8 -8:-1];
    s = floor(size(response_hat,1)/2);
    x = [0:s -s:-1];
    n_org = 2*Korg+1;

    n_new = 2*Kg+1;
    s_inv = sqrt(n_org.^2.*n_new.^2./(n_org.^2-n_new.^2));
    s_hat = s_inv/(2*pi);
    f_hat = exp(-0.5 * bsxfun(@rdivide,x(:),s_hat).^2);
    responseAtOrder = bsxfun(@times,response_hat,f_hat);
end
