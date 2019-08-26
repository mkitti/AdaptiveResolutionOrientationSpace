function vq = horner_vec_real_freq_simpler(v_h,xq)
    % v represents the coefficients of the polynomial
    %   D x N
    %   D = degree of the polynomial - 1
    %   N = number of polynomials
    % xq represents the query points
    %   Q x N
    %   Q = number of query points per polynomial
    %   N = number of polynomials
    % vq will be a Q x N matrix of the value of each polynomial
    %    evaluated at Q query points
    %
    % SIMPLIFICATIONS
    % v_h is a 2-D Matrix   D x N
    % xq is a 2-D Matrix    Q x N
    
%     org_sz = size(xq);
%     v_h = v_h(:,:);
%     xq  = xq(:,:);

    nQ = size(xq,1);
    
    s = size(v_h);
    scale_factor = s(1);

    % Calculate fft and nyquist frequency
    nyquist = ceil((s(1)+1)/2);
    
    % If there is an even number of fourier coefficients, split the nyquist frequency
    if(~rem(s(1),2))
        % even number of coefficients
        % split nyquist frequency
        v_h(nyquist,:) = real(v_h(nyquist,:))/2;
%        v_h = v_h([1:nyquist nyquist nyquist+1:end],:);
%        v_h = reshape(v_h,[s(1)+1 s(2:end)]);
    end
    
    % z is Q x N
    z = exp(1i*pi*xq);
    % vq starts as 1 x N
%     colon = {':'};
%     v_h_colon = colon(ones(ndims(v_h)-1,1));
       
    vq = v_h(nyquist,:);
    vq = repmat(vq,nQ,1);
    for j = nyquist-1:-1:2
        vq = z.*vq;
        vq = repmat(v_h(j,:),nQ,1)+vq;
%         vq = bsxfun(@times,z,vq);
%         vq = bsxfun(@plus,v_h(j,v_h_colon{:}),vq);
    end
       
    % Last multiplication
%     vq = bsxfun(@times,z,vq); % We only care about the real part
    vq = z.*vq;
    vq = real(vq);
%     vq = bsxfun(@times,real(z),real(vq))-bsxfun(@times,imag(z),imag(vq));
    % Add Constant Term and Scale
%     vq = bsxfun(@plus,v_h(1,v_h_colon{:}),vq*2);
    vq = repmat(v_h(1,:),nQ,1)+vq*2;
%     vq = real(vq); % We already selected the real part above
    vq = vq./scale_factor;
    
%     vq = reshape(vq,org_sz);
end
