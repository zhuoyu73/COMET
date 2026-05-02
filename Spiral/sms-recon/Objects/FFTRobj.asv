% test 2 with soft-threhsolding and Rbeta = 1
classdef FFTRobj
    properties
        beta
        mask
    end

    methods
        function obj = FFTRobj(mask, varargin)
            if ~islogical(mask)
                mask = logical(mask);
            end
            p = inputParser;
            addOptional(p, 'beta', 1);
            parse(p, varargin{:});
            obj.beta = p.Results.beta;
            obj.mask = mask;
        end

        function out = penal(obj, x)
            x_img = embed(x, obj.mask);
            Xf = fft2(x_img);
            out = obj.beta * sum(abs(Xf(:)));
        end

        function out = cgrad(obj, x)
            x_img = embed(x, obj.mask);
            Xf = fft2(x_img);
            Xf_thresh = soft_thresh(Xf, obj.beta);
            x_rec = ifft2(Xf_thresh);
            grad_img = x_img - x_rec;
            out = col(real(grad_img(obj.mask)));
        end

        function out = denom(obj, ddir, x)
            dx_img = embed(ddir, obj.mask);
            x_img = embed(x, obj.mask);
            FX = fft2(x_img);
            Fdx = fft2(dx_img);
            W = 1 ./ (abs(FX) + eps);
            out = obj.beta * sum((W(:).^2) .* abs(Fdx(:)).^2);
        end
    end
end

function out = soft_thresh(x, lambda)
    out = max(abs(x) - lambda, 0) .* sign(x);
end



% test 1 with Rbeta = 1e-3
%{
classdef FFTRobj
    properties
        beta
        mask
    end

    methods
        function obj = FFTRobj(mask, varargin)
            if ~islogical(mask)
                mask = logical(mask);
            end
            p = inputParser;
            addOptional(p, 'beta', 1);
            parse(p, varargin{:});
            obj.beta = p.Results.beta;
            obj.mask = mask;
        end

        function out = cgrad(obj, x)
            x_img = embed(x, obj.mask);
            Xf = fft2(x_img);
            grad_Xf = soft_thresh_grad(Xf, obj.beta);  % subgradient
            x_rec = ifft2(grad_Xf);
            x_rec = real(x_rec);
            out = col(x_rec(obj.mask));
        end
        
        function out = penal(obj, x)
            x_img = embed(x, obj.mask);
            Xf = fft2(x_img);
            out = obj.beta * sum(abs(Xf(:)));  % L1 norm
        end
        
        function out = denom(obj, ddir, x)
            dx_img = embed(ddir, obj.mask);
            x_img = embed(x, obj.mask);
            FX = fft2(x_img);
            Fdx = fft2(dx_img);
            W = 1 ./ (abs(FX) + eps);  % adaptive weighting
            out = obj.beta * sum((W(:).^2) .* abs(Fdx(:)).^2);
        end


    end
end

function X_grad = soft_thresh_grad(X, beta)
    % Subgradient: sign(X) if |X| > beta, 0 otherwise
    X_grad = sign(X) .* (abs(X) > beta);
end
%}