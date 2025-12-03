function b_hat = qpsk_demod(y)

    Nsym = length(y);
    b_pairs = zeros(Nsym,2);

    for i = 1:Nsym
        xr = real(y(i));
        xi = imag(y(i));

        if xr>=0 && xi>=0
            b_pairs(i,:) = [1 1];
        elseif xr<0 && xi>=0
            b_pairs(i,:) = [1 0];
        elseif xr<0 && xi<0
            b_pairs(i,:) = [0 0];
        else
            b_pairs(i,:) = [0 1];
        end
    end

    b_hat = reshape(b_pairs.',[],1);

end
