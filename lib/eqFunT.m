%u = EQFUNT(k, [t, w])
%       equlibirium function u = @(w, t), 'w' is CPL mass fraction in aqueous phase, 'u' in organic, 't' is temperature [°C]
%       it is a 3rd grade polynomial with zero intercept for 'w', which is multiplied by linear function for 't', giving a total of six parameters:
%              u = (a*w^3 + b*w^2 + c*w) * (d*t + e)
%                = a*t*w^3 + b*w^3 + c*t*w^2 + d*w^2 + e*t*w + f*w
%       the [t, w] matrix can have multiple rows, in that case 'u' will be a column vector
%       'k' is vector of five parameters, which are then scaled

function u = eqFunT(k, xdata)
	t = xdata(:,1);
	w = xdata(:,2);
	k = k .* [1e3, 1e5, 10, 1e3, 1, 1]; % scaling factors of parameters

	u = k(1)*(t.*w.^3) + k(2)*(w.^3) + k(3)*(t.*w.^2) + k(4)*(w.^2) + k(5)*(t.*w) + k(6)*(w);
end
