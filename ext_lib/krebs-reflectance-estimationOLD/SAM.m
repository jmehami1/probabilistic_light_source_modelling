function alpha = SAM(X, Y)
%SAM compute the spectral angle mapper
%http://www.isprs.org/proceedings/XXXV/congress/comm4/papers/432.pdf for
%example

alpha = acos(sum(X.*Y,1)./sqrt(sum(X.^2,1).*sum(Y.^2,1)));

end

