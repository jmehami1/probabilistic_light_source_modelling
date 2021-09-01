function X = Conjugate_gradient(multiply_Hx,f,X0,cpt_max)
%Conjugate_gradient solves H*x-f = 0 for x
% https://en.wikipedia.org/wiki/Conjugate_gradient_method
%
% H must be symmetric
% Arguments : X0 is the starting point
%             multiply_Hx is a function that returns the product Hx
%             f is a matrix of the same size as X0
%             cpt_max is the maximum number of iteration allowed

X = X0;
Hx = feval(multiply_Hx, X);
r = Hx-f;
p = -r;
cpt = 0;

while(cpt<cpt_max)
    Hp = feval(multiply_Hx, p);
    sumr2 = sum(r(:).^2);
    hpp = sum(Hp(:).*p(:));
    if(abs(hpp)<10e-10)
        break;
    end
    if(hpp<0.0)
        disp('Alert hessian is not positive definite')
    end
    alpha = sumr2/hpp;
    X = X+alpha*p;
    r = r+alpha*Hp;
    
    if(abs(sumr2)<10e-10)
        break;
    end
    
    beta = sum(r(:).^2)/sumr2;
    p = -r+beta*p;
    cpt = cpt+1;
end