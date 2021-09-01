function X = Gradient_projection_QP(multiply_Hx,f,lb,ub,X0,cpt_max)
%Gradient_projection_QP find the minimum of 0.5*x'*H*x-f'*x for x subject to
%simple bounds
% 
% subject to lb<=x<=ub
% H must be symmetric
% Arguments : X0 is the starting point
%             lb and ub are respectively the lower and the upper bounds acceptable for X
%             multiply_Hx is a function that returns the product Hx
%             f is a matrix of the same size as X0
%             cpt_max is the maximum number of iteration allowed
% Example : minimize x^2-2*x between 10 and 20
%           Gradient_projection_QP(@(x)(2*x),2,10,20,15)

X = X0;
cpt = 0;
X = min(X,ub);
X = max(X,lb);
while(cpt<cpt_max)
    
    Hx = feval(multiply_Hx, X);
    grad = Hx-f;%compute gradient
    p = -grad;%direction negative gradient
    
    %identify what constraints are active
    at_lb = (X-lb).^2<10e-12;
    at_ub = (X-ub).^2<10e-12;
    
    % force the direction p to be zero on active constraints
    p(at_lb & p<0) = 0;
    p(at_ub & p>0) = 0;
    p(abs(p)<10e-12) = 0;
    
    %binary masks
    positive_p = p>10e-12;
    negative_p = p<-10e-12;
    
    %if no acceptable direction found
    if(sum(p(:).^2)<10e-12)
        break;
    end
    
    %compute ideal alpha for he update x <- x+alpha*p
    H_p = feval(multiply_Hx, p);
    hpp = sum(p(:).*H_p(:));
    if(abs(hpp)<10e-12)
        break;
    end
    ideal_alpha = -sum(grad(:).*p(:))/hpp;% here need only pixelwise multiplication and division
    
    %determine minimum and maximum acceptable values of alpha 
    temp1 = (ub-X)./p;
    temp2 = (lb-X)./p;
    alpha_max = inf(size(X,1),size(X,2));
    alpha_min = -inf(size(X,1),size(X,2));
    alpha_max(positive_p) = temp1(positive_p);
    alpha_max(negative_p) = temp2(negative_p);
    alpha_min(positive_p) = temp2(positive_p);
    alpha_min(negative_p) = temp1(negative_p);

    %final alpha
    alpha = max(min(ideal_alpha,min(alpha_max(:))),max(alpha_min(:)));
    
    %update
    X = X+alpha*p;
    cpt = cpt+1;
    
end