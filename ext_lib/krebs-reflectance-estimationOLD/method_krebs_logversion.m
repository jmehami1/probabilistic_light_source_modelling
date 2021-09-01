function [S, g, k] = method_krebs_logversion(I, mask)
%returns S g and k maps of the dichromatic model
%I is a multispectral image Nx*Ny*Nc
%Nx is the number of rows, Ny is the number of columns and Nc the number of
%wavelength

[Nx, Ny, Nc] = size(I);
Np = Nx*Ny; % number of pixels

%reshape the images
I_reshaped = reshape(double(I), Np, Nc);
mask_reshaped = reshape(double(mask>0), Np, 1);

% CSTE
filter_width = 11;
filter_half_width = floor(filter_width/2);
q = 3;%exp spatiale
r = (pi*pi/3600);%exp spectrale sam
r2 = (pi*pi/400);

%precompute mean and std
minI = min(I_reshaped,[],2);
meanI = mean(I_reshaped,2);
stdI = std(I_reshaped,[],2);
stdI(stdI==0) = min(stdI(stdI>0));% to avoid division by zero after

% precompute log(stdI) and log(meanI)
logstdI = log(stdI);
logstdI(~isfinite(logstdI)) = min(logstdI(isfinite(logstdI)));%avoid log(0)
logmeanI = log(meanI);
logmeanI(~isfinite(logmeanI)) = min(logmeanI(isfinite(logmeanI)));%avoid log(0)

zetaI = SAM(I_reshaped', ones(Nc, Np))';%spectral angle mapper : http://ijiset.com/v1s4/IJISET_V1_I4_27.pdf
zetaI = exp(-(zetaI.^2)/r2);%got value between 0 and 1
zetaI(~isfinite(zetaI))=0;% if NaN
zetaI = reshape(zetaI,Nx*Ny,1);
izetaI = 1-zetaI;% by definition

%prepare A % for f1
I_sparse = cell(1,filter_width*filter_width);
J_sparse = cell(1,filter_width*filter_width);
V_sparse = cell(1,filter_width*filter_width);
cur_idx = 0;
for di = -filter_half_width:filter_half_width
    for dj = -filter_half_width:filter_half_width
        cur_idx = cur_idx+1;
        
        if(di>0)
            i = 1:(Nx-di);
        else
            i = (1-di):Nx;
        end
        if(dj>0)
            j = 1:(Ny-dj);
        else
            j = (1-dj):Ny;
        end
        
        [X, Y] = meshgrid(j,i);
        P1 = (X(:)-1)*Nx+Y(:);
        P2 = P1+dj*Nx+di;
        IP1 = I_reshaped(P1,:);
        IP2 = I_reshaped(P2,:);
        diff_spe = SAM(IP2', IP1')';
        diff_spe = exp(-(diff_spe.^2)/r);
        diff_spe(~isfinite(diff_spe))=0;
        
        diff = diff_spe.*mask_reshaped(P1).*mask_reshaped(P2).*izetaI(P1).*izetaI(P2);        
        
        I_sparse{cur_idx} = P1';
        J_sparse{cur_idx} = P2';
        V_sparse{cur_idx} = diff';
    end
end
clear P1 P2 IP1 IP2 X Y i j diff_spe diff
I_sparse = cell2mat(I_sparse);
J_sparse = cell2mat(J_sparse);
V_sparse = cell2mat(V_sparse);
A = sparse(I_sparse, J_sparse, V_sparse,Np,Np);
clear I_sparse J_sparse V_sparse
%figure, hist(A(A~=0),200),title('A')
%drawnow


%prepare B % for f3
I_sparse = cell(1,filter_width*filter_width);
J_sparse = cell(1,filter_width*filter_width);
V_sparse = cell(1,filter_width*filter_width);
cur_idx = 0;
for di = -filter_half_width:filter_half_width
    for dj = -filter_half_width:filter_half_width
        cur_idx = cur_idx+1;
        
        if(di>0)
            i = 1:(Nx-di);
        else
            i = (1-di):Nx;
        end
        if(dj>0)
            j = 1:(Ny-dj);
        else
            j = (1-dj):Ny;
        end

        [X, Y] = meshgrid(j,i);
        P1 = (X(:)-1)*Nx+Y(:);
        P2 = P1+dj*Nx+di;
        IP1 = I_reshaped(P1,:);
        IP2 = I_reshaped(P2,:);
        diff_spe = SAM(IP2', IP1')';
        diff_spe = exp(-(diff_spe.^2)/r);
        diff_spe(~isfinite(diff_spe))=0;
        
        diff = diff_spe.*mask_reshaped(P1).*mask_reshaped(P2).*zetaI(P1).*zetaI(P2);
        
        I_sparse{cur_idx} = P1';
        J_sparse{cur_idx} = P2';
        V_sparse{cur_idx} = diff';
    end
end
clear P1 P2 IP1 IP2 X Y i j diff_spe diff
I_sparse = cell2mat(I_sparse);
J_sparse = cell2mat(J_sparse);
V_sparse = cell2mat(V_sparse);
B = sparse(I_sparse, J_sparse, V_sparse,Np,Np);
clear I_sparse J_sparse V_sparse


%prepare C % for f2
I_sparse = cell(1,filter_width*filter_width);
J_sparse = cell(1,filter_width*filter_width);
V_sparse = cell(1,filter_width*filter_width);
cur_idx = 0;
for di = -filter_half_width:filter_half_width
    for dj = -filter_half_width:filter_half_width
        cur_idx = cur_idx+1;
        
        if(di>0)
            i = 1:(Nx-di);
        else
            i = (1-di):Nx;
        end
        if(dj>0)
            j = 1:(Ny-dj);
        else
            j = (1-dj):Ny;
        end
        
        [X, Y] = meshgrid(j,i);
        P1 = (X(:)-1)*Nx+Y(:);
        P2 = P1+dj*Nx+di;
        IP1 = I_reshaped(P1,:);
        IP2 = I_reshaped(P2,:);
        diff_spe = SAM(IP2', IP1')';
        diff_spe = exp(-(diff_spe.^2)/r);
        diff_spe(~isfinite(diff_spe))=0;
        
        diff = mask_reshaped(P1).*mask_reshaped(P2).*(1-diff_spe.*izetaI(P1).*izetaI(P2)-diff_spe.*zetaI(P1).*zetaI(P2));
        
        I_sparse{cur_idx} = P1';
        J_sparse{cur_idx} = P2';
        V_sparse{cur_idx} = diff';
    end
end
clear P1 P2 IP1 IP2 X Y i j diff_spe diff
I_sparse = cell2mat(I_sparse);
J_sparse = cell2mat(J_sparse);
V_sparse = cell2mat(V_sparse);
C = sparse(I_sparse, J_sparse, V_sparse,Np,Np);
clear I_sparse J_sparse V_sparse

% formulation of H
%disp('formulation H')
H1 = spdiags(sum(A,2),0,Np,Np)-A;
H2 = spdiags(sum(C,2),0,Np,Np)-C;
H3 = spdiags(sum(B,2),0,Np,Np)-B;
H = H1+H2+H3;

%[v,p] = chol(H,'lower');
%[V,D] = eig(full(H));
%figure, plot(sort(diag(D))),title('D matrix aig');
%figure, imagesc(V),title('V matrix aig');
% [U,S,V] = svd(full(H));
% figure, plot(diag(S)),title('S matrix SVD');
% figure, imagesc((U.*V)>0);
% figure, plot(U(:,1)),title('U1 matrix SVD');
% figure, plot(U(:,2)),title('U2 matrix SVD');
% figure, plot(U(:,3)),title('U3 matrix SVD');
% figure, plot(U(:,4)),title('U4 matrix SVD');
% figure, plot(U(:,5)),title('U5 matrix SVD');
% figure, plot(U(:,6)),title('U6 matrix SVD');
% [V,D]=eig(full(H));
% figure, plot(sort(diag(D))),title('D matrix eig');
% figure, imagesc(full(H));
% figure, imagesc(inv(full(H)));



% formulation of f
%disp('formulation f')
f_g = H1*logstdI+H3*logmeanI;
f_k = H1*(meanI./stdI);
f_gS = H1*logstdI+(H2+H3)*logmeanI;
clear H1 H1 H3 zetaI izetaI

multiply_H = @(x) (H*x);

%debug
% [B, d] = spdiags(H);
% for i=1:size(B,2);
%     temp = B(:,i);
%     figure, imagesc(reshape(temp,Nx,Ny)),title(strcat('diag',num2str(i)));
% end
% 
% fake_img1 = reshape(1:Np,Nx,Ny);
% fake_img2 = ones(Nx,Ny);
% fake_img3 = zeros(Nx,Ny);
% fake_img3(round(Nx/2),round(Ny/2))=1;
% fake_img4 = ones(Nx,1)*[1:Ny];
% fake_img5 = [1:Ny]'*ones(1,Ny);
% figure, imagesc(reshape(H*fake_img1(:),Nx,Ny)),title('fake1');
% figure, imagesc(reshape(H*fake_img2(:),Nx,Ny)),title('fake2');
% figure, imagesc(reshape(H*fake_img3(:),Nx,Ny)),title('fake3');
% figure, imagesc(reshape(H*fake_img4(:),Nx,Ny)),title('fake4');
% figure, imagesc(reshape(H*fake_img5(:),Nx,Ny)),title('fake5');
% figure, imagesc(reshape(f_g,Nx,Ny));
% figure, imagesc(reshape(f_k,Nx,Ny));

%recover g
ginit = (ones(Nx*Ny,1)*5);% random initialisation of g
g = Conjugate_gradient(multiply_H,f_g,ginit,150);%conjugate gradient algorithm
g = exp(g);
g(~mask) = 0;

%recover k
% kinit = zeros(Nx*Ny,1);% init k with zeros
% mink = zeros(Nx*Ny,1);% lower bound acceptable for k
% maxk = minI./stdI;% upper bound acceptable for k
% k = Gradient_projection_QP(multiply_H,f_k,mink,maxk,kinit,150);% projected gradient algorithm
% k = k.*stdI;
% k(~mask) = 0;

%recover gS
gSinit = logmeanI;
mingS = meanI-minI;
mingS(mingS<=0)=0;
mingS = log(mingS);

%mingS(:)=-inf;
maxgS = logmeanI;
%gS = gSinit;
gS = Gradient_projection_QP(multiply_H,f_gS,mingS,maxgS,gSinit,150);% projected gradient algorithm
%figure, imagesc(gS)
%figure, imagesc(mingS)
%figure, imagesc(logmeanI)
%figure, imagesc((gS-mingS)<0)
%figure, imagesc((maxgS-gS)<0)

offset = min(logmeanI(:)-gS(:))
%pause()
gS = exp(gS+offset);
gS(~mask) = 0;
k = meanI-gS;
k(~mask) = 0;

%plot results
%figure, imagesc(g);
%figure, imagesc(k);
%figure, imagesc(gS);
%pause()

%reshape
S = (I_reshaped-k*ones(1,Nc))./(g*ones(1,Nc));
S(~isfinite(S))=0;
S = reshape(S, Nx, Ny, Nc);
g = reshape(g, Nx, Ny);
k = reshape(k, Nx, Ny);
%figure, imagesc(g);
%figure, imagesc(k);
%figure, imagesc(k>reshape(minI,Nx,Ny))
%figure, imagesc(reshape(minI,Nx,Ny))
%figure, imagesc(S);
%figure, plot3(I_reshaped(:,1)/(2^16),I_reshaped(:,2)/(2^16),I_reshaped(:,3)/(2^16));
%figure, imagesc(I/(2^16));
%pause()

clearvars -except g k S

end