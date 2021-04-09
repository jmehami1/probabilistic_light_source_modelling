function symmPnt = SymmetricPointLinescanLight(pnt, lightPos, lightDir, vpNormal)
%Calculate the symmetrical point about a principal axis of apoint-light source 
%where the original point is measured on the view-plane of a line-scan camera

pnt2lightPos = lightPos - pnt;

% norm(pnt2lightPos)

lightPos2SymmPnt = pnt2lightPos - 2.*(dot(lightDir, pnt2lightPos)).*lightDir;

% norm(lightPos2SymmPnt)


%symmetrical point on line-scan viewplane
symmPnt = lightPos2SymmPnt + lightPos;

% norm(vpNormal)

% %project light Position onto LS view-plane
% proj_lightPos = lightPos - dot(lightPos, vpNormal).*vpNormal;
% 
% %vector from point to projected light source position on view-plane
% proj_pnt2lightPos = proj_lightPos - pnt;
% 
% % norm(proj_pnt2lightPos)
% 
% 
% %Projection of "vector from point to light source position" onto line-scan view-plane
% %proj_pnt2lightPos = proj_pnt2lightPos - dot(proj_pnt2lightPos, vpNormal).*vpNormal;
% 
% %projection of light principal direction onto line-scan view-plane
% proj_lightDir = lightDir - dot(lightDir, vpNormal).*vpNormal;
% proj_lightDir = proj_lightDir./norm(proj_lightDir);
% 
% % norm(proj_lightDir)
% 
% % proj_lightDir = lightDir;
% 
% %vector from projected light source position to symmetrical point about projected light source
% %direction
% proj_lightPos2SymmPnt = proj_pnt2lightPos - 2.*(dot(proj_lightDir, proj_pnt2lightPos)).*proj_lightDir;
% 
% % norm(proj_lightPos2SymmPnt)
% 
% 
% %symmetrical point on line-scan viewplane
% symmPnt = proj_lightPos + proj_lightPos2SymmPnt;
% 
% norm(pnt-lightPos)
% 
% norm(symmPnt-lightPos)
% 
% dot(pnt-lightPos, lightDir)
% 
% dot(symmPnt-lightPos, lightDir)

end

