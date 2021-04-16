function poseSym = SymmetricTargetPose_LightSrc_Linescan(targetPose, lPos, lDir, symTarPos, vpNorm)
%Find the symmetrical pose of the target about the light source direction
%vector.

tarPosLS = targetPose(1:3, 4);
dirTarZ = targetPose(1:3, 3); %this is surface normal
dirTarX = targetPose(1:3, 1);
% dirTarY = targetPose(1:3, 2);


%tip of the surface normal vector on target from line-scan camera
tipTarZ_LS = tarPosLS + dirTarZ;

%tip of the direction vector in X direction from line-scan camera
% tipTarX_LS = tarPosLS + dirTarX;

%tip of the direction vector in Y direction from line-scan camera
% tipTarY_LS = tarPosLS + dirTarY;

%incident vector from target origin to source
incid = lPos - tarPosLS;

%reflected vector from source to symmetric point
reflect = symTarPos - lPos;

%tip of the surface normal vector to source
tipTarZ_S = lPos - tipTarZ_LS;

%tip of the direction vector in X direction to source
% tipTarX_S = lPos - tipTarX_LS;

%tip of the direction vector in X direction to source
% tipTarY_S = lPos - tipTarY_LS;


rotmZ = FindSymmetricDirectionVector(lDir, incid, tipTarZ_S, reflect, vpNorm);

rotmY = cross(rotmZ, dirTarX);

rot = [dirTarX, rotmY, rotmZ];

poseSym = [rot, symTarPos; 0, 0, 0, 1];

end

