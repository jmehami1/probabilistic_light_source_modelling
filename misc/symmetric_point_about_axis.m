    untitled
    %find the origin of the symmetric target about the light source direction vector
    targetSymmPtsLS = zeros(3, length(vtargetPix));
    targetCentSymmPos = SymmetricPoint_LightSrc_Linescan(targetCentPos, locLightSrcLS, rotLightSrcDirLS, vpNormal);
    
    
    
    
    %calculate the symmetric normal, which is actually finding the rotation
    %of the symmetrical target and taking its z-directional axis
    %y direction vector found using edge points
    ydir = targetSymmPtsLS(:,end) - targetSymmPtsLS(:,1);
    ydir = ydir./norm(ydir);
    
    %x direction vector is same as the original
    xdir = T_target_2_LS(1:3,1);
    zdir = cross(xdir, ydir);
    
    rot = [xdir, ydir, zdir];
    targetPoseSymLS = [rot, targetCentSymmPos; 0, 0, 0, 1];
    targetSymmPoseFrame = T_LS_2_F*targetPoseSymLS;
    
    targetSymmPtsFrame = T_LS_2_F*[targetSymmPtsLS; ones(1,length(targetSymmPtsLS(1,:)))];
    targetSymmPtsFrame = targetSymmPtsFrame(1:3,:);
    targetSampleSymmNormal(:,sample) = targetSymmPoseFrame(1:3,3);
    
