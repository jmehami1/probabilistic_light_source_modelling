classdef LightSimulator
    %Class used to simulate a near-field non-isotropic disk light source
    %which is placed an an environment relative to a frame camera.
    
    properties(Access = private)
        ligSourLoc
        ligSourOrien
        ligSourDirVec
        maxRadiantInt
        mu
        rAtt
    end
    
    methods
        function obj = LightSimulator(ligSourLoc_,ligSourOrien_, maxRadiantInt_, mu_, rAtt_)
            %Constuctor of class
            obj.ligSourLoc = ligSourLoc_;
            obj.ligSourOrien = ligSourOrien_;
            obj.ligSourDirVec = ligSourOrien_*[0;0;1];
            obj.maxRadiantInt = maxRadiantInt_;
            obj.mu = mu_;
            obj.rAtt = rAtt_;
        end
        
        function PlotLightSource(obj, fig)
            %Plot light source and arrow in passed figure
            figure(fig);
            scatter3(obj.ligSourLoc(1), obj.ligSourLoc(2), obj.ligSourLoc(3), 200, [0,1,0], 'filled');
            arrow3(obj.ligSourLoc', (obj.ligSourLoc + 0.5*obj.ligSourDirVec)', 'v', 5);
            
            axis equal;
            drawnow();
        end
        
        function PlotRadiantSlice(obj, fig, S, sliceDist)
            %plot a slice of the light source along the principal axis into
            %the passed in figure
            
            figure(fig);
            
            %create slice mesh points
            x = linspace(-0.5, 0.5, 100);
            y = linspace(-0.5, 0.5, 100);
            z = sliceDist;
            
            [X, Y, Z] = meshgrid(x,y,z);
            
            %rotate mesh using light source principal direction
            XYZrot = [X(:),Y(:),Z(:)]*obj.ligSourOrien';
            
            rows = size(X);
            
            X = reshape(XYZrot(:,1),rows) + obj.ligSourLoc(1);
            Y = reshape(XYZrot(:,2),rows) + obj.ligSourLoc(2);
            Z = reshape(XYZrot(:,3),rows) + obj.ligSourLoc(3);
            
            radIntMag = RadiantIntensityMesh(obj, X, Y, Z);
            
            S.XData = X;
            S.YData = Y;
            S.ZData = Z;
            S.CData = radIntMag;
            S.Visible = 'on';
            S.EdgeColor = 'none';
            
            
            colormap(jet(100));
            colorbar; caxis([0, obj.maxRadiantInt]);
        end
        
        function [radIntMag, radIntVec] = RadiantIntensityMesh(obj, X, Y, Z)
            %calculate the radiant intensity of the light source for a
            %given mesh of points
            [rows, cols] = size(X);
            
            radIntMag = zeros(rows, cols);
            radIntVec = zeros(rows, cols, 3);
            
            %pre-subtract for creating light vector
            X = X - obj.ligSourLoc(1);
            Y = Y - obj.ligSourLoc(2);
            Z = Z - obj.ligSourLoc(3);
            
            for i = 1:rows
                for j = 1:cols
                    ligVec = [ X(i,j); Y(i,j); Z(i,j)];
                    
                    S_P =  obj.NonIsotropicDiskLightModel(obj.ligSourDirVec, obj.maxRadiantInt, obj.mu, ligVec, obj.rAtt);
                    
                    radIntMag(i,j) = norm(S_P);
                    
                    radIntVec(i,j,:) = S_P;
                end
            end
        end
        
        function [radIntMag, radIntVec] = RadiantIntensityAtPoint(obj, pnt)
            %calculate the radiant intensity at a given point from the
            %source
            
            %negative of light vector
            ligVec = [
                pnt(1) - obj.ligSourLoc(1);
                pnt(2) - obj.ligSourLoc(2);
                pnt(3) - obj.ligSourLoc(3);
                ];
            
            radIntVec =  obj.NonIsotropicDiskLightModel(obj.ligSourDirVec, obj.maxRadiantInt, obj.mu, ligVec, obj.rAtt);
            
            radIntMag = norm(radIntVec);
        end
        
        function RadiIn = RadianceInMaterialPoint(obj, pnt, normal)
            %calculate the in going radiance from the source to some point
            %on the material given its normal and location.
            
            %Radiant intensity at point
            [~, radIntVec] = RadiantIntensityAtPoint(obj, pnt);
            
            %Radiance received at point on surface
            RadiIn = dot(-radIntVec, normal);
            
            %incoming radiance cannot be negative.
            if RadiIn < 0
               RadiIn = 0; 
            end           
        end
        
        function RadiOut = RadianceOutMaterialPoint(obj, pnt, normal, reflectance)
            %calculate the out going radiance that is reflected/emitted
            %at a point on a surface with a given normal and reflectance
            %which receives incoming radiance from the source.
            
            %Radiant intensity at point
            [~, radIntVec] = RadiantIntensityAtPoint(obj, pnt);
            
            %Radiance received at point on surface
            RadiIn = dot(-radIntVec, normal);
            
            if RadiIn < 0
               RadiIn = 0; 
            end
            
            %Radiance outgoing at point on surface
            RadiOut = (reflectance/pi) * RadiIn;
            
        end
        
        function lightDir = get_SourDirVec(obj)
            lightDir = obj.ligSourDirVec;
            return
        end
    end
    
    methods (Static, Access = private)
        function S_P = NonIsotropicDiskLightModel(vDir, Phi0, mu, l, rD)
            %The model for a non-istotropic near field disk light model which
            %includes a cosine RID, and light attenutation modelled by an equation which
            %does not contain a singularity.
            %Inputs:
            %       vDir - principal direction vector of the light source
            %       Phi0 - maximum radiant intensity along the principal direction
            %       mu - parameter which alters of the shape of the RID to suit the
            %       required drop-off with angle theta
            %       l - light source direction from the source to some point in 3D
            %       space
            %       r - effective radius of the disk which alters the shape of the
            %       attenuation drop-off with distance
            %Output:
            %       S_P -  radiant intensity recieved at point P from the light source
            %       as a vector. Norm of this will give you the radiant intensity
            
            %Author: Jasprabhjit Mehami, 13446277
            
            lnorm = norm(l);
            
            l_dot_vDir = dot(l,vDir);
            
            if (l_dot_vDir < 0)
                l_dot_vDir = 0;
            end
            
            %RID
            PhiTheta = Phi0 * (l_dot_vDir/lnorm)^mu;
            
            %attenuation
            f_att = 1/((lnorm/rD + 1)^2);
            
            %normalise light vector (direction vector of l)
            dir_l = l./lnorm;
            
            S_P = PhiTheta .* f_att .* dir_l;
            
        end
    end
end
