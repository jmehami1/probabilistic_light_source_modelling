classdef LightSimulator < handle
    %Class used to simulate a near-field non-isotropic disk light source
    %which is placed an an environment relative to a frame camera origin
    
    properties(Access = private)
        %Pose of source w.r.t to frame camera
        ligSourLoc
        ligSourOrien
        ligSourDirVec
        
        maxRadiantInt
        mu
        rAtt
        distFromSrc
        fig
        surfDir
        targetPts
        targetbandInten
        numBands
        numPts
        targetScat
        RIDspline
        realRID
    end
    
    methods
        function obj = LightSimulator(varargin)
            
            if nargin == 9
                %Constuctor of class
                obj.ligSourLoc = varargin{1};
                obj.ligSourOrien = varargin{2};
                obj.ligSourDirVec = varargin{2}*[0;0;1];
                obj.maxRadiantInt = varargin{3};
                obj.mu = varargin{4};
                obj.rAtt = varargin{5};
                obj.distFromSrc = varargin{6};
                obj.fig = varargin{7};
                obj.surfDir = varargin{8};
                obj.RIDspline = varargin{9};
                
                %plot light source and its directional planar slice (slice will
                %not be visible on start up)
                obj.SetupSimulatorFigure();
                obj.PlotLightSource();
%                 obj.PlotRadiantSlice();

                obj.realRID = true;
            
            elseif nargin == 8
                %Constuctor of class
                obj.ligSourLoc = varargin{1};
                obj.ligSourOrien = varargin{2};
                obj.ligSourDirVec = varargin{2}*[0;0;1];
                obj.maxRadiantInt = varargin{3};
                obj.mu = varargin{4};
                obj.rAtt = varargin{5};
                obj.distFromSrc = varargin{6};
                obj.fig = varargin{7};
                obj.surfDir = varargin{8};
                
                %plot light source and its directional planar slice (slice will
                %not be visible on start up)
                obj.SetupSimulatorFigure();
                obj.PlotLightSource();
%                 obj.PlotRadiantSlice();
            elseif nargin == 6
                %Constuctor of class
                obj.ligSourLoc = varargin{1};
                obj.ligSourOrien = varargin{2};
                obj.ligSourDirVec = varargin{2}*[0;0;1];
                obj.maxRadiantInt = varargin{3};
                obj.mu = varargin{4};
                obj.rAtt = varargin{5};
                obj.fig = varargin{6};
                
                obj.PlotLightSource();

            elseif nargin == 5
                %Constuctor of class
                obj.ligSourLoc = varargin{1};
                obj.ligSourOrien = varargin{2};
                obj.ligSourDirVec = varargin{2}*[0;0;1];
                obj.maxRadiantInt = varargin{3};
                obj.mu = varargin{4};
                obj.rAtt = varargin{5};
                
            elseif nargin == 7
                %Constuctor of class
                obj.ligSourLoc = varargin{1};
                obj.ligSourOrien = varargin{2};
                obj.ligSourDirVec = varargin{2}*[0;0;1];
                obj.fig = varargin{3};
                obj.targetPts = varargin{4};
                obj.targetbandInten = varargin{5};
                obj.numBands = varargin{6};
                obj.numPts = varargin{7};
                
                obj.PlotLightSource();
                obj.PlotTargetPts();
            else
                error("Wrong number of inputs into class");
            end
        end
        
        function PlotTargetPts(obj)
            figure(obj.fig);
            
            X = obj.targetPts(:,1);
            Y = obj.targetPts(:,2);
            Z = obj.targetPts(:,3);
            intenFirstBand = obj.targetbandInten(:,1);
            
            obj.targetScat = scatter3(X, Y, Z, 5, intenFirstBand, 'filled'); hold on;
            caxis([0, max(intenFirstBand, [], 'all')])        
            colormap(jet(1000)); colorbar;
            
            %band slider
            slider_band = uicontrol('Parent',obj.fig,'Style','slider','Position',[81,120,300,23],...
                'value',1, 'min', 1, 'max',obj.numBands, 'SliderStep', [1/(obj.numBands+1), 1/(obj.numBands+1)]);
            bgcolor = obj.fig.Color;
            uicontrol('Parent',obj.fig,'Style','text','Position',[50,90,23,23],...
                'String','1','BackgroundColor',bgcolor);
            uicontrol('Parent',obj.fig,'Style','text','Position',[350,90,40,23],...
                'String',num2str(obj.numBands),'BackgroundColor',bgcolor);
            uicontrol('Parent',obj.fig,'Style','text','Position',[150,90,200,23],...
                'String','Band','BackgroundColor',bgcolor);
%             disp_dist = uicontrol('Parent',obj.fig,'Style','text','Position',[200,145,50,20],...
%                 'String', num2str(obj.distFromSrc),'BackgroundColor', [1,1,1]);
            %callback function at the end of the script
            slider_band.Callback = @(src, eventData) obj.band_callback(src, eventData);
            
            drawnow();
        end
            
        
        
        %         function obj = LightSimulator(ligSourLoc_,ligSourOrien_, maxRadiantInt_, mu_, rAtt_, distFromSrc_, fig_, surfDir_)
        %             %Constuctor of class
        %             obj.ligSourLoc = ligSourLoc_;
        %             obj.ligSourOrien = ligSourOrien_;
        %             obj.ligSourDirVec = ligSourOrien_*[0;0;1];
        %             obj.maxRadiantInt = maxRadiantInt_;
        %             obj.mu = mu_;
        %             obj.rAtt = rAtt_;
        %             obj.distFromSrc = distFromSrc_;
        %             obj.fig = fig_;
        %             obj.surfDir = surfDir_;
        %
        %             %plot light source and its directional planar slice (slice will
        %             %not be visible on start up)
        %             obj.PlotLightSource();
        %             obj.PlotRadiantSlice();
        %         end
        
        function PlotLightSource(obj)
            %Plot light source and arrow in passed figure
            figure(obj.fig);
            scatter3(obj.ligSourLoc(1), obj.ligSourLoc(2), obj.ligSourLoc(3), 250, [0,1,0], 'filled');
            arrow3(obj.ligSourLoc', (obj.ligSourLoc + 0.5*obj.ligSourDirVec)', 'v', 5);
            
            poseLightSrc = [obj.ligSourOrien, obj.ligSourLoc; 0,0,0,1];
            trplot(poseLightSrc, 'rviz', 'length', 0.1);
            
            axis equal; drawnow();
        end
        
        function PlotRadiantSlice(obj)
            %plot a slice of the light source along the principal axis into
            %the passed in figure
            
            figure(obj.fig);
            
            %create slice mesh points
            x = linspace(-0.5, 0.5, 100);
            y = linspace(-0.5, 0.5, 100);
            z = obj.distFromSrc;
            
            [X, Y, Z] = meshgrid(x,y,z);
            
            %rotate mesh using light source principal direction
            XYZrot = [X(:),Y(:),Z(:)]*obj.ligSourOrien';
            
            rows = size(X);
            
            X = reshape(XYZrot(:,1),rows) + obj.ligSourLoc(1);
            Y = reshape(XYZrot(:,2),rows) + obj.ligSourLoc(2);
            Z = reshape(XYZrot(:,3),rows) + obj.ligSourLoc(3);
            
            radIntMag = RadiantIntensityMesh(obj, X, Y, Z);
            
            obj.surfDir.XData = X;
            obj.surfDir.YData = Y;
            obj.surfDir.ZData = Z;
            obj.surfDir.CData = radIntMag;
            obj.surfDir.EdgeColor = 'none';
            
            colormap(hot);
             colorbar; caxis([0, obj.maxRadiantInt]);
            
            axis equal; drawnow();
        end
        
        function PlotHemisphereMesh(obj, noPtsHem, orien, locCentre, radius, S)
            %INPUTS:
            %       noPtsHem - number of points used for hemisphere mesh
            %       orien - rotation of hemisphere relative to frame camera
            %       locCentre -  centre of hemisphere relative to frame camera
            %       radius - radius of hemisphere
            %       S - surface mesh object
            
            %create unit hemisphere mesh
            [X,Y,Z] = sphere(noPtsHem);
            [rows, ~] = size(Z);
            
            %scale hemisphere using radius
            X = radius.*X((rows-1)/2:end,:);
            Y = radius.*Y((rows-1)/2:end,:);
            Z = radius.*Z((rows-1)/2:end,:);
            
            %rotate
            XYZrot = [X(:),Y(:),Z(:)]*orien;
            
            rows_cols = size(X);
            
            %translate
            X = reshape(XYZrot(:,1),rows_cols) + locCentre(1);
            Y = reshape(XYZrot(:,2),rows_cols) + locCentre(2);
            Z = reshape(XYZrot(:,3),rows_cols) + locCentre(3);
            
            radianceInSphere = zeros(size(X));
            
            for i = 1:rows_cols(1)
                for j = 1:rows_cols(2)
                    pntSph = [X(i,j); Y(i,j); Z(i,j)];
                    
                    %calculate normal
                    normSph = pntSph - locCentre';
                    %normalise vector to get direction vector
                    normSph = normSph/norm(normSph);
                    
                    radianceInSphere(i,j) = RadianceInMaterialPoint(obj, pntSph, normSph);
                end
            end
            
            %             radIntMag = RadiantIntensityMesh(obj, X, Y, Z);
            
            figure(obj.fig);
            
            S.XData = X;
            S.YData = Y;
            S.ZData = Z;
            S.CData = radianceInSphere;
            %             obj.surfDir.Visible = 'on';
            S.EdgeColor = 'none';
            
            
            axis equal; drawnow();
        end
        
        function [radIntMag, radIntVec] = RadiantIntensityMesh(obj, X, Y, Z)
            %calculate the radiant intensity of the light source for a
            %given mesh of points
            [rows, cols] = size(X);
            
            radIntMag = zeros(rows, cols);
            radIntVec = zeros(rows, cols, 3);
            
            %pre-subtract for creating light vector
            X = obj.ligSourLoc(1) - X;
            Y = obj.ligSourLoc(2) - Y;
            Z = obj.ligSourLoc(3) - Z;
            
            for i = 1:rows
                for j = 1:cols
                    ligVec = [ X(i,j); Y(i,j); Z(i,j)];
                    
                    
                    if obj.realRID
                        S_P =  obj.NonIsotropicDiskLightModelRID(obj.ligSourDirVec, obj.maxRadiantInt, obj.mu, ligVec, obj.rAtt, obj.RIDspline);
                    else
                        S_P =  obj.NonIsotropicDiskLightModel(obj.ligSourDirVec, obj.maxRadiantInt, obj.mu, ligVec, obj.rAtt);
                    end
                    
%                     S_P =  obj.NonIsotropicDiskLightModel(obj.ligSourDirVec, obj.maxRadiantInt, obj.mu, ligVec, obj.rAtt);
                    
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
                obj.ligSourLoc(1) - pnt(1);
                obj.ligSourLoc(2) - pnt(2);
                obj.ligSourLoc(3) - pnt(3);
                ];
            
            if obj.realRID
                radIntVec =  obj.NonIsotropicDiskLightModelRID(obj.ligSourDirVec, obj.maxRadiantInt, obj.mu, ligVec, obj.rAtt, obj.RIDspline);
            else
                radIntVec =  obj.NonIsotropicDiskLightModel(obj.ligSourDirVec, obj.maxRadiantInt, obj.mu, ligVec, obj.rAtt);
            end
            
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
        
        function loc = get_ligSourLoc(obj)
            loc = obj.ligSourLoc;
            return;
        end
        
        function band_callback(obj, src, ~)
            i = round(src.Value);
            curBandInt = obj.targetbandInten(:,i);
            obj.targetScat.CData = curBandInt;
            caxis([0, max(curBandInt, [], 'all')])
        end
        
        
        function distFromSrc_callback(obj, src, ~, texboxHandle)
            obj.distFromSrc = src.Value;
            obj.PlotRadiantSlice();
            texboxHandle.String = num2str(src.Value);
        end
        
        function mu_callback(obj, src, ~, texboxHandle)
            obj.mu = src.Value;
            obj.PlotRadiantSlice();
            texboxHandle.String = num2str(src.Value);
        end
        
        function surfDir_visible_callback(obj, src, ~)
            value = get(src,'Value');
            
            if value
                obj.surfDir.Visible = 'on';
            else
                obj.surfDir.Visible = 'off';
            end
        end
        
        
        function SetupSimulatorFigure(obj)
            figure(obj.fig);
            %distFromSrc slider
            slider_dist = uicontrol('Parent',obj.fig,'Style','slider','Position',[81,120,300,23],...
                'value',obj.distFromSrc, 'min', 0, 'max',1, 'SliderStep', [0.1/1, 0.1/1]);
            bgcolor = obj.fig.Color;
            uicontrol('Parent',obj.fig,'Style','text','Position',[50,90,23,23],...
                'String','0','BackgroundColor',bgcolor);
            uicontrol('Parent',obj.fig,'Style','text','Position',[350,90,23,23],...
                'String','1','BackgroundColor',bgcolor);
            uicontrol('Parent',obj.fig,'Style','text','Position',[150,90,200,23],...
                'String','Distance from the source','BackgroundColor',bgcolor);
            disp_dist = uicontrol('Parent',obj.fig,'Style','text','Position',[200,145,50,20],...
                'String', num2str(obj.distFromSrc),'BackgroundColor', [1,1,1]);
            %callback function at the end of the script
            slider_dist.Callback = @(src, eventData) obj.distFromSrc_callback(src, eventData, disp_dist);
            
            %mu slider
            slider_mu = uicontrol('Parent',obj.fig,'Style','slider','Position',[81,200,300,23],...
                'value', obj.mu, 'min',0, 'max',10, 'SliderStep', [0.5/10, 0.5/10]);
            bgcolor = obj.fig.Color;
            uicontrol('Parent',obj.fig,'Style','text','Position',[50,170,23,23],...
                'String','0','BackgroundColor',bgcolor);
            uicontrol('Parent',obj.fig,'Style','text','Position',[350,170,23,23],...
                'String','10','BackgroundColor',bgcolor);
            uicontrol('Parent',obj.fig,'Style','text','Position',[150,170,200,23],...
                'String','mu','BackgroundColor',bgcolor);
            disp_mu = uicontrol('Parent',obj.fig,'Style','text','Position',[200,225,50,20],...
                'String', num2str(obj.mu),'BackgroundColor', [1,1,1]);
            %callback function at the end of the script
            slider_mu.Callback = @(src, eventData) obj.mu_callback(src, eventData, disp_mu);
            
            
            %checkerbox to turn surface plane visibility on/off
            surfDirDispCheck = uicontrol('Parent',obj.fig,'Style','checkbox', 'String', 'Display Direction Plane', 'Position', [20,45,200,20] );
            surfDirDispCheck.Callback = @(src, eventData) obj.surfDir_visible_callback(src, eventData);
            
        end
        
        
    end
    
    
    
    
    
    
    methods (Static, Access = private)
         function S_P = NonIsotropicDiskLightModelRID(vDir, Phi0, mu, l, rD, ridSpline)
            %The model for a non-istotropic near field disk light model which
            %includes a cosine RID, and light attenutation modelled by an equation which
            %does not contain a singularity.
            %Inputs:
            %       vDir - principal direction vector of the light source
            %       Phi0 - maximum radiant intensity along the principal direction
            %       mu - parameter which alters of the shape of the RID to suit the
            %       required drop-off with angle theta
            %       l - light vector from the point-of-interest to the
            %       source
            %       rd - effective radius of the disk which alters the shape of the
            %       attenuation drop-off with distance
            %Output:
            %       S_P -  radiant intensity from the source to the
            %       point-of-interest along the light vector.
            
            %Author: Jasprabhjit Mehami, 13446277
            
            lnorm = norm(l);
            
            %transform light vector to light source c.f.
            
            
            
            %normalise light vector (direction vector of l)
             dir_l = -l./lnorm;
            
                        
            l_dot_vDir = dot(-l,vDir);
            
            if (l_dot_vDir < 0)
                l_dot_vDir = 0;
            end
            
            theta_l = acos(l_dot_vDir ./ lnorm);

            if dir_l(1) < vDir(1)
                theta_l = -theta_l;
            end
            
            
            %RID
            PhiTheta = Phi0 * fnval(ridSpline, theta_l);
            
            %attenuation
            f_att = 1/((lnorm/rD + 1)^2);
            
            %normalise light vector (direction vector of l)
            dir_l = -l./lnorm;
            
            %radiant intensity vector from the source to the
            %point-of-interest. For plotting we just use the vector, but
            %when the ray interacts with a surface, the negative of the
            %vector must be taken to ensure the dot-product between them is
            %positive.
            S_P = PhiTheta .* f_att .* dir_l;
            
        end
        
        
        function S_P = NonIsotropicDiskLightModel(vDir, Phi0, mu, l, rD)
            %The model for a non-istotropic near field disk light model which
            %includes a cosine RID, and light attenutation modelled by an equation which
            %does not contain a singularity.
            %Inputs:
            %       vDir - principal direction vector of the light source
            %       Phi0 - maximum radiant intensity along the principal direction
            %       mu - parameter which alters of the shape of the RID to suit the
            %       required drop-off with angle theta
            %       l - light vector from the point-of-interest to the
            %       source
            %       rd - effective radius of the disk which alters the shape of the
            %       attenuation drop-off with distance
            %Output:
            %       S_P -  radiant intensity from the source to the
            %       point-of-interest along the light vector.
            
            %Author: Jasprabhjit Mehami, 13446277
            
            lnorm = norm(l);
            
            l_dot_vDir = dot(-l,vDir);
            
            if (l_dot_vDir < 0)
                l_dot_vDir = 0;
            end
            
            %RID
            PhiTheta = Phi0 * (l_dot_vDir/lnorm)^mu;
            
            %attenuation
            f_att = 1/((lnorm/rD + 1)^2);
            
            %normalise light vector (direction vector of l)
            dir_l = -l./lnorm;
            
            %radiant intensity vector from the source to the
            %point-of-interest. For plotting we just use the vector, but
            %when the ray interacts with a surface, the negative of the
            %vector must be taken to ensure the dot-product between them is
            %positive.
            S_P = PhiTheta .* f_att .* dir_l;
            
        end
        
    end
end
