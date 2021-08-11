classdef Hypercube_Visualiser < handle
    % Visualiser hypercubes onto a provided figure with a slider to flick
    % through bands. Each hypercube will be plotted onto a separate subplot
    %
    
    properties(Access = private)
        numCubes_
        fig_
        numBands_
        hCubeCell_
        bandDisp_
        plotTitle_
        imgCell_
    end
    
    methods
        function obj = Hypercube_Visualiser(fig)
            obj.fig_ = fig;
            obj.numCubes_ = 0;
            obj.bandDisp_ = 1;
            obj.hCubeCell_ = {};
            obj.plotTitle_ = "";
            obj.imgCell_ = {};
        end
        
        function AddHypercube(obj, hCube, title)
            % Adds a hypercube to be visualised. Hypercube should have
            % dimensions [pixels, lines, bands]
            
            [~,~,bands] = size(hCube);
   
            %first hypercube
            if obj.numCubes_ < 1
                %all following hypercubes must have this many bands
                obj.numBands_ = bands;
            else
                if bands ~= obj.numBands_
                    warning(['Input hypercube should have ', num2str(obj.numBands_), ' number of bands']);
                    return;
                end
            end
            
            obj.numCubes_ = obj.numCubes_ + 1;
            obj.hCubeCell_(obj.numCubes_) = {mat2gray(hCube)};
            obj.plotTitle_(obj.numCubes_) = title;
            obj.PlotHypercubes();
        end
        
        function AddImage(obj, img, title)
            
        end
        
        function PlotHypercubes(obj)
            % Plot hypercube onto figure with slider
            figure(obj.fig_);
            clf(obj.fig_);
            
            %band slider
            slider_band = uicontrol('Parent',obj.fig_,'Style','slider','Position',[81,120,300,23],...
                'value',1, 'min', 1, 'max',obj.numBands_, 'SliderStep', [1/(obj.numBands_+1), 1/(obj.numBands_+1)]);
            bgcolor = obj.fig_.Color;
            uicontrol('Parent',obj.fig_,'Style','text','Position',[50,90,23,23],...
                'String','1','BackgroundColor',bgcolor);
            uicontrol('Parent',obj.fig_,'Style','text','Position',[350,90,40,23],...
                'String',num2str(obj.numBands_),'BackgroundColor',bgcolor);
            uicontrol('Parent',obj.fig_,'Style','text','Position',[150,90,200,23],...
                'String','Band','BackgroundColor',bgcolor);
            %             disp_dist = uicontrol('Parent',obj.fig,'Style','text','Position',[200,145,50,20],...
            %                 'String', num2str(obj.distFromSrc),'BackgroundColor', [1,1,1]);
            %callback function at the end of the script
            slider_band.Callback = @(src, eventData) obj.band_callback(src, eventData);
            
            obj.PlotHypercubeLoop();
        end
        
        function band_callback(obj, src, ~)
            band = round(src.Value);
            obj.bandDisp_ = band;
            obj.PlotHypercubeLoop();
        end
    end
    
    methods(Access = private)
        function PlotHypercubeLoop(obj)
            n = obj.numCubes_;
            
            for i = 1:n
                subplot(1, n, i, 'replace');
                curCube = obj.hCubeCell_{i};
                img = curCube(:,:,obj.bandDisp_);
                imshow(img);
                title(obj.plotTitle_(i))
            end
        end
    end
end

