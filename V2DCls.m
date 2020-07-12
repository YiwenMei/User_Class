classdef V2DCls
% V2DCls is a class for spatial raster image.

properties
  Fnm % Name of the file(s)
  vtp % Type of variable
  ndv % No data value for the variable
  Ulm % Physical upper limit
  Llm % Physical lower limit
  GIC % Convention of the geographic information ('R'/'I' for regular/irregular grid)
  GIf % Geographic information of the variable ([xl yt;xr yb;Rx Ry] where x/y/R is
      %  the horizontal/vertical/resolution, l/r/b/t stands for left/right/bottom/top)

  Vnm % Name of the variable
  unt % Unit of the variable
end

methods
%% Object building
  function obj=V2DCls(Fnm,vtp,ndv,Ulm,Llm,GIC,GIf,varargin)
    narginchk(7,9);
    ips=inputParser;
    ips.FunctionName=mfilename;

    addRequired(ips,'Fnm',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'Fnm'));
    addRequired(ips,'vtp',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'vtp'));
    addRequired(ips,'ndv',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'ndv'));
    addRequired(ips,'Ulm',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'Ulm'));
    addRequired(ips,'Llm',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'Llm'));
    addRequired(ips,'GIC',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'GIC'));
    if strcmp(GIC,'Bound')
      addRequired(ips,'GIf',@(x) validateattributes(x,{'double'},{'size',[3,2]},mfilename,'GIf'));
    elseif strcmp(GIC,'Grids')
      addRequired(ips,'GIf',@(x) validateattributes(x,{'struct'},{'nonempty'},mfilename,'GIf'));
    else
      error('Spatial extend convention must be "Bound/Grids" for boundary/grids');
    end

    addOptional(ips,'Vnm','',@(x) validateattributes(x,{'char'},{},mfilename,'Vnm'));
    addOptional(ips,'unt','-',@(x) validateattributes(x,{'char'},{},mfilename,'unt'));

    parse(ips,Fnm,vtp,ndv,Ulm,Llm,GIC,GIf,varargin{:});
    Vnm=ips.Results.Vnm;
    unt=ips.Results.unt;
    clear ips varargin

    obj.Fnm=Fnm;
    obj.vtp=vtp;
    obj.ndv=ndv;
    obj.Ulm=Ulm;
    obj.Llm=Llm;
    obj.GIC=GIC;
    obj.GIf=GIf;

    obj.Vnm=Vnm;
    obj.unt=unt;
  end

%% Forcing variable reading
  function v2d=readCls(obj)
    [~,nm,fex]=fileparts(obj.Fnm);
    switch fex
      case {'.tif','tiff'} % compatable for .tiff
        v2d=double(imread(obj.Fnm));
        nm=[nm fex];

      case {'.nc4','nc'} % compatable for .nc
        v2d=double(ncread(obj.Fnm,obj.Vnm))';
        nm=[nm fex ':' obj.Vnm];

      case {'.hdf','hdf5'} % compatable for .hdf5
        v2d=double(hdfread(obj.Fnm,obj.Vnm));
        nm=[nm fex ':' obj.Vnm];

      case {'.asc','.txt'}
        v2d=double(readmatrix(obj.Fnm,'Delimiter',',','NumHeaderLines',5));
        nm=[nm fex];

      case '.mat'
        v2d=matfile(obj.Fnm);
        eval(sprintf('v2d=v2d.%s;',obj.Vnm));
        nm=[nm fex ':' obj.Vnm];
    end
    if ~isnan(obj.ndv)
      v2d(v2d==obj.ndv)=NaN;
    end

% Check the boundary
    k=v2d>obj.Ulm | v2d<obj.Llm;
    N=length(find(k))/length(find(~isnan(v2d)));
    if N~=0
      if N>.05
        warning('%s has more than 5%% of data points out of bound\n',nm);
      else
        fprintf('%d points found for %s\n',length(find(k)),nm);
      end
      [X,Y]=meshgrid(1:size(v2d,2),1:size(v2d,1));
      id=knnsearch([X(~isnan(v2d) & ~k) Y(~isnan(v2d) & ~k)],[X(k) Y(k)],'K',4);
      v=v2d(~isnan(v2d))';
      v2d(k)=mean(v(id),2);
    end
  end

%% Grids of variable
  function [X,Y,sz,rsn]=GridCls(obj)
    switch obj.GIC
      case 'Bound'
        X=obj.GIf(1,1)+obj.GIf(3,1)/2:obj.GIf(3,1):obj.GIf(2,1)-obj.GIf(3,1)/2;
        Y=obj.GIf(1,2)-obj.GIf(3,2)/2:-obj.GIf(3,2):obj.GIf(2,2)+obj.GIf(3,2)/2;
        rsn=obj.GIf(3,:);
      case 'Grids'
        X=obj.GIf.x;
        Y=obj.GIf.y;        
        rsn=round(abs([mean(diff(obj.GIf.x)) mean(diff(obj.GIf.y))]));
    end
    [X,Y]=meshgrid(X,Y);
    sz=size(X);
  end
end
end
