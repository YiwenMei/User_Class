classdef Wind2DCls
% Wind2DTCls is a class for a stack of spatial raster image of wind with a time line.

properties
  Fnm % Cell array stores the names of the input files (Each cell has one or
      %  two rows for, if Typ is 'Total', the name of the total wind file or
      %  both the names of the total wind and wind direction file; if Typ is
      %  'Component', the two rows store file names for the U and V wind)
  vtp % Type of variable
  ndv % No data value for the variables
  Ulm % Physical upper limit of the total wind
  Llm % Physical lower limit of the total wind
  Gtg % Geographic information of the variables ([xl yt;xr yb;Rx Ry] where x/y/R
      %  is the horizontal/vertical/resolution, l/r/b/t stands for left/right/bottom/top)

  Vnm % Cell array stores the names of the variables (must have the same dimension
      %  as Fnm within each cell for the corresponding variables except the case
      %  of geotiff file)
  unt % Unit of the variable (if it is not specified, m/s is used)
end

methods
%% Object building
  function obj=Wind2DCls(Fnm,ndv,Ulm,Llm,Gtg,varargin)
    narginchk(5,7);
    ips=inputParser;
    ips.FunctionName=mfilename;

    addRequired(ips,'Fnm',@(x) validateattributes(x,{'char'},{'nrows',2},mfilename,'Fnm'));
    addRequired(ips,'ndv',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'ndv'));
    addRequired(ips,'Ulm',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'Ulm'));
    addRequired(ips,'Llm',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'Llm'));
    addRequired(ips,'Gtg',@(x) validateattributes(x,{'double'},{'size',[3,2]},mfilename,'Gtg'));

    addOptional(ips,'Vnm',{'Uspd','Vspd'},@(x) validateattributes(x,{'cell'},...
        {'numel',2},mfilename,'Vnm'));
    addOptional(ips,'unt','m/s',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'unt'));

    parse(ips,Fnm,ndv,Ulm,Llm,Gtg,varargin{:});
    Vnm=ips.Results.Vnm;
    unt=ips.Results.unt;
    clear ips varargin

    obj.Fnm=Fnm;
    obj.vtp='Wspd';
    obj.ndv=ndv;
    obj.Ulm=Ulm;
    obj.Llm=Llm;
    obj.Gtg=Gtg;

    obj.Vnm=Vnm;
    obj.unt=unt;
  end

%% Forcing variable reading
  function [ws,wa,U,V]=readCls(obj)
    U=[];
    for v=1:size(obj.Fnm,1)
      [~,~,fex]=fileparts(obj.Fnm(v,:));
      switch fex
        case {'.tif','tiff'} % compatable for .tiff
          U=cat(3,U,double(imread(obj.Fnm(v,:))));

        case {'.nc4','nc'} % compatable for .nc
          U=cat(3,U,double(ncread(obj.Fnm(v,:),obj.Vnm{v}))');

        case {'.hdf','hdf5'} % compatable for .hdf5
          U=cat(3,U,double(hdfread(obj.Fnm(v,:),obj.Vnm{v})));

        case {'.asc','.txt'}
          U=cat(3,U,double(readmatrix(obj.Fnm(v,:),'Delimiter',',','NumHeaderLines',5)));

        case '.mat'
          V=matfile(obj.Fnm(v,:));
          eval(sprintf('V=V.%s;',obj.Vnm{v}));
          U=cat(3,U,V);
      end
    end
    U(U==obj.ndv)=NaN;
    V=U(:,:,2);
    U=U(:,:,1);

% Calculat the total wind
    wa=atan2d(V,U); % E is 0, counter-clock's wise is +
    wa(wa<0)=wa(wa<0)+360;
    ws=hypot(U,V);

% Check the output
    validateattributes(ws(~isnan(ws)),{'double'},{'<=',obj.Ulm,'>=',obj.Llm},'','Total wind');
  end

%% Grids of variable
  function [X,Y]=GridCls(obj)
    X=obj.Gtg(1,1)+obj.Gtg(3,1)/2:obj.Gtg(3,1):obj.Gtg(2,1)-obj.Gtg(3,1)/2;
    Y=obj.Gtg(1,2)-obj.Gtg(3,2)/2:-obj.Gtg(3,2):obj.Gtg(2,2)+obj.Gtg(3,2)/2;
    [X,Y]=meshgrid(X,Y);
  end

%% Extract the time line
  function Ttg=TimeCls(obj)
    [~,ds,~]=cellfun(@(X) fileparts(X(1,:)),obj.Fnm,'UniformOutput',false);
    ds=cellfun(@(X) X(length(obj.TmF{1})+1:end),ds,'UniformOutput',false);
    Ttg=datenum(cell2mat(ds),obj.TmF{2});
  end
end
end
