classdef Wind2DTCls
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
  GIf % Geographic information of the variables ([xl yt;xr yb;Rx Ry] where x/y/R
      %  is the horizontal/vertical/resolution, l/r/b/t stands for left/right/bottom/top)
  ofs % offset to UTC in hour
  TmC % Time-window convention
  TmR % Time resolution in number of hours
  TmF % Format of the file names

  Vnm % Cell array stores the names of the variables (must have the same dimension
      %  as Fnm within each cell for the corresponding variables except the case
      %  of geotiff file)
  unt % Unit of the variable (if it is not specified, m/s is used)
end

methods
%% Object building
  function obj=Wind2DTCls(Fnm,ndv,Ulm,Llm,GIf,ofs,TmC,TmR,TmF,varargin)
    narginchk(9,11);
    ips=inputParser;
    ips.FunctionName=mfilename;

    addRequired(ips,'Fnm',@(x) validateattributes(x,{'cell'},{'nonempty'},mfilename,'Fnm'));
    addRequired(ips,'ndv',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'ndv'));
    addRequired(ips,'Ulm',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'Ulm'));
    addRequired(ips,'Llm',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'Llm'));
    addRequired(ips,'GIf',@(x) validateattributes(x,{'double'},{'size',[3,2]},mfilename,'GIf'));
    addRequired(ips,'ofs',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'ofs'));
    addRequired(ips,'TmC',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'TmC'));
    addRequired(ips,'TmR',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'TmR'));
    addRequired(ips,'TmF',@(x) validateattributes(x,{'cell'},{'numel',2},mfilename,'TmF'));

    addOptional(ips,'Vnm',{},@(x) validateattributes(x,{'cell'},{},mfilename,'Vnm'));
    addOptional(ips,'unt','m/s',@(x) validateattributes(x,{'char'},{},mfilename,'unt'));

    parse(ips,Fnm,ndv,Ulm,Llm,GIf,ofs,TmC,TmR,TmF,varargin{:});
    Vnm=ips.Results.Vnm;
    unt=ips.Results.unt;
    clear ips varargin

    obj.Fnm=Fnm;
    obj.vtp='Wspd';
    obj.ndv=ndv;
    obj.Ulm=Ulm;
    obj.Llm=Llm;
    obj.GIf=GIf;
    obj.ofs=ofs;
    obj.TmC=TmC;
    obj.TmR=TmR;
    obj.TmF=TmF;

    obj.Vnm=Vnm;
    obj.unt=unt;
  end

%% Forcing variable reading
  function [ws,wa,U,V]=readCls(obj,n)
    U=[];
    for v=1:size(obj.Fnm{n},1)
      [~,~,fex]=fileparts(obj.Fnm{n}(v,:));
      switch fex
        case {'.tif','tiff'} % compatable for .tiff
          U=cat(3,U,double(imread(obj.Fnm{n}(v,:))));

        case {'.nc4','nc'} % compatable for .nc
          U=cat(3,U,double(ncread(obj.Fnm{n}(v,:),obj.Vnm{v}))');

        case {'.hdf','hdf5'} % compatable for .hdf5
          U=cat(3,U,double(hdfread(obj.Fnm{n}(v,:),obj.Vnm{v})));

        case {'.asc','.txt'}
          U=cat(3,U,double(readmatrix(obj.Fnm{n}(v,:),'Delimiter',obj.Vnm,'NumHeaderLines',5)));

        case '.mat'
          V=matfile(obj.Fnm{n}(v,:));
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
    X=obj.GIf(1,1)+obj.GIf(3,1)/2:obj.GIf(3,1):obj.GIf(2,1)-obj.GIf(3,1)/2;
    Y=obj.GIf(1,2)-obj.GIf(3,2)/2:-obj.GIf(3,2):obj.GIf(2,2)+obj.GIf(3,2)/2;
    [X,Y]=meshgrid(X,Y);
  end

%% Extract the time line
  function Tutc=TimeCls(obj,Cflg)
% Convert to the UTC time line
    [~,ds,~]=cellfun(@(X) fileparts(X(1,:)),obj.Fnm,'UniformOutput',false);
    ds=cellfun(@(X) X(length(obj.TmF{1})+1:end),ds,'UniformOutput',false);
    Tutc=datenum(cell2mat(ds),obj.TmF{2});
    if obj.ofs~=0
      Tutc=Tutc+obj.ofs;
    end

% Adjust to user specified time convention 
    switch Cflg
      case 'begin' % To
        switch obj.TmC
          case 'center' % From
            Tutc=Tutc-obj.TmR/24/2;
          case 'end' % From
            Tutc=Tutc-obj.TmR/24;
        end

      case 'center'
        switch obj.TmC
          case 'begin'
            Tutc=Tutc+obj.TmR/24/2;
          case 'end'
            Tutc=Tutc-obj.TmR/24/2;
        end

      case 'end'
        switch obj.TmC
          case 'center'
            Tutc=Tutc+obj.TmR/24/2;
          case 'begin'
            Tutc=Tutc+obj.TmR/24;
        end
    end
  end
end
end
