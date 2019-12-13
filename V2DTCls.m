classdef V2DTCls
% V2DTCls is a class for a stack of spatial raster image with a time line.

properties
  Fnm % Name of the file(s)
  vtp % Type of variable
  ndv % No data value for the variable
  Ulm % Physical upper limit
  Llm % Physical lower limit
  GIf % Geographic information of the variable ([xl yt;xr yb;Rx Ry] where x/y/R is
      %  the horizontal/vertical/resolution, l/r/b/t stands for left/right/bottom/top)
  ofs % offset to UTC in hour
  TmC % Time-window convention
  TmR % Time resolution in number of hours
  TmF % Format of the file names

  Vnm % Name of the variable
  unt % Unit of the variable
end

methods
%% Object building
  function obj=V2DTCls(Fnm,vtp,ndv,Ulm,Llm,GIf,ofs,TmC,TmR,TmF,varargin)
    narginchk(10,12);
    ips=inputParser;
    ips.FunctionName=mfilename;

    addRequired(ips,'Fnm',@(x) validateattributes(x,{'cell'},{'nonempty'},mfilename,'Fnm'));
    addRequired(ips,'vtp',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'vtp'));
    addRequired(ips,'ndv',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'ndv'));
    addRequired(ips,'Ulm',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'Ulm'));
    addRequired(ips,'Llm',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'Llm'));
    addRequired(ips,'GIf',@(x) validateattributes(x,{'double'},{'size',[3,2]},mfilename,'GIf'));
    addRequired(ips,'ofs',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'ofs'));
    addRequired(ips,'TmC',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'TmC'));
    addRequired(ips,'TmR',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'TmR'));
    addRequired(ips,'TmF',@(x) validateattributes(x,{'cell'},{'numel',2},mfilename,'TmF'));

    addOptional(ips,'Vnm','',@(x) validateattributes(x,{'char'},{},mfilename,'Vnm'));
    addOptional(ips,'unt','-',@(x) validateattributes(x,{'char'},{},mfilename,'unt'));

    parse(ips,Fnm,vtp,ndv,Ulm,Llm,GIf,ofs,TmC,TmR,TmF,varargin{:});
    Vnm=ips.Results.Vnm;
    unt=ips.Results.unt;
    clear ips varargin

    obj.Fnm=Fnm;
    obj.vtp=vtp;
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
  function v2d=readCls(obj,n)
    [~,nm,fex]=fileparts(obj.Fnm{n});
    switch fex
      case {'.tif','tiff'} % compatable for .tiff
        v2d=double(imread(obj.Fnm{n}));
        nm=[nm fex];

      case {'.nc4','nc'} % compatable for .nc
        v2d=double(ncread(obj.Fnm{n},obj.Vnm))';
        nm=[nm fex ':' obj.Vnm];

      case {'.hdf','hdf5'} % compatable for .hdf5
        v2d=double(hdfread(obj.Fnm{n},obj.Vnm));
        nm=[nm fex ':' obj.Vnm];

      case {'.asc','.txt'}
        v2d=double(readmatrix(obj.Fnm{n},'Delimiter',obj.Vnm,'NumHeaderLines',5));
        nm=[nm fex];

      case '.mat'
        v2d=matfile(obj.Fnm{n});
        eval(sprintf('v2d=v2d.%s;',obj.Vnm));
        nm=[nm fex ':' obj.Vnm];
    end
    v2d(v2d==obj.ndv)=NaN;

% Check the boundary
    validateattributes(v2d(~isnan(v2d)),{'double'},{'<=',obj.Ulm,'>=',obj.Llm},'',nm);
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
    [~,ds,~]=cellfun(@(X) fileparts(X),obj.Fnm,'UniformOutput',false);
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
