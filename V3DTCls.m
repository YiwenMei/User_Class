classdef V3DTCls
% V3DTCls is a class for a series of 3D image with a time line.

properties
  Fnm % Name of the file(s)
  vtp % Type of variable
  ndv % No data value for the variable
  Ulm % Physical upper limit
  Llm % Physical lower limit
  GIC % Convention of the geographic information (can be supplied as a retagualr domain "Boundary" or
      %  "Coordinate" of horizontal and vertical grid center)
  GIf % Geographic information of the variable; for 'Boundary', it is specified as a 3-by-2 matrix
      %  [xl yt;xr yb;Rx Ry] where x/y/R is the horizontal/vertical/resolution, l/r/b/t stands for
      %  left/right/bottom/top; for 'Coordinate', it is a structure array with two field, x & y,
      %  for the horizontal and vertical coordinate of grid center
  VHt % Height(-)/depth(+) of vertical layer
  ofs % offset to UTC in hour
  TmC % Time-window convention
  TmR % Time resolution in number of hours
  TmF % Format of the file names

  Vnm % Name of the variable
  unt % Unit of the variable
end

methods
%% Object building
  function obj=V3DTCls(Fnm,vtp,ndv,Ulm,Llm,GIC,GIf,VHt,ofs,TmC,TmR,TmF,varargin)
    narginchk(12,14);
    ips=inputParser;
    ips.FunctionName=mfilename;

    addRequired(ips,'Fnm',@(x) validateattributes(x,{'cell'},{'nonempty'},mfilename,'Fnm'));
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
    addRequired(ips,'VHt',@(x) validateattributes(x,{'double'},{'vector'},mfilename,'VHt'));
    addRequired(ips,'ofs',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'ofs'));
    addRequired(ips,'TmC',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'TmC'));
    addRequired(ips,'TmR',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'TmR'));
    addRequired(ips,'TmF',@(x) validateattributes(x,{'cell'},{'numel',4},mfilename,'TmF'));

    addOptional(ips,'Vnm','',@(x) validateattributes(x,{'char'},{},mfilename,'Vnm'));
    addOptional(ips,'unt','-',@(x) validateattributes(x,{'char'},{},mfilename,'unt'));

    parse(ips,Fnm,vtp,ndv,Ulm,Llm,GIC,GIf,VHt,ofs,TmC,TmR,TmF,varargin{:});
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
    obj.VHt=VHt;
    obj.ofs=ofs;
    obj.TmC=TmC;
    obj.TmR=TmR;
    obj.TmF=TmF;

    obj.Vnm=Vnm;
    obj.unt=unt;
  end

%% Read one/all layer of a time step
  function v2d=readCls(obj,n,varargin)
    l=cell2mat(varargin);
    [~,nm,fex]=fileparts(obj.Fnm{n});
    nm=[nm fex];

    switch obj.TmF{3}
      case {'nc4','nc'} % compatable for .nc
        if ~isempty(l)
          L=ones(1,4);
          L(obj.TmF{4}(3))=l;
          v2d=ncread(obj.Fnm{n},obj.Vnm,L,[Inf 1 Inf 1]);
          v2d=double(flipud(permute(v2d,obj.TmF{4})));
          nm=sprintf('%s-%s-L%i',nm,obj.Vnm,l);
        else
          v2d=ncread(obj.Fnm{n},obj.Vnm,[1 1 1 1],[Inf Inf Inf 1]);
          v2d=double(flipud(permute(v2d,obj.TmF{4})));
          nm=sprintf('%s-%s-All',nm,obj.Vnm);
        end

      case 'mat'
        v2d=matfile(obj.Fnm{n});
        eval(sprintf('v2d=v2d.%s;',obj.Vnm));
        if ~isempty(l)
          v2d=v2d(:,:,l);
          nm=sprintf('%s-%s-L%i',nm,obj.Vnm,l);
        else
          nm=sprintf('%s-%s-All',nm,obj.Vnm);
        end
    end
    v2d(v2d==obj.ndv)=NaN;

% Check the boundary
    validateattributes(v2d(~isnan(v2d)),{'double'},{'<=',obj.Ulm,'>=',obj.Llm},'',nm);
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
    sz=[size(X) length(obj.VHt)];
  end

%% Time line of file
  function Tutc=TimeCls(obj,Cflg)
% Convert to the UTC time line
    [~,ds,~]=cellfun(@(X) fileparts(X),obj.Fnm,'UniformOutput',false);
    k=strfind(obj.TmF{2},obj.TmF{1});
    ds=cellfun(@(X) X(k:k+length(obj.TmF{1})-1),ds,'UniformOutput',false);
    Tutc=datenum(cell2mat(ds),obj.TmF{1});
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
