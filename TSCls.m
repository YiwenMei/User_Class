classdef TSCls
% TSCls is a class for a time series pair with one or multiple vertical profies.

properties
  vtp % Type of variable
  Ulm % Physical upper limit
  Llm % Physical lower limit
  Gtg % Geographic information of the time series (Gtg.Sid Gtg.Xcor Gtg.Ycor where Sid is the
      %  name of the location and Xcor/Ycor is the horizontal/vertical coordinate)
  gid % Group ID of the file
  os1 % offset in hours to UTC
  TS1 % Values of the variable
  TL1 % Time line of the variable
  TR1 % Time resolution
  TC1 % Time-window convention
  os2 % offset in hours to UTC
  TS2 % Values of the variable
  TL2 % Time line of the variable
  TR2 % Time resolution
  TC2 % Time-window convention
  Vh1 % Height(-)/depth(+) of vertical layer
  Vh2 % Height(-)/depth(+) of vertical layer

  unt % Unit of the variable
end

methods
%% Object building
  function obj=TSCls(vtp,Ulm,Llm,Gtg,gid,os1,TS1,TL1,TR1,TC1,os2,TS2,TL2,TR2,TC2,Vh1,Vh2,varargin)
% Check inputs
    narginchk(17,18);
    ips=inputParser;
    ips.FunctionName=mfilename;

    addRequired(ips,'vtp',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'vtp'));
    addRequired(ips,'Ulm',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'Ulm'));
    addRequired(ips,'Llm',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'Llm'));
    addRequired(ips,'Gtg',@(x) validateattributes(x,{'struct'},{'nonempty'},mfilename,'Gtg'));
    addRequired(ips,'gid',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'gid'));
    addRequired(ips,'os1',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'os1'));
    addRequired(ips,'TS1',@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'TS1'));
    addRequired(ips,'TL1',@(x) validateattributes(x,{'double'},{'vector'},mfilename,'TL1'));
    addRequired(ips,'TR1',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'TR1'));
    addRequired(ips,'TC1',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'TC1'));
    addRequired(ips,'os2',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'os2'));
    addRequired(ips,'TS2',@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'TS2'));
    addRequired(ips,'TL2',@(x) validateattributes(x,{'double'},{'vector'},mfilename,'TL2'));
    addRequired(ips,'TR2',@(x) validateattributes(x,{'double','char'},{'nonempty'},mfilename,'TR2'));
    addRequired(ips,'TC2',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'TC2'));
    addRequired(ips,'Vh1',@(x) validateattributes(x,{'double'},{'vector','numel',size(TS1,2)},...
        mfilename,'Vh1'));
    addRequired(ips,'Vh2',@(x) validateattributes(x,{'double'},{'vector','numel',size(TS2,2)},...
        mfilename,'Vh2'));

    addOptional(ips,'unt','-',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'unt'));

    parse(ips,vtp,Ulm,Llm,Gtg,gid,os1,TS1,TL1,TR1,TC1,os2,TS2,TL2,TR2,TC2,Vh1,Vh2,varargin{:});
    unt=ips.Results.unt;
    clear ips varargin

% Assign values
    obj.vtp=vtp;
    obj.Ulm=Ulm;
    obj.Llm=Llm;
    obj.Gtg=Gtg;
    obj.gid=gid;
    obj.os1=os1;
    obj.TS1=TS1;
    obj.TL1=TL1;
    obj.TR1=TR1;
    obj.TC1=TC1;
    obj.os2=os2;
    obj.TS2=TS2;
    obj.TL2=TL2;
    obj.TR2=TR2;
    obj.TC2=TC2;
    obj.Vh1=Vh1;
    obj.Vh2=Vh2;

    obj.unt=unt;

% Check the boundary
    ts=[reshape(obj.TS1,numel(obj.TS1),1);reshape(obj.TS2,numel(obj.TS2),1)];
    ts(isnan(ts))=[];
    validateattributes(ts,{'double'},{'<=',obj.Ulm,'>=',obj.Llm,'nonnan'});
  end

%% Unify the time line
  function [TL1,TL2]=UniTL(obj,CTp)
% Align the time series
    switch CTp
      case 'UTC'
        TL1=obj.TL1+obj.os1/24;
        TL2=obj.TL2+obj.os2/24;
      case 'tg'
        TL1=obj.TL1;
        TL2=obj.TL2+obj.os2/24-obj.os1/24;
      case 'rf'
        TL1=obj.TL1+obj.os1/24-obj.os2/24;
        TL2=obj.TL2;
    end

% Convert to center time convention
    if ischar(obj.TR1)
      obj.TR1=0;
    end
    if ischar(obj.TR2)
      obj.TR2=0;
    end

    if strcmp(obj.TC1,'begin')
      TL1=obj.TL1+obj.TR1/24/2;
    elseif strcmp(obj.TC1,'end')
      TL1=obj.TL1-obj.TR1/24/2;
    end

    if strcmp(obj.TC2,'begin')
      TL2=obj.TL2+obj.TR2/24/2;
    elseif strcmp(obj.TC2,'end')
      TL2=obj.TL2-obj.TR2/24/2;
    end
  end
end
end
