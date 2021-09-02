% Yiwen Mei (ymei2@gmu.edu)
% CEIE, George Mason University
% Last update: 10/5/2019

%% Functionality
% This function convert geo-referenced 2D matlab variable to geotif.

%% Input
% tfn : full name of the output geotiff file;
% matV: input matlab variable;
% ndv : no data value of the image;
% obj : geo-reference information supplied as a Matlab structure or a file (for
%        the case of a Matlab structure, it must have three fields - xll, yll, and
%        rs - for the lower-left corner coordinates and the resolution; for the
%        case of a file, this must be the filename of a geotiff file; 
% ors : coordinate system of the image;
% wpth: working directory of the code.

%% Output
% The output image is tfn.

function matV2tif(tfn,matV,ndv,obj,ors,wpth)
%% Check inputs
narginchk(6,6);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'tfn',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'tfn'));
addRequired(ips,'matV',@(x) validateattributes(x,{'double'},{'nonempty'},mfilename,'matV'));
addRequired(ips,'ndv',@(x) validateattributes(x,{'double'},{'scalar'},mfilename,'ndv'));
addRequired(ips,'obj',@(x) validateattributes(x,{'struct','char'},{'nonempty'},mfilename,'obj'));
addRequired(ips,'ors',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'ors'));
addRequired(ips,'wpth',@(x) validateattributes(x,{'char'},{'nonempty'},mfilename,'wpth'));

parse(ips,tfn,matV,ndv,obj,ors,wpth);
clear ips varargin

%% Extract geoinfo
if isa(obj,'struct')
  xll=obj.xll;
  yll=obj.yll;
  rs=obj.rs;
elseif isa(obj,'char')
  [~,~,ext]=fileparts(obj);
  if strcmp(ext,'.tif')
    I=geotiffinfo(obj);
    xll=I.BoundingBox(1,1);
    yll=I.BoundingBox(1,2);
    rs=I.PixelScale;
    if abs((rs(1)-rs(2))/mean(rs(1:2)))<.001
      rs=rs(1);
    else
      error('X and Y resolution must be the same');
    end
  else
    error('Geo-reference file must be geotiff (.tif)');
  end
end

%% Write Matlab variable to .asc
[~,nm,~]=fileparts(tfn);
afn=fullfile(wpth,sprintf('%s.asc',nm));
fid=fopen(afn,'w');
fprintf(fid,'%s\n%s\n%s\n%s\n%s\n%s\n',['ncols ' num2str(size(matV,2))],...
    ['nrows ' num2str(size(matV,1))],['xllcorner ' num2str(xll,12)],['yllcorner '...
    num2str(yll,12)],['cellsize ' num2str(rs)],['NODATA_value ' num2str(ndv)]);
if ~isempty(find(isnan(matV), 1))
  matV(isnan(matV))=ndv; % Assign ndv to NaN
end
dlmwrite(afn,matV,'delimiter',' ','-append');
fclose(fid);

%% Convert .asc to geotiff
fun='gdal_translate -q';
pr1=sprintf('-a_srs %s',ors);
pr2=sprintf('-a_nodata %i',ndv);
cmstr=sprintf('%s %s %s "%s" "%s"',fun,pr1,pr2,afn,tfn);
system(cmstr);
delete(afn);
end
