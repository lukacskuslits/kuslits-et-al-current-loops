function [l_array,h_array]=lh_array(array,varargin)

p=inputParser;
p.addParamValue('gridpoints',@(x) isnumeric (x)); % load gridpoints if necessary
p.addParamValue('abs',5,@(x) isscalar (x) && x>0 && x<360); % set absolute changes if necessary
p.addParamValue('rel',0,@(x) isscalar (x) && x>0 && x<1); % determine relative change
p.addParamValue('domx',45,@(x) isscalar(x)); % determine x coord constraints
p.addParamValue('domy',90,@(x) isscalar(x)); % determine y coord constraints
p.addParamValue('const_l',[],@(x) isnumeric (x)); % determine lower boundaries
p.addParamValue('const_h',[],@(x) isnumeric(x)); % determine upper boundaries
p.addParamValue('rows',0,@(x) isnumeric(x)); % determine rows for the constraints to apply
p.addParamValue('mode','nul',@(x) ischar(x));
%Parse input arguments
p.parse(varargin{:});
p.Results

gridpoints=p.Results.gridpoints;
abs=p.Results.abs;
rel=p.Results.rel;
domx=p.Results.domx;
domy=p.Results.domy;
const_l=p.Results.const_l;
const_h=p.Results.const_h;
rows=p.Results.rows;
mode=p.Results.mode;

if exist('abs','var')==1
if rel~=0
l_array(1:2,:)=gridpoints(1:2,:)-abs;
l_array(length(array(:,1))+3,:)=round(gridpoints(3,:)-gridpoints(3,:)*rel);

h_array(1:2,:)=gridpoints(1:2,:)+abs;
h_array(length(array(:,1))+3,:)=round(gridpoints(3,:)+gridpoints(3,:)*rel);
else
l_array(1:2,:)=gridpoints(1:2,:)-abs;

h_array(1:2,:)=gridpoints(1:2,:)+abs;
end
else
l_array(1:2,:)=round(gridpoints(1:2,:)-gridpoints(1:2,:)*rel);
l_array(length(array(:,1))+3,:)=round(gridpoints(3,:)-gridpoints(3,:)*rel);

h_array(1:2,:)=round(gridpoints(1:2,:)+gridpoints(1:2,:)*rel);
h_array(length(array(:,1))+3,:)=round(gridpoints(3,:)+gridpoints(3,:)*rel);
end

h_array(h_array>domy) = domy;
hlat = h_array(1,:);
hlat(hlat>domx) = domx;
h_array(1,:) = hlat;

l_array(l_array<1)=1;

% disp(l_array)
% disp(h_array)

end