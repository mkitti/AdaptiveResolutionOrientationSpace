function [ varargout ] = dealcell( C )
%dealcell Deals a cell to the output
%
% INPUT
% C - a cell array
%
% OUTPUT
% a comma separated array list
%
% EXAMPLE
%
% mystruct(5) = struct();
% [mystruct.a] = dealcell(num2cell(1:5));
%
% Copyright (C) 2019, Jaqaman Lab - UT Southwestern, Goldman Lab - Northwestern 
%
% This file is part of AdaptiveResolutionOrientationSpace.
% 
% AdaptiveResolutionOrientationSpace is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% AdaptiveResolutionOrientationSpace is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with AdaptiveResolutionOrientationSpace.  If not, see <http://www.gnu.org/licenses/>.
% 
% 

% Mark Kittisopikul, October 2018
% Goldman Lab
% Northwestern University

    [varargout{1:nargout}] = deal(C{:});
    
end

