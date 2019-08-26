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

% Mark Kittisopikul, October 2018
% Goldman Lab
% Northwestern University

    [varargout{1:nargout}] = deal(C{:});
    
end

