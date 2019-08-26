function logic= isvalid(aDouble)
% Dirty trick to handle isvalid absence in <2013b
global debuggingMode__;
if(debuggingMode__)
    warning(['isvalid:' class(aDouble) ':undefined'], ...
        ['isvalid is not defined for ' class(aDouble) ...
        '. Using always true function at ' mfilename('fullpath') '.']);
end
logic=true;
