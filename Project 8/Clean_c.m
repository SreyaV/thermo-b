function cout = Clean_c(cin)
% Test for errors and clean up a composition.
% Complain and return all zeros if it is not viable.
% C.F. Edwards, 2/1/09

global toler

% Set the output to zeros--to be returned if unsuccessful.
cout = zeros(1,length(cin));

% Check to be sure the composition is valid.
if(sum(isnan(cin)) ~= 0)
    cin
    disp('NaN mole fraction in Clean_c')
    return
elseif(max(cin) > 1)
    cin
    disp('Mole fraction exceeds unity in Clean_c')
    return
elseif(min(cin) < 0)
    cin
    disp('Negative mole fraction in Clean_c')
    return
elseif(abs(sum(cin)-1) > toler)
    cin
    disp('Sum of mole fractions is non-unity in Clean_c')
    return
end

% If we got here, we are basically OK.  Just clean up a bit.
cout = cin/sum(cin);    % Remove any trace of nonunity sum.
cout = real(cout);      % Remove any trace of imaginary.
