function [X_i, Y_i] = my_polybool(Xa,Ya,Xb,Yb)

% We assume that the coordinates are CW, and that the polygons are perfect
% boxes orthogonal to the lat lon grid:
Xa1 = min(Xa);
Xa2 = max(Xa);
Xb1 = min(Xb);
Xb2 = max(Xb);
Ya1 = min(Ya);
Ya2 = max(Ya);
Yb1 = min(Yb);
Yb2 = max(Yb);

X_i = [];
Y_i = [];


if Xa2 < Xb1 || Xb2 < Xa1 || Ya2 < Yb1 || Yb2 < Ya1
    % No overlap!
    return;
else
    % Go CW from NW corner:
    X_i = [max(Xa1,Xb1);...
        min(Xa2,Xb2);...
        min(Xa2,Xb2);...
        max(Xa1,Xb1)];
    Y_i = [min(Ya2,Yb2);...
        min(Ya2,Yb2);...
        max(Ya1,Yb1);...
        max(Ya1,Yb1)];
end

end