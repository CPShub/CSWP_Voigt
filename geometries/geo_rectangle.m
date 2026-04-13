function square = geo_rectangle(pts, lengthA, lengthB, varargin)
    % set default values
    if nargin < 3
        pts = [0,0];  
        lengthA = 1.0; 
        lengthB = 1.0;
    end

    % Initialisiere Matrix: 4 Koordinaten (x,y,z,w) x 2 Punkte U x 2 Punkte V
    coefs = zeros(4,2,2);

    % Eckpunkte berechnen (Zentriert um pts)
    x1 = pts(1) - lengthA/2; x2 = pts(1) + lengthA/2;
    y1 = pts(2) - lengthB/2; y2 = pts(2) + lengthB/2;

    % V-Richtung 1 (Unten)
    coefs(:,1,1) = [x1; y1; 0; 1]; % Unten links
    coefs(:,2,1) = [x2; y1; 0; 1]; % Unten rechts

    % V-Richtung 2 (Oben)
    coefs(:,1,2) = [x1; y2; 0; 1]; % Oben links
    coefs(:,2,2) = [x2; y2; 0; 1]; % Oben rechts

    knots{1} = [0 0 1 1];
    knots{2} = [0 0 1 1];

    % NURBS Fläche erstellen
    square = nrbmak(coefs, knots);

    % Grad erhöhen (von linear auf quadratisch)
    square = nrbdegelev(square, [2, 2]);

    % Verfeinerung (Knoten einfügen)
    RefinementX = 4;
    RefinementY = 4;
    iuknots = (1:RefinementX) / (RefinementX + 1);
    ivknots = (1:RefinementY) / (RefinementY + 1);
    square = nrbkntins(square, {iuknots, ivknots});

    % Plotten
    nrbplot(square, [20, 20]);
    axis equal;
    view(2);
end