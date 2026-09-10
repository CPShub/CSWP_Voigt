function square = geo_square(pts, length, varargin)
% This function generates a square geometry.
% Input:
    % pts           - (1,2) vector describing the center location
    % length        - length of the square
% Output:
    % square        - geometry structure with defined nurb vectors
% ------------------------------------------------------------------------
% Copyright (C) 2026 Tobias Henkels and Juan C. Alzate Cobo.
%
% This code is an extension and modification of the NLIGA framework
% originally developed by Du et al. (2020).
%
% ------------------------------------------------------------------------
% CITATION:
% If you use this code for your research, please cite:
%
% (1) J.C. Alzate Cobo, T. Henkels and O. Weeger, "The cross-sectional
% warping problem for hyperelastic beams: A compact formulation in
% Voigt notation", DOI: 10.48550/arXiv.2604.12886
% (2) X. Du, G. Zhao, W. Wang, M. Guo, R. Zhang, J. Yang, "NLIGA: A MATLAB
% framework for nonlinear isogeometric analysis", Computer Aided
% Geometric Design, 80, 101869, 2020.
% https://doi.org/10.1016/j.cagd.2020.101869
% ------------------------------------------------------------------------
% LICENSE:
% This function is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. (GPL-3.0-or-later)
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
% for more details.
% ------------------------------------------------------------------------
% CONTACT:
% - Tobias Henkels (tobias.henkels@stud.tu-darmstadt.de)
% - Juan C. Alzate Cobo (alzate@cps.tu-darmstadt.de)
% Technische Universität Darmstadt, Germany
% ------------------------------------------------------------------------

    if nargin < 2    % default parameters
        pts = [0,0];
        length = 2.0;
        show_plot_bool = 1;
    end

    % Handle additional geometry generation and visualization options
    if ~isempty(varargin)
        options = varargin{1};

        if isfield(options, "RefinementX")
            RefinementX = options.RefinementX;
        else
            RefinementX = 4;
        end

        if isfield(options, "RefinementY")
            RefinementY = options.RefinementY;
        else
            RefinementY = 4;
        end

        if isfield(options, "degelev")
            degelev = options.degelev;
        else
            degelev = [2,2];
        end

        if isfield(options, "show_plot")
            show_plot_bool = options.show_plot;
        else
            show_plot_bool = 0;
        end
    else
        RefinementX = 4;
        RefinementY = 4;
        degelev = [2,2];
        show_plot_bool = 0;
    end

    L = length;

    % Define control points and knot vectors
    coefs = zeros(4,2,2);

    coefs(:,:,1) = [pts(1)-L/2, pts(2)-L/2, 0, 1; ...
                    pts(1)+L/2, pts(2)-L/2, 0, 1]';
    coefs(:,:,2) = [pts(1)-L/2, pts(2)+L/2, 0, 1; ...
                    pts(1)+L/2, pts(2)+L/2, 0, 1]';

    knots{1} = [0 0 1 1];
    knots{2} = [0 0 1 1];

    % Build NURBS geometry using the control points and knot vectors
    square = nrbmak(coefs, knots);

    % Degree elevation
    square = nrbdegelev(square,degelev);

    % Insert knots
    iuknots = 1/(RefinementX+1):1/(RefinementX+1):RefinementX/(RefinementX+1);
    ivknots = 1/(RefinementY+1):1/(RefinementY+1):RefinementY/(RefinementY+1);
    square = nrbkntins(square, {iuknots ivknots});

    % Plot NURBS geometry
    if show_plot_bool
        figure
        plot_nurbs(square, 0,1);
        axis equal;
        view(2);
    end

end