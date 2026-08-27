function square = geo_square( pts, length, varargin)
    %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
    % Build a plane square
    % Input:
    %   pts - coordinates of the bottom-left corner vertex, [x y]
    %   length - length of the square
    %
    %  ---------------------------------------
    %  Please feel free to contact us with any questions! 
    %  - Xiaoxiao Du, Beihang University
    %  - duxxhf@gmail.com / duxiaoxiao@buaa.edu.cn
    %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
    
    if nargin < 2    % default parameters
        pts = [0,0];  
        length = 2.0; 
        show_plot_bool = 1;
    end
    
    % Handle Additional Geometry generation and visualization options
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
    
    
    
    x = pts(1);
    y = pts(2);
    L = length;
    % define control points and knots vector
    coefs = zeros(4,2,2);
    % I changed all the coefficients manually, so it is different that the
    % original file!!!
    coefs(:,:,1) = [-L/2, -L/2, 0, 1; L/2, -L/2, 0, 1;]';
    coefs(:,:,2) = [-L/2, L/2, 0, 1; L/2, L/2, 0, 1;]';
    knots{1} = [0 0 1 1];
    knots{2} = [0 0 1 1];
    
    % build nurbs solid by using control points and knots vector
    square = nrbmak(coefs, knots);
    
    % degeree elevate
    square = nrbdegelev(square,degelev);
    
    % insert knots
    %RefinementX = 4;    % the number of knots inseted in u direction 
    %RefinementY = 4;    % the number of knots inseted in v direction
    iuknots = 1/(RefinementX+1):1/(RefinementX+1):RefinementX/(RefinementX+1);
    ivknots = 1/(RefinementY+1):1/(RefinementY+1):RefinementY/(RefinementY+1);
    square = nrbkntins(square, {iuknots ivknots});
    
    % plot nurbs solid
    if show_plot_bool
        figure
        plot_nurbs(square, 0,1);
        axis equal;
        view(2);
    end

end

