clear;clf;hold on

% Original data in 4D
X_Fit = rand(30,3);
plot3(X_Fit(:,1),X_Fit(:,2),X_Fit(:,3),'r.','markersize',20)

% Generate 3 layers, each delta=0.05 outside hull
[X_penalty, penalty] = convexHullLayers(X_Fit, 3, 0.05);

plot3(X_penalty(:,1),X_penalty(:,2),X_penalty(:,3),'b.','markersize',20)
axis equal
size(X_penalty)   % points generated
penalty(1:10)     % corresponding penalties




function [X_penalty, penalty] = convexHullLayers(X_Fit, numLayers, delta)
% Generate layers of points outside the convex hull with penalties.
%
% X_Fit      : N × d data points (here d = 4)
% numLayers  : number of outer layers to generate
% delta      : step size for each outward layer
%
% X_penalty  : matrix of generated outer points
% penalty    : corresponding penalty values

    d = size(X_Fit,2);       % dimension (4D here)
    K = convhulln(X_Fit);    % convex hull facets
    
    % storage
    X_penalty = [];
    penalty = [];
    
    % Loop over facets
    for i = 1:size(K,1)
        % facet vertices
        pts = X_Fit(K(i,:), :);
        
        % normal vector (null space of facet)
        N = null(pts(2:end,:) - pts(1,:));
        n = N(:,1); 
        n = n / norm(n);
        
        % check orientation: ensure outward normal
        c = mean(X_Fit,1); % centroid of hull
        if dot(n, pts(1,:) - c) < 0
            n = -n;
        end
        
        % generate layers by offsetting vertices
        for k = 1:numLayers
            shift = k * delta * n';
            newPts = pts + shift;
            
            X_penalty = [X_penalty; newPts]; %#ok<AGROW>
            penalty   = [penalty; repmat(k*delta, size(newPts,1), 1)]; %#ok<AGROW>
        end
    end
end
