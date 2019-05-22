classdef Zippit < handle
%ZIPPIT A minimal class for Zipper-like algorithms
% Usage:
%   mz = zippit([0 1 1i].', 'Geodesic');
%
% Fields:
%   polygon
%       Complex vector
%       Vertices of a Jordan polygon in counterclockwise order
%   zn
%       Complex vector
%       Sequence of algorithm-required parameters
%   nv
%       Integer scalar
%       Number of polygon vertices
%   forward
%       Function handle
%       Algorithm-approximated conformal map from polygon interior to upper
%       half-plane
%   inverse
%       Function handle
%       Algorithm-approximated conformal map from upper half-plane to
%       polygon interior
%       
% Methods:
%   Zippit()
%       Construct a new Zippit instance
%   geodesic(), slit(), and zipper()
%       Run the geodesic, slit or Zipper algorithm on the polygon's vertices
%       in order to solve its parameters. Set up corresponding forward and
%       inverse maps.
%
% References:
% "Convergence of the Zipper algorithm for conformal mapping"
% Donald E. Marshall & Steffen Rohde
% http://arxiv.org/abs/math/0605532

    properties(Access = public)
        polygon
        zn
        nv
        forward
        inverse
    end

    methods
        function obj = Zippit(varargin)
            obj.polygon     = [];
            obj.zn          = [];
            obj.nv          = 0;
            
            obj.forward     = @(z)(error('Forward map is not defined.'));
            obj.inverse     = @(w)(error('Inverse map is not defined.'));
        end
        
        geodesic(obj)
        slit(obj)
        zipper(obj)
        plot(obj, ctr)
        
    end
    
end