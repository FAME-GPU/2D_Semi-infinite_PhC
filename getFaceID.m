function faceID = getFaceID(pdeModel,refCoord)
% getfaceID(pdeModel,refCoord)
%   getFaceID(pdeModel,refCoord) gives the ID of the Face closest to the input coordinate.
%   	pdeModel should have a valid Discrete Geometry and Mesh associated with it.
%   	refCoord is the input coordinate. 
%          For 3D geometry, three coordinates should be provided. [x y z] 
%          For 2D geometry, two coordinates should be provided. [x y]
%
%   Example 1: Find Face ID of a multi-cuboid object
%   	pdem = createpde;
%   	pdem.Geometry = multicuboid(1,1,1);
%       pdem.generateMesh();
%   	coords = [0 0 0]; 
%   	a=getFaceID(pdem,coords);

%   Copyright 2019 The Mathworks, Inc.

    % Get the total number Faces of the geometry
    numFaces = pdeModel.Geometry.NumFaces;
    % Get closest Mesh Node to the coordinate of choice
    closestNode = pdeModel.Mesh.findNodes('nearest',refCoord');

    % Split the mesh node IDs based on Face numbers. Each row represents a
    % Face ID. 
    for ii = 1:numFaces
        fList = pdeModel.Mesh.findNodes('region','Face',ii);
        if (find(fList == closestNode))
            faceID = ii;
            break;
        end
    end

end