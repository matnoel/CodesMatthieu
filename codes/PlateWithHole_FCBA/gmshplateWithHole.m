function [model] = gmshplateWithHole(domain, circle, clD, clC, split, optimesh, pathname)

filename = 'plateWithHole';
coefZone = 5;

if not(optimesh)

    model = gmshdomainwithhole(domain,circle,clD,clC,fullfile(pathname,filename));
    
else
    
    % Recovered pts from the domain
    PtsDomaine = getvertices(domain);
    
    % Calculates the height and width of the domain
    height = PtsDomaine{3}(2);
    width = PtsDomaine{2}(1)*2;    
    
    % Recover pts from circle
    PtsCercle = getvertices(circle);
    
    % Get the diameter and vertical position of the hole
    holeDiameter = abs(PtsCercle{1}(1)-PtsCercle{3}(1));
    holeHeight = abs(PtsCercle{1}(2)-height);

%     clD = height/30; % Mesh sizes for undamaged area

    % Creation of all the coordinates of the points
    P{1} = [-width/2,0];
    P{2} = [-holeDiameter*coefZone/2,0];
    P{3} = [holeDiameter*coefZone/2,0];
    P{4} = [width/2,0];
    P{5} = [width/2,height-holeHeight-holeDiameter*coefZone/2];
    P{6} = [width/2,height-holeHeight+holeDiameter*coefZone/2];
    P{7} = [width/2,height];
    P{8} = [holeDiameter*coefZone/2,height];
    P{9} = [-holeDiameter*coefZone/2,height];
    P{10} = [-width/2,height];
    P{11} = [-width/2,height-holeHeight+holeDiameter*coefZone/2];
    P{12} = [-width/2,height-holeHeight-holeDiameter*coefZone/2];
    
    % Initialize the Geo pt in which we will write
    G = GMSHFILE();
    
    % Create the points and then the lines in the .geo
    if split == "AnisotropicHe" || split == "AnisotropicMiehe"
        % Vertical zone
        for i=[1,4,5,6,7,10,11,12]
            G = createpoints(G,P{i},clD,i);
        end
        for i=[2,3,8,9]
            G = createpoints(G,P{i},clC,i);
        end

        % Create lines
        points1 = [1,2,9,10,11,12];   lignes1 = [12,29,910,1011,1112,121];
        points2 = [2,3,8,9];      lignes2 = [23,38,89,92];
        points3 = [3,4,5,6,7,8];  lignes3 = [34,45,56,67,78,83];
    else
        % Horizontal zone
        for i=[1,2,3,4,7,8,9,10]
            G = createpoints(G,P{i},clD,i);
        end
        for i=[5,6,11,12]
            G = createpoints(G,P{i},clC,i);
        end
        
        % Create lines
        points1 = [1,2,3,4,5,12];   lignes1 = [12,23,34,45,512,121];
        points2 = [5,6,11,12];  lignes2 = [56,611,1112,125];
        points3 = [6,7,8,9,10,11];  lignes3 = [67,78,89,910,1011,116];
    end
    
    % Create the 3 outlines
    G = createcontour(G,points1,lignes1,1);
    G = createcontour(G,points2,lignes2,2);
    G = createcontour(G,points3,lignes3,3);
    
    % Create the circle
    circle = CIRCLE(0.0,height-holeHeight,holeDiameter/2);
    nCenter = 13;
    ptsCercle = 14:17;
    lignesCercle = [1415,1516,1617,1714];
    
    % .geo of the circle
    GC = gmshfile(circle,clC,nCenter,ptsCercle,lignesCercle,2);
    
    % Merges the 2 .geo's
    G = G+GC;
    
    % Create the loop that contains the hole
    G = createlineloop(G,[lignes2,-lignesCercle],22);
    
    % Create the 3 planes
    G = createplanesurface(G,1,1);
    G = createplanesurface(G,22,2);
    G = createplanesurface(G,3,3);

    % Creer une surface avec les 3 plans
    G = createphysicalsurface(G,1:3,1);

    % Builds .geo and .msh
    G = setfile(G,fullfile(pathname,filename));
    
    % Builds the model
    model = gmsh2femobject(2,G);
    
    % Combines the 3 surfaces
    model = concatgroupelem(model);
end