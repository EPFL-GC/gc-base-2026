import igl
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components, depth_first_order
import triangle as tr
import warnings

def triangulate_mesh(Ps, Holes, area_density):
    '''
    Triangulates a set of points given by P using Delaunay scheme

    Input:
    - P : list of N lists of 2 elements giving the vertices positions
    - Holes: list of M lists of 2 elements giving the holes positions
    - area_density: the maximum triangle area

    Output:
    - V : array of shape (nNodes, 3) containing the vertices position
    - F : array of vertices ids (list of length nTri of lists of length 3). The ids are such that the algebraic area of
          each triangle is positively oriented
    '''

    num_points = 0
    points = []
    segments = []
    for k, _ in enumerate(Ps):
        P = np.array(Ps[k])
        N = P.shape[0]
        points.append(P)
        index = np.arange(N)
        seg = np.stack([index, index + 1], axis=1) % N + num_points
        segments.append(seg)
        num_points += N

    points = np.vstack(points)
    segments = np.vstack(segments)

    data = []
    if Holes == [] or Holes == [[]] or Holes == None:
        data = dict(vertices=points, segments=segments)
    else:
        data = dict(vertices=points, segments=segments, holes = Holes)

    tri = tr.triangulate(data, 'qpa{}'.format(area_density))

    V = np.array([[v[0], v[1], 0] for v in tri["vertices"]])
    F = np.array([[f[0], f[1], f[2]] for f in tri["triangles"]])
    return [V, F]

def simplify_boundary_add_slit(scaled_V, F, lengthSlit=2.9, heightSlit=5.0, chamfer=0.5):
    '''
    Simplifies the boundary of the mesh by removing the vertices that are on the contact line, then add a slit for each contact line
    
    Args:
        scaled_V: the scaled vertices of the mesh
        F: the faces of the mesh
        lengthSlit: the length of the slit [mm]
        heightSlit: the height of the slit [mm]
        chamfer: the chamfer of the slit [mm]
        
    Returns:
        newScaledV: the new scaled vertices of the mesh
        simplifiedBoundaryFacets: the new boundary facets of the mesh
    '''
    
    # First find the contact lines
    contactPaths = []
    tol = 1.0e-6
    for Es in igl.boundary_facets(F)[0].tolist():
        
        if abs(scaled_V[Es[0], 1]) < tol and abs(scaled_V[Es[1], 1]) < tol:
            es0IdPath = None
            es0AtBeginning = None
            es1IdPath = None
            es1AtBeginning = None
            for idPath, path in enumerate(contactPaths):
                if Es[0] in path:
                    es0IdPath = idPath
                    es0AtBeginning = Es[0] == path[0]
                if Es[1] in path:
                    es1IdPath = idPath
                    es1AtBeginning = Es[1] == path[0]
            # Merge the two paths
            if es0IdPath is not None and es1IdPath is not None:
                if es0AtBeginning:
                    if es1AtBeginning:
                        newPath = contactPaths[es0IdPath][::-1] + contactPaths[es1IdPath] 
                    else:
                        newPath = contactPaths[es1IdPath] + contactPaths[es0IdPath] 
                else:
                    if es1AtBeginning:
                        newPath = contactPaths[es0IdPath] + contactPaths[es1IdPath] 
                    else:
                        newPath = contactPaths[es0IdPath] + contactPaths[es1IdPath][::-1]
                contactPaths.pop(max(es0IdPath, es1IdPath))
                contactPaths.pop(min(es0IdPath, es1IdPath))
            # Append to either paths
            elif es0IdPath is not None:
                if es0AtBeginning:
                    newPath = [Es[1]] + contactPaths[es0IdPath]
                else:
                    newPath = contactPaths[es0IdPath] + [Es[1]]
                contactPaths.pop(es0IdPath)
            elif es1IdPath is not None:
                if es1AtBeginning:
                    newPath = [Es[0]] + contactPaths[es1IdPath]
                else:
                    newPath = contactPaths[es1IdPath] + [Es[0]]
                contactPaths.pop(es1IdPath)
            # Add a new path
            else:
                newPath = Es.copy()
            contactPaths.append(newPath)
    
    # For each contact path, simplify first, then add the slit
    nVertices = scaled_V.shape[0]
    simplifiedContactPaths = [[path[0], path[-1]] for path in contactPaths]
    newScaledV = scaled_V.copy()

    for i, path in enumerate(simplifiedContactPaths):
        if scaled_V[path[0], 0] > scaled_V[path[1], 0]:
            path = path[::-1]
        simplifiedContactPaths[i] = [path[0]] + list(range(nVertices, nVertices+4)) + [path[1]]
        assert (scaled_V[path[1], 0] - scaled_V[path[0], 0]) > lengthSlit, "The slit is too long {:.2f} (line length) < {:.2f} (slit width)".format(scaled_V[path[1], 0] - scaled_V[path[0], 0], lengthSlit)
        if scaled_V[path[1], 0] - scaled_V[path[0], 0] < lengthSlit + 6.0:
            warnings.warn("The remaining contact line is too small on either side {:.2f}mm (remaining length) < 3mm (recommended)".format((scaled_V[path[1], 0] - scaled_V[path[0], 0] - lengthSlit) / 2.0))
        centerContact = np.mean(scaled_V[path], axis=0)
        newPoints = np.stack([
            centerContact + np.array([-lengthSlit/2, 0.0, 0.0]),
            centerContact + np.array([-lengthSlit/2 + chamfer, heightSlit, 0.0]),
            centerContact + np.array([lengthSlit/2 - chamfer, heightSlit, 0.0]),
            centerContact + np.array([lengthSlit/2, 0.0, 0.0]),
        ], axis=0)
        newScaledV = np.concatenate([newScaledV, newPoints], axis=0)
        nVertices += 4
        
    # Get the new boundary facets
    simplifiedBoundaryFacets = []
    for Es in igl.boundary_facets(F)[0].tolist():
        if abs(scaled_V[Es[0], 1]) > tol or abs(scaled_V[Es[1], 1]) > tol:
            simplifiedBoundaryFacets.append(Es)

    simplifiedBoundaryFacets = simplifiedBoundaryFacets + [[path[i], path[i+1]] for path in simplifiedContactPaths for i in range(len(path)-1)]
    
    return newScaledV, simplifiedBoundaryFacets

def edges_to_closed_path(vertices, edges):
    # get boundary edges in the correct, "sorted" order by doing a DFS graph traversal
    E = np.array(edges)
    N = vertices.shape[0]
    entries = np.ones(E.shape[0]), (E[:,0], E[:,1])
    graph = coo_matrix(entries, shape=(N, N)).tocsr()
    numLoops, vtxLoopLabels = connected_components(graph, directed=False, return_labels=True)

    paths = []

    for loopIdx in range(numLoops):
        startVtxIdx = np.argmax(vtxLoopLabels >= loopIdx)
        # get boundary edges in the correct, "sorted" order by doing a DFS graph traversal
        bV = depth_first_order(graph, startVtxIdx, directed=False, return_predecessors=False)
        if len(bV) < 2: continue
        bV = bV.tolist() + [bV[0]]
        paths.append(bV)

    return paths
