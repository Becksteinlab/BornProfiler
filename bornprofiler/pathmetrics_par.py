import numpy as np
import logging

# NOTE: need to set axis=(1,2) to function correctly with real trajectories in
# N-dimensional space, where N = natoms*d and d = 3 for 3-dimensional real space; for
# toy trajectories, using axis=1, since trajectories are not composed of atoms.
# ---> FIXED
def d_H((P, Q)):
    '''
        Given two piecewise-linear paths path and testPath, compute the discrete
        Hausdorff distance between each path.
    '''
    axis = (1,2) if len(P.shape) == 3 else 1
    def vsqnorm(v, axis=axis):
        '''
           Compute the conventional norm of an N-dimensional vector.
        '''
        return np.sum(v*v, axis=axis)

    try:
        # alongP = np.amax( np.array([ np.amin( np.sum((pt - Q)**2, axis=axis) ) for pt in P ]) )
        # alongQ = np.amax( np.array([ np.amin( np.sum((pt - P)**2, axis=axis) ) for pt in Q ]) )
        alongP = np.max( [ np.min( vsqnorm(pt - Q, axis=axis) ) for pt in P ] )
        alongQ = np.max( [ np.min( vsqnorm(pt - P, axis=axis) ) for pt in Q ] )
        out = alongP if alongP > alongQ else alongQ
        # return ( out )**0.5
        return ( out / P.shape[1] )**0.5
    except (KeyboardInterrupt, SystemExit):
        return


def d_DF((P, Q)):
  
    p = len(P)
    q = len(Q)
    ca = -np.ones((p, q))
    
    def vsqnorm(v):
        '''
           Compute the conventional norm of an N-dimensional vector.
        '''
        return np.sum(v*v)

    P0 = P[0]
    Q0 = Q[0]
    startdist = vsqnorm(P0-Q0)
    
    def c(i, j):
        if ca[i,j] > -1 : return ca[i,j]
        elif i == 1 and j == 1 : ca[i,j] = startdist
        elif i > 1 and j == 1 : ca[i,j] = max(c(i-1,j), vsqnorm(P[i]-Q0))
        elif i == 1 and j > 1 : ca[i,j] = max(c(i,j-1), vsqnorm(P0-Q[j]))
        elif i > 1 and j > 1 : ca[i,j] = \
            max( min(c(i-1,j), c(i-1,j-1), c(i,j-1)), vsqnorm(P[i]-Q[j]))
        else: ca[i,j] = np.inf
        return ca[i,j]

    # return ( c(p-1, q-1) )**0.5
    return ( c(p-1, q-1) / P.shape[1] )**0.5


    # for i in xrange(1, p):
    #     for j in xrange(1, q):
    #         ca = max


# def linearinfimum( path, testpath, axis=(1,2) ):
#     '''
#         For each point in path compute the minimum distance to testPath. Return
#         the square of the greatest of all the minimum distances.
#     '''
#     '''
#         Option 1 --> for-loop
#         minima = np.zeros(len(path))
#         for i in xrange(0, len(path)) :
#             minima[i] = ptwiseminimum(path[i], testpath)
#         return max(minima)

#         Option 2 --> list comprehension
#     '''
#     return np.amax( np.array([ptwisemin(point, testpath, axis=axis) for point in path]) )


# def ptwisemin( point, path, axis=(1,2) ):
#     '''
#         Given a point and a piecewise-linear path, return the minimum squared
#         distance between the point and path, where distances are computed between
#         the point and each point comprising the path.
#     '''
#     '''
#         Option 1 --> for-loop
#         sqdist = np.zeros(len(path))
#         for i in xrange(0, len(path)) :
#             sqdist[i] = vsqnorm(point-path[i])
#         return np.amin(sqdist)

#         Option 2 --> list comprehension (slightly faster)
#         return np.amin( np.array([vsqnorm(point-path_pt) for path_pt in path]) )

#         Option 3 --> array broadcasting (fastest by far)
#     '''
#     return np.amin( vsqnorm(point-path, axis=axis) )


# def d_max( P, Q ):
#     '''
#         Given two piecewise-linear paths path and testPath, compute the discrete
#         maximum distance between each path.
#     '''
#     return ( max(linearmax(P, Q), linearmax(Q, P)) / len(P[0]) )**0.5


# def linearmax( path, testpath ):
#     '''
#         Determine the maximum of all point-wise maxima, i.e., for each maximum
#         value computed from a point in path to testPath.
#     '''
#     return np.amax( np.array([ptwisemax(point, testpath) for point in path]) )


# def ptwisemax( point, path ):
#     '''
#         Given a point and a piecewise-linear path, return the maximum squared
#         distance between the point and path, where distances are computed between
#         the point and each point in path.
#     '''
#     return np.amax( vsqnorm(point-path, axis=(1,2)) )


# def d_RMS( P, Q ):
#     '''
#         Given two piecewise-linear paths path and testPath, compute the discrete
#         Hausdorff distance between each path.
#     '''
#     p = len(P)
#     q = len(Q)
#     return ( (totalsum(P, Q)/(p*q)) / len(P[0]) )**0.5


# def totalsum( path, testpath ):
#     '''
#         For each point in path compute the minimum distance to testPath. Return
#         the square of the greatest of all the minimum distances.
#     '''
#     return np.sum( np.array([ptwisesum(point, testpath) for point in path]) )


# def ptwisesum( point, path ):
#     '''
#         Given a point and a piecewise-linear path, return the minimum squared
#         distance between the point and path, where distances are computed between
#         the point and each point comprising the path.
#     '''
#     return np.sum( vsqnorm(point-path, axis=(1,2)) )



def vsqnorm(v, axis=None):
    '''
        Compute the conventional norm of an N-dimensional vector.
    '''
    return np.sum(v*v, axis=axis)