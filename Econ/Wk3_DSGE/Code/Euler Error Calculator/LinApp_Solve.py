"""
MATLAB version 1.1 written by Kerk Phillips, April 2014

PYTHON version adapted by Yulong Li, November 2015 

This PYTHON version was also based on the previous adaptations of 
Uhlig's Toolkit (1999) by Spencer Lyon in May 2012 and later Chase Coleman

"""
from __future__ import division
import scipy as sp
import numpy as np
from numpy import hstack, vstack, zeros, dot, eye, kron
from scipy import linalg as la
from numpy import linalg as npla

def _nullSpaceBasis(A):
    """
    This funciton will find the basis of the null space of the matrix A.

    Parameters
    ----------
    A : array_like, dtype=float
        The matrix you want the basis for

    Returns
    -------
    vecs : array_like, dtype=float
        A numpy matrix containing the vectors as row vectors.

    Notes
    -----
    If A is an empty matrix, an empty matrix is returned.

    """
    if A.any():
        U, s, Vh = la.svd(A)
        vecs = np.array([])
        toAppend = A.shape[1] - s.size
        s = np.append(s, zeros((1, toAppend)))
        for i in range(0, s.size):
            if s[i] == 0:
                vecs = Vh[-toAppend:, :]
        if vecs.size == 0:
            vecs = zeros((1, A.shape[1]))
        return np.mat(vecs)
    else:
        return zeros((0, 0))

def qzswitch(i, A, B, Q, Z):
    '''
    Takes U.T. matrices A, B, orthonormal matrices Q,Z, interchanges
    diagonal elements i and i+1 of both A and B, while maintaining 
    Q'AZ' and Q'BZ' unchanged.  Does nothing if ratios of diagonal elements
    in A and B at i and i+1 are the same.  Aborts if diagonal elements of
    both A and B are zero at either position.

    Parameters
    ----------
    i : number, dtype=int
        Index (>=1) of the diagonal element to be interchanged
    A : array_like, dtype=float
        The upper triangular matrix of which some diagonal elements are to
        be interchanged
    B : array_like, dtype=float
        The other upper triangular matrix of which some diagonal elements are
        to be interchanged
    Q : array_like, dtype=float
        An orthonormal matrix from the QZ decomposition
    Z : array_like, dtype=float
        An orthonormal matrix from the QZ decomposition
    Returns
    -------
    A : array_like, dtype=float
        Altered A matrix
    B : array_like, dtype=float
        Altered A matrix
    Q : array_like, dtype=float
        Altered Q matrix
    Z : array_like, dtype=float
        Altered Z matrix
    Notes
    -----
    Copyright: C.A. Sims, 1996, Yale University.
    '''
    a = A[i-1,i-1]
    d = B[i-1,i-1]
    b = A[i-1,i]
    e = B[i-1,i]
    c = A[i,i]
    f = B[i,i]
  
    wz = hstack((dot(c,e)-dot(f,b), (dot(c,d)-dot(f,a)).conj().T))
    xy = hstack(((dot(b,d)-dot(e,a)).conj().T, (dot(c,d)-dot(f,a)).conj().T))

    n = np.sqrt(dot(wz,wz.conj().T))
    m = np.sqrt(dot(xy,xy.conj().T))
    
    if n == 0:
        print ("qzswitch(): Inputs unchanged!")
        return A, B, Q, Z
    else:
       wz = wz/n
       xy = xy/m
       wz = vstack(( wz, hstack((-wz[1].conj().T, wz[0].conj().T)) ))
       xy = vstack(( xy, hstack((-xy[1].conj().T, xy[0].conj().T)) ))
       A[i-1:i+1,:] = xy.dot(A[i-1:i+1,:])
       B[i-1:i+1,:] = xy.dot(B[i-1:i+1,:])
       A[:,i-1:i+1] = A[:,i-1:i+1].dot(wz)
       B[:,i-1:i+1] = B[:,i-1:i+1].dot(wz)
       Z[:,i-1:i+1] = Z[:,i-1:i+1].dot(wz)
       Q[i-1:i+1,:] = xy.dot(Q[i-1:i+1,:])
    return A, B, Q, Z

def qzdiv(stake, A, B, Q, Z):
    '''
    Takes U.T. matrices A, B, orthonormal matrices Q,Z, rearranges them
    so that all cases of abs(B(i,i)/A(i,i))>stake are in lower right 
    corner, while preserving U.T. and orthonormal properties and Q'AZ' and Q'BZ'.
    
    Parameters
    ----------
    stake : number, dtype=float
    A : array_like, dtype=float
        An upper triangular matrix
    B : array_like, dtype=float
        An upper triangular matrix
    Q : array_like, dtype=float
        An orthonormal matrix from the QZ decomposition
    Z : array_like, dtype=float
        An orthonormal matrix from the QZ decomposition
    Returns
    -------
    A : array_like, dtype=float
        Rearranged A matrix
    B : array_like, dtype=float
        Rearranged B matrix
    Q : array_like, dtype=float
        Rearranged Q matrix
    Z : array_like, dtype=float
        Rearranged Z matrix
    Notes
    -----
    Copyright: C.A. Sims, 1996, Yale University.
    '''
    n, jnk = A.shape

    root = abs(vstack((np.diag(A), np.diag(B))).T)
    tmp = (root[:,0]<1.e-13).astype(int)
    root[:,0] = root[:,0]- tmp *(root[:,0]+root[:,1])
    root[:,1] = root[:,1]/root[:,0]
    for i in range(n,0,-1):
        m=0
        for j in range(i,0,-1):
            if (root[j-1,1] > stake or root[j-1,1] < -.1):
                m=j
                break
        if m==0:
            print ("qzdiv(): Inputs unchanged!")
            return A, B, Q, Z
        for k in range(m,i,1):
            A, B, Q, Z = qzswitch(k,A,B,Q,Z)
            tmp = root[k-1,1]
            root[k-1,1] = root[k,1]
            root[k,1] = tmp
    return A, B, Q, Z


def LinApp_Solve(AA, BB, CC, DD, FF, GG, HH, JJ, KK, LL, MM, NN, Z0, Sylv):
    """
    This code takes Uhlig's original code and puts it in the form of a
    function.  This version outputs the policy function coefficients: PP,
    QQ and UU for X, and RR, SS and VV for Y.

    Inputs overview:
    The matrices of derivatives: AA - TT.
    The autoregression coefficient matrix NN from the law of motion for Z.
    Z0 is the Z-point about which the linearization is taken.  For
    linearizing about the steady state this is Zbar and normally Zbar = 0.
    Sylv is an indicator variable telling the program to use the built-in
    function sylvester() to solve for QQ and SS, if possible.  Default is
    to use Sylv=1.

    Parameters
    ----------
    AA : array_like, dtype=float, shape=(ny, nx)
        The matrix represented above by :math:`A`. It is the matrix of
        derivatives of the Y equations with repsect to :math:`X_t`
    BB : array_like, dtype=float, shape=(ny, nx)
        The matrix represented above by :math:`B`. It is the matrix of
        derivatives of the Y equations with repsect to
        :math:`X_{t-1}`.
    CC : array_like, dtype=float, shape=(ny, ny)
        The matrix represented above by :math:`C`. It is the matrix of
        derivatives of the Y equations with repsect to :math:`Y_t`
    DD : array_like, dtype=float, shape=(ny, nz)
        The matrix represented above by :math:`C`. It is the matrix of
        derivatives of the Y equations with repsect to :math:`Z_t`
    FF : array_like, dtype=float, shape=(nx, nx)
        The matrix represetned above by :math:`F`. It is the matrix of
        derivatives of the model's characterizing equations with
        respect to :math:`X_{t+1}`
    GG : array_like, dtype=float, shape=(nx, nx)
        The matrix represetned above by :math:`G`. It is the matrix of
        derivatives of the model's characterizing equations with
        respect to :math:`X_t`
    HH : array_like, dtype=float, shape=(nx, nx)
        The matrix represetned above by :math:`H`. It is the matrix of
        derivatives of the model's characterizing equations with
        respect to :math:`X_{t-1}`
    JJ : array_like, dtype=float, shape=(nx, ny)
        The matrix represetned above by :math:`J`. It is the matrix of
        derivatives of the model's characterizing equations with
        respect to :math:`Y_{t+1}`
    KK : array_like, dtype=float, shape=(nx, ny)
        The matrix represetned above by :math:`K`. It is the matrix of
        derivatives of the model's characterizing equations with
        respect to :math:`Y_t`
    LL : array_like, dtype=float, shape=(nx, nz)
        The matrix represetned above by :math:`L`. It is the matrix of
        derivatives of the model's characterizing equations with
        respect to :math:`Z_{t+1}`
    MM : array_like, dtype=float, shape=(nx, nz)
        The matrix represetned above by :math:`M`. It is the matrix of
        derivatives of the model's characterizing equations with
        respect to :math:`Z_t`
    NN : array_like, dtype=float, shape=(nz, nz)
        The autocorrelation matrix for the exogenous state vector z.
    Z0 : array, dtype=float, shape=(nz,)
        the Z-point about which the linearization is taken.  For linearizing 
        about the steady state this is Zbar and normally Zbar = 0.
        QQ if true.
    Sylv: binary, dtype=int 
        an indicator variable telling the program to use the built-in
        function sylvester() to solve for QQ and SS, if possible.  Default is
        to use Sylv=1.

    Returns
    -------
    P : 2D-array, dtype=float, shape=(nx, nx)
        The matrix :math:`P` in the law of motion for endogenous state
        variables described above.
    Q : 2D-array, dtype=float, shape=(nx, nz)
        The matrix :math:`Q` in the law of motion for exogenous state
        variables described above.
    R : 2D-array, dtype=float, shape=(ny, nx)
        The matrix :math:`R` in the law of motion for endogenous state
        variables described above.
    S : 2D-array, dtype=float, shape=(ny, nz)
        The matrix :math:`S` in the law of motion for exogenous state
        variables described above.
        
    References
    ----------
    .. [1] Uhlig, H. (1999): "A toolkit for analyzing nonlinear dynamic
       stochastic models easily," in Computational Methods for the Study
       of Dynamic Economies, ed. by R. Marimon, pp. 30-61. Oxford
       University Press.

    """
    #The original coding we did used the np.matrix form for our matrices so we
    #make sure to set our inputs to numpy matrices.
    AA = np.matrix(AA)
    BB = np.matrix(BB)
    CC = np.matrix(CC)
    DD = np.matrix(DD)
    FF = np.matrix(FF)
    GG = np.matrix(GG)
    HH = np.matrix(HH)
    JJ = np.matrix(JJ)
    KK = np.matrix(KK)
    LL = np.matrix(LL)
    MM = np.matrix(MM)
    NN = np.matrix(NN)
    Z0 = np.array(Z0)
    #Tolerance level to use
    TOL = .000001

    # Here we use matrices to get pertinent dimensions.
    nx = FF.shape[1]
    l_equ = CC.shape[0]
    ny = CC.shape[1]
    nz = min(NN.shape)

    # The following if and else blocks form the
    # Psi, Gamma, Theta Xi, Delta mats
    if l_equ == 0:
        if CC.any():
            # This blcok makes sure you don't throw an error with an empty CC.
            CC_plus = la.pinv(CC)
            CC_0 = _nullSpaceBasis(CC.T)
        else:
            CC_plus = np.mat([])
            CC_0 = np.mat([])
        Psi_mat = FF
        Gamma_mat = -GG
        Theta_mat = -HH
        Xi_mat = np.mat(vstack((hstack((Gamma_mat, Theta_mat)),
                        hstack((eye(nx), zeros((nx, nx)))))))
        Delta_mat = np.mat(vstack((hstack((Psi_mat, zeros((nx, nx)))),
                           hstack((zeros((nx, nx)), eye(nx))))))

    else:
        CC_plus = la.pinv(CC)
        CC_0 = _nullSpaceBasis(CC.T)
        if l_equ != ny:
            Psi_mat = vstack((zeros((l_equ - ny, nx)), FF \
                            - dot(dot(JJ, CC_plus), AA)))
            Gamma_mat = vstack((dot(CC_0, AA), dot(dot(JJ, CC_plus), BB) \
                        - GG + dot(dot(KK, CC_plus), AA)))
            Theta_mat = vstack((dot(CC_0, BB), dot(dot(KK, CC_plus), BB) - HH))
        else:
            CC_inv = la.inv(CC)
            Psi_mat = FF - dot(JJ.dot(CC_inv), AA)
            Gamma_mat = dot(JJ.dot(CC_inv), BB) - GG + dot(dot(KK, CC_inv), AA)
            Theta_mat = dot(KK.dot(CC_inv), BB) - HH
        Xi_mat = vstack((hstack((Gamma_mat, Theta_mat)), \
                            hstack((eye(nx), zeros((nx, nx))))))
        Delta_mat = vstack((hstack((Psi_mat, np.mat(zeros((nx, nx))))),\
                                hstack((zeros((nx, nx)), eye(nx)))))

    # Now we need the generalized eigenvalues/vectors for Xi with respect to
    # Delta. That is eVals and eVecs below.

    eVals, eVecs = la.eig(Xi_mat, Delta_mat)
    if npla.matrix_rank(eVecs) < nx:
        print("Error: Xi is not diagonalizable, stopping...")

    # From here to line 158 we Diagonalize Xi, form Lambda/Omega and find P.
    else:
        Xi_sortabs = np.sort(abs(eVals))
        Xi_sortindex = np.argsort(abs(eVals))
        Xi_sortedVec = np.array([eVecs[:, i] for i in Xi_sortindex]).T
        Xi_sortval = eVals[Xi_sortindex]
        Xi_select = np.arange(0, nx)
        if np.imag(Xi_sortval[nx - 1]).any():
            if (abs(Xi_sortval[nx - 1] - sp.conj(Xi_sortval[nx])) < TOL):
                drop_index = 1
                cond_1 = (abs(np.imag(Xi_sortval[drop_index-1])) > TOL)
                cond_2 = drop_index < nx
                while cond_1 and cond_2:
                    drop_index += 1
                if drop_index >= nx:
                    print("There is an error. Too many complex eigenvalues."
                          +" Quitting...")
                else:
                    print("Droping the lowest real eigenvalue. Beware of" +
                          " sunspots!")
                    Xi_select = np.array([np.arange(0, drop_index - 1),\
                                          np.arange(drop_index, nx + 1)])
        # Here Uhlig computes stuff if user chose "Manual roots" I skip it.
        if max(abs(Xi_sortval[Xi_select])) > 1 + TOL:
            print("It looks like we have unstable roots. This might not work...")
        if abs(max(abs(Xi_sortval[Xi_select])) - 1) < TOL:
            print("Check the model to make sure you have a unique steady" +
                  " state we are having problems with convergence.")
        Lambda_mat = np.diag(Xi_sortval[Xi_select])
        Omega_mat = Xi_sortedVec[nx:2 * nx, Xi_select]

        if npla.matrix_rank(Omega_mat) < nx:
            print("Omega matrix is not invertible, Can't solve for P; we" +
                    " proceed with QZ-method instead.")

            #~~~~~~~~~ QZ-method codes from SOLVE_QZ ~~~~~~~~#
            Delta_up,Xi_up,UUU,VVV=la.qz(Delta_mat,Xi_mat, output='complex')
            UUU=UUU.T
            Xi_eigval = np.diag( np.diag(Xi_up)/np.maximum(np.diag(Delta_up),TOL))
            Xi_sortabs= np.sort(abs(np.diag(Xi_eigval)))
            Xi_sortindex= np.argsort(abs(np.diag(Xi_eigval)))
            Xi_sortval = Xi_eigval[Xi_sortindex, Xi_sortindex]
            Xi_select = np.arange(0, nx)
            stake = max(abs(Xi_sortval[Xi_select])) + TOL

            Delta_up, Xi_up, UUU, VVV = qzdiv(stake,Delta_up,Xi_up,UUU,VVV)
                    
            #Check conditions from line 49-109
            if np.imag(Xi_sortval[nx - 1]).any():
                if (abs(Xi_sortval[nx - 1] - sp.conj(Xi_sortval[nx])) < TOL):
                    print("Problem: You have complex eigenvalues! And this means"+
                        " PP matrix will contain complex numbers by this method." )
                drop_index = 1
                cond_1 = (abs(np.imag(Xi_sortval[drop_index-1])) > TOL)
                cond_2 = drop_index < nx
                while cond_1 and cond_2:
                    drop_index += 1
                if drop_index >= nx:
                    print("There is an error. Too many complex eigenvalues."
                              +" Quitting...")
                else:
                    print("Dropping the lowest real eigenvalue. Beware of" +
                          " sunspots!")
                    for i in range(drop_index,nx+1):
                        Delta_up,Xi_up,UUU,VVV = qzswitch(i,Delta_up,Xi_up,UUU,VVV)
                    Xi_select1 = np.arange(0,drop_index-1)
                    Xi_select = np.append(Xi_select1, np.arange(drop_index,nx+1))

            if Xi_sortval[max(Xi_select)] < 1 - TOL:
                print('There are stable roots NOT used. Proceeding with the' +
                        ' smallest root.')
            if max(abs(Xi_sortval[Xi_select])) > 1 + TOL:
                print("It looks like we have unstable roots. This might not work...")
            if abs(max(abs(Xi_sortval[Xi_select])) - 1) < TOL:
                print("Check the model to make sure you have a unique steady" +
                          " state we are having problems with convergence.")
            #End of checking conditions
            #Lambda_mat = np.diag(Xi_sortval[Xi_select]) # to help sol_out.m
            
            VVV=VVV.conj().T
            VVV_2_1 = VVV[nx : 2*nx, 0 : nx]
            VVV_2_2 = VVV[nx : 2*nx, nx :2*nx]
            UUU_2_1 = UUU[nx : 2*nx, 0 : nx]
            VVV = VVV.conj().T
            
            if abs(la.det(UUU_2_1))< TOL:
                print("One necessary condition for computing P is NOT satisfied,"+
                    " but we proceed anyways...")
            if abs(la.det(VVV_2_1))< TOL:
                print("VVV_2_1 matrix, used to compute for P, is not invertible; we"+
                    " are in trouble but we proceed anyways...")
            
            PP = np.matrix( la.solve(- VVV_2_1, VVV_2_2) )
            PP_imag = np.imag(PP)
            PP = np.real(PP)
            if (sum(sum(abs(PP_imag))) / sum(sum(abs(PP))) > .000001).any():
                print("A lot of P is complex. We will continue with the" +
                      " real part and hope we don't lose too much information.")
            #~~~~~~~~~ End of QZ-method ~~~~~~~~~#

        #This follows the original uhlig.py file
        else:
            PP = dot(dot(Omega_mat, Lambda_mat), la.inv(Omega_mat))
            PP_imag = np.imag(PP)
            PP = np.real(PP)
            if (sum(sum(abs(PP_imag))) / sum(sum(abs(PP))) > .000001).any():
                print("A lot of P is complex. We will continue with the" +
                      " real part and hope we don't lose too much information.")
    # The code from here to the end was from he Uhlig file calc_qrs.m.
    # I think for python it fits better here than in a separate file.

    # The if and else below make RR and VV depending on our model's setup.
    if l_equ == 0:
        RR = zeros((0, nx))
        VV = hstack((kron(NN.T, FF) + kron(eye(nz), \
            (dot(FF, PP) + GG)), kron(NN.T, JJ) + kron(eye(nz), KK))) 

    else:
        RR = - dot(CC_plus, (dot(AA, PP) + BB))
        VV = sp.vstack((hstack((kron(eye(nz), AA), \
                        kron(eye(nz), CC))), hstack((kron(NN.T, FF) +\
                        kron(eye(nz), dot(FF, PP) + dot(JJ, RR) + GG),\
                        kron(NN.T, JJ) + kron(eye(nz), KK)))))

    # Now we use LL, NN, RR, VV to get the QQ, RR, SS, VV matrices.
    # first try using Sylvester equation solver
    if Sylv:
        if ny>0:
            PM = (FF-la.solve(JJ.dot(CC),AA))
            if npla.matrix_rank(PM)< nx+ny:
                Sylv=0
                print("Sylvester equation solver condition is not satisfied;"\
                        +" proceed with the original method...")
        else:
            if npla.matrix_rank(FF)< nx:
                Sylv=0
                print("Sylvester equation solver condition is not satisfied;"\
                        +" proceed with the original method...")
        print("Using Sylvester equation solver...")
        if ny>0:
            Anew = la.solve(PM, (FF.dot(PP)+GG+JJ.dot(RR)-\
                    la.solve(KK.dot(CC), AA)) )
            Bnew = NN
            Cnew1 = la.solve(JJ.dot(CC),DD.dot(NN))+la.solve(KK.dot(CC), DD)-\
                    LL.dot(NN)-MM
            Cnew = la.solve(PM, Cnew1)
            QQ = la.solve_sylvester(Anew,Bnew,Cnew)
            SS = la.solve(-CC, (AA.dot(QQ)+DD))
        else:
            Anew = la.solve(FF, (FF.dot(PP)+GG))
            Bnew = NN
            Cnew = la.solve(FF, (-LL.dot(NN)-MM))
            QQ = la.solve_sylvester(Anew,Bnew,Cnew)
            SS = np.zeros((0,nz)) #empty matrix
    
    # then the Uhlig's way
    else:
        '''
        # This code is from Spencer Lypn's 2012 version
        q_eqns = sp.shape(FF)[0]
        m_states = sp.shape(FF)[1]
        l_equ = sp.shape(CC)[0]
        n_endog = sp.shape(CC)[1]
        k_exog = min(sp.shape(sp.mat(NN))[0], sp.shape(sp.mat(NN))[1])
        sp.mat(LL.T)
        sp.mat(NN)
        sp.dot(sp.mat(LL),sp.mat(NN))
        LLNN_plus_MM = sp.dot(sp.mat(LL),sp.mat(NN)) + sp.mat(MM.T)
        QQSS_vec = sp.dot(la.inv(sp.mat(VV)), sp.mat(LLNN_plus_MM))
        QQSS_vec = -QQSS_vec
        if max(abs(QQSS_vec)) == sp.inf:
            print("We have issues with Q and S. Entries are undefined. Probably because V is no inverible.")
        
        QQ = sp.reshape(QQSS_vec[0:m_states*k_exog],(m_states,k_exog))
        SS = sp.reshape(QQSS_vec[(m_states*k_exog):((m_states+n_endog)*k_exog)],(n_endog,k_exog))
        
        # The vstack and hstack's below are ugly, but safe. If you have issues with WW uncomment the first definition and comment out the second one.
        #WW = sp.vstack((\
        #sp.hstack((sp.eye(m_states), sp.zeros((m_states,k_exog)))),\
        #sp.hstack((sp.dot(RR,sp.pinv(PP)), (SS-sp.dot(sp.dot(RR,sp.pinv(PP)),QQ)))),\
        #sp.hstack((sp.zeros((k_exog,m_states)),sp.eye(k_exog)))))
        
        WW = sp.array([[sp.eye(m_states), sp.zeros((m_states,k_exog))],\
        [sp.dot(RR,la.pinv(PP)), (SS-sp.dot(sp.dot(RR,la.pinv(PP)),QQ))],\
        [sp.zeros((k_exog,m_states)),sp.eye(k_exog)]])
    
        '''
        
        # this code is from Yulong Li's 2015 version
        if (npla.matrix_rank(VV) < nz * (nx + ny)):
            print("Sorry but V is not invertible. Can't solve for Q and S;"+
                     " but we proceed anyways...")
        
        LL = sp.mat(LL)
        NN = sp.mat(NN)
        LLNN_plus_MM = dot(LL, NN) + MM

        if ny>0:
            impvec = vstack([DD, LLNN_plus_MM])
        else:
            impvec = LLNN_plus_MM

        impvec = np.reshape(impvec, ((nx + ny) * nz, 1), 'F')
        
        QQSS_vec = np.matrix(la.solve(-VV, impvec))

        if (max(abs(QQSS_vec)) == sp.inf).any():
            print("We have issues with Q and S. Entries are undefined." +
                      " Probably because V is no inverible.")

        #Build QQ SS
        QQ = np.reshape(np.matrix(QQSS_vec[0:nx * nz, 0]),
                            (nx, nz), 'F')

        SS = np.reshape(QQSS_vec[(nx * nz):((nx + ny) * nz), 0],\
                            (ny, nz), 'F')
        

    return np.array(PP), np.array(QQ), np.array(RR), np.array(SS)