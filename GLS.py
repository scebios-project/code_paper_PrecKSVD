import numpy
import scipy
import scipy.linalg

def greedy_least_squares(data, dictionary, stopval):
    """
    Greedy Least Squares Pursuit algorihtm
    :param data: 2D array containing the data to decompose, columnwise
    :param dictionary: dictionary containing the atoms, columnwise
    :param stopval: stopping criterion
    :return: coefficients
    """

    if len(data.shape) == 1:
        data = numpy.atleast_2d(data)
        if data.shape[0] < data.shape[1]:
            data = numpy.transpose(data)
    coef         = numpy.zeros((dictionary.shape[1], data.shape[1]))
    coef_support = numpy.zeros((dictionary.shape[1], data.shape[1]))
    coef_debias  = numpy.zeros((dictionary.shape[1], data.shape[1]))
    final_error  = numpy.zeros((1,data.shape[1]))

    # Prepare the pseudoinverse of the dictionary, used for all signals
    Dpinv_original = scipy.linalg.pinv(dictionary)    # tall matrix
    
    for i in range(data.shape[1]):
        
        # Solve GLSP for single vector data[i]
        y = data[:,i]
        N = dictionary.shape[1]
        T = numpy.array((), dtype=int)      # Set of selected atoms
        Tc = numpy.arange(N, dtype =int)    # Set of remaining atoms
        gamma = scipy.linalg.lstsq(dictionary,y)[0]
        Dpinv = Dpinv_original.copy()       # make fresh copy, Dpinv will be destroyed in the process
        
        niter = 0
        finished = False
        while not finished:
            
            # 1.Select new atom
            # Make a copy of gamma and zero the chosen atoms, 
            # so the selected atom preserves the original numbering
            gammatemp = gamma.copy()
            gammatemp[T] = 0    # skip already chosen elements
            sel = numpy.argmax(numpy.abs(gammatemp))
            # Update sets
            T = numpy.append(T,sel)
            Tc = numpy.setdiff1d(Tc, [sel], assume_unique=False)
            
            #print 'GLSP step %d: selected atom %d'%(niter, sel)
            
            # 2. Update gamma
            # Compute new null space direction and make column vector 2D
            newNS = numpy.dot(Dpinv, dictionary[:,sel])
            newNS = numpy.atleast_2d(newNS)
            if newNS.shape[0] < newNS.shape[1]:
                newNS = newNS.T
            # Update gamma
            newNSProj = numpy.dot(newNS, newNS.T) / numpy.linalg.norm(newNS, 2)**2
            gamma = gamma - numpy.dot(newNSProj, gamma)
            #gamma = gamma - numpy.dot(newNS, numpy.dot(newNS.T, gamma)) / numpy.dot(newNS.T, newNS)
            #gamma = numpy.squeeze(gamma) # make vector again, for easier access to entries later
            # Update Dpinv (QR orthogonalization)
            Dpinv = Dpinv - numpy.dot(newNSProj, Dpinv)
            
            # Check against whole NS calculation
            #gamma1 = solveGammaWholeNS(y, dictionary, T)
            # Check against system solve GAP-style
            #gamma2 = solveGammaGAPSystem(y, dictionary, T)
            
            niter = niter + 1            
            
            # 3. Check termination conditions
            if (stopval < 1) and (numpy.linalg.norm(gamma[Tc],2) < stopval):
                finished = True
            elif (stopval >= 1) and (niter == stopval):
                finished = True
         
        #print 'GLSP final: T = %s'%(T)
        
        # Prepare outputs for i'th data vector
        # Recompute gamma on support T, since our algorithm only preserved values in Tc (later edit: true, double-checked)
        coef[:,i] = gamma.copy()  # only the Tc coefficients are in gamma
        coef[T,i] = scipy.linalg.lstsq(dictionary[:,T], y-dictionary@gamma)[0]
        coef_support[T,i] = 1
        final_error[0,i] = numpy.linalg.norm(gamma[Tc],2) # final remaining norm of the cosupport atoms
        # Ouptut debiased solution as well

        gamma[T] = scipy.linalg.lstsq(dictionary[:,T], y)[0]
        gamma[Tc] = 0
        coef_debias[:,i] = gamma
    
    return coef, coef_support, coef_debias, final_error
