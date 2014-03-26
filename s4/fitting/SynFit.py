import s4
import numpy as np
import re
from itertools import product 



def break_args(*args):
    """
    Break *args* and return a dictionary with parameters as keys and 
    iterator vector as values. 
    
    Parameters
    ----------
    
    parameter: str;
        The parameter in which will be looped.
        
    min_value: float, int, str;
        Minimum value.
    
    max_value: float, int, str;
        Maximum value.
        
    step: float, int, str;
        Step for each loop.
        
    Return
    ------
    
    dic: dict;
        Dictionary with parameter as key and a vector of each value
        of the loop as the dictionary value.    
    
    
    *args* is a variable length argument,
    allowing for multiple sets of parameters, inferior limit, superior 
    limit and step. For example, to return a dictionary with parameter
    *teff*, minimum value fo 15000 K, maximum value of 20000 K and 
    1000K as step::

      break_args('teff', 15000, 20000, 1000)

    An arbitrary number of *parameters*, *max_value*, *min_value* and 
    *step* groups can be specified::

      break_args(param1, min_value1, max_value1, step1, param2, 
                 min_value2, max_value2)
                 
    In case of abundance, only supports one chemical element.
    """
    
    dic = {}
       
    if len(args) < 4:
        raise TypeError, 'It should be at least 4 input parameters.'
    else:
        for i in range(len(args)/4):
            # Support for abundance iteration
            if args[4*i] == 'abund':
                ptrn = r'(\d+)[,\s]+\1[,\s]+(\d+\.\d+|\d+)'
                #check args
                for j in args[4*i+1:4*i+3]:
                    if not re.search(ptrn, j):
                        raise TypeError, 'Abundance %s has the wrong format' % j
                # Get min, max and step        
                atomic_number = re.findall(ptrn, args[4*i+1])[0][0]
                min_value = float(re.findall(ptrn, args[4*i+1])[0][-1])
                max_value = float(re.findall(ptrn, args[4*i+2])[0][-1])
                
            else: # all other cases
                min_value = float(args[4*i+1])
                max_value = float(args[4*i+2])
                
            step = float(args[4*i+3])
            n_values = np.rint((max_value - min_value)/step + 1)
            if args[4*i] == 'abund':
                vector = ['[{}, {}, {:.2f}]'.format(atomic_number, atomic_number, val)
                          for val in np.linspace(min_value, max_value, n_values)]
            else:
                vector = np.linspace(min_value, max_value, n_values)
            
            dic[args[4*i]] = vector
                
    return dic




def iterator(*args):
    """Create the iterator vector."""
    
    if len(args) == 1:
        return list(*args)
    else:
        return list(product(*args))





def synfit_sketch(*args, **kwargs):
    """
    Fit a spectral line by iterating on user defined parameter 
    an returns the best fit by minimizing the $\chi^2$ of the 
    synthetic spectrum and the observed spectrum.
    
    Parameters
    ----------
    
    synplot_path: str (optional);
        Path to synplot and synspec. 
        If not set, it wil use the `S4` default.
        
    args: float, int, str;
        Parameter to be fitted:
        
        parameter: str;
            The parameter in which will be looped.
            
        min_value: float, int, str;
            Minimum value.
        
        max_value: float, int, str;
            Maximum value.
            
        step: float, int, str;
            Step for each loop.
            
    kwargs: float, int, str;
        Any Synplot parameters, including `teff` and `logg`, 
        `synplot_path`, `idl`, `noplot`

        You MUST give the data as a Nx2 numpy array with the
        keyword argument 'observ'
        
    """
    
    # Check if there is a an observed spectrum.
    # If not quit.
    if 'observ' not in kwargs:
        raise IOError, 'There is not any observed spectrum.'
    
    #Prepare kwargs
    if 'synplot_path' in kwargs:
        synplot_path = kwargs.pop('synplot_path')
    else:
        synplot_path = None
        
    if 'idl' in kwargs:
        idl = kwargs.pop('idl')
    else:
        idl = True
        
    if 'noplot' in kwargs:
        noplot = True
        del kwargs['noplot']
    else:
        noplot = False
    
    #Obtain the parameters and vectors for the loop.
    fit_params = break_args(*args)
    
    #Obtain the number of varying params
    n_params = len(fit_params)
    
    # Obtain teff and logg if set on kwargs
    if 'teff' in kwargs:
        teff = kwargs.pop('teff')
        
    if 'logg' in kwargs:
        logg = kwargs.pop('logg')
    
    # Create the list to iterate
    iterate_list = iterator(*fit_params.values())
    
    # array to store the values of each parameter an the chisquare
    mdtype = []
    for key, value in fit_params.iteritems():
        if type(value[0]) is type(''):
            mdtype.append((key, 'S14'))
        else:
            mdtype.append((key, type(value[0])))
    mdtype.append(('chisquare',float))
   
    shape = np.shape(iterate_list)

    store_values = np.ones([shape[0], 1], dtype=mdtype)
      
    # Loop it!
    for n, it in enumerate(iterate_list):

        if n_params == 1:
            # this case correspond when it is fitting only one parameter
            syn_params = {fit_params.keys()[0]:it}
        else:
            syn_params = {key:val for key, val in zip(fit_params, it)}
            
        #make plot title before removing teff and logg
        plot_title = ', '.join(['{}={}'.format(key, val) 
                                 for key, val in syn_params.iteritems()])
        
        # Check if teff and logg were selected to be fitted.
        # If yes, set a variable to them.
        if 'teff' in syn_params:
            teff = syn_params.pop('teff')
            
        if 'logg' in syn_params:
            logg = syn_params.pop('logg')
                      
        # Set parameters for synplot
        synplot_params = kwargs.copy()
        synplot_params.update(syn_params)
            
        # Synthesize spectrum    
        syn = s4.synthesis.Synplot(teff, logg, synplot_path, idl,  **synplot_params)

        syn.run()

        #import matplotlib.pyplot as plt
        #plt.plot(syn.spectrum[:,0], syn.spectrum[:,1], label='synthetic')
        #plt.plot(kwargs['observ'][:,0], kwargs['observ'][:,1], label='observed')
        #plt.legend(loc='best')
        #plt.show()
        
        #if not noplot:
        #    figname = plot_title.pdf
        #    syn.plot(title=plot_title, file_name=figname)
            #plt.show()
            

        # Apply scale and radial velocity if needed
        if 'scale' in syn.parameters:
            syn.apply_scale()
            
        if 'rv' in syn.parameters:
            syn.observation[:, 0] *= s4.spectools.rvcorr(syn.parameters['rv'])
            
        #Do an interpolation

        flm = np.interp(syn.observation[:,0], syn.spectrum[:, 0], 
                        syn.spectrum[:, 1])#/max(syn.observation[:,1])
        
        #Some kind of normalization on the observed flux?
        fobm = syn.observation[:,1]#/max(syn.observation[:,1])
         
        # Calculate the chi**2
        
        chisq = sum((fobm - flm)**2/flm)# * weights)
        #chisq = chisq * max(fobs)                 #????            

        # store the values of the parameters
        if n_params == 1: 
            store_values[fit_params.keys()[0]][n] = it
        else:    
            for j, v in enumerate(fit_params):
                store_values[v][n] = it[j]
        store_values['chisquare'][n] = chisq

    # print the best result
    print 'Best fit:'
    best_fit_vector = store_values[np.argmin(store_values['chisquare'])]
    best_fit_string = ''
    for key in fit_params.keys():
        best_fit_string += '{} = {}, '.format(key, best_fit_vector[key][0])
    best_fit_string += 'chi^2 = {:.6f}.'.format(best_fit_vector['chisquare'][0])    
    print best_fit_string
