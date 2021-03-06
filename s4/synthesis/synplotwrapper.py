#=============================================================================
# Modules
import numpy as np
import os
import matplotlib.pyplot as plt
import re
import warnings
import lineid_plot
from ..spectools import rvcorr
from ..utils import File
from ..plottools import plot_windows
from ..io import specio, wrappers
#=============================================================================



class Synplot:
    """Add docstring"""

    def __init__(self, teff, logg, synplot_path = None, idl = True,
                 **kwargs):
        if synplot_path is None:
            self.spath = os.getenv('HOME')+'/.s4/synthesis/synplot/'
        else:
            self.spath = synplot_path

        # Set software to run Synplot.pro
        if isinstance(idl, str):
            self.software = idl
        elif idl:
            #self.software = '/Applications/itt/idl/bin/idl'  #Only works for me now!
            self.software = 'idl'
        else:
            self.software = 'gdl'

        # Setting teff and logg on the dictionary
        kwargs['teff'] = teff
        kwargs['logg'] = logg

        #Check if some params were defined
        if 'wstart' not in kwargs.keys():
            kwargs['wstart'] = float(File(self.spath + 'fort.19').\
                               head()[0].split()[0]) * 10
            print 'wstart not defined.'
            print 'Setting as {:.2f} Angstrons.\n'.format(kwargs['wstart'])

        if 'wend' not in kwargs.keys():
            kwargs['wend'] = float(File(self.spath + 'fort.19').\
                               tail()[0].split()[0]) * 10
            print 'wend not defined.'
            print 'Setting as {:.2f} Angstrons.\n'.format(kwargs['wend'])

        self.parameters = kwargs

        # Check if a observation spectrum is available
        if 'observ' in self.parameters:
            if isinstance(self.parameters['observ'], str):
              self.observation = specio.load_spectrum(self.parameters['observ'])
            elif isinstance(self.parameters['observ'], np.ndarray):
              self.observation = self.parameters['observ']
            #Delete entry to not input in IDL
            del self.parameters['observ']

        #Override IDL plotting
        self.parameters['noplot'] = '1'

        # check for normalization
        if 'relative' not in self.parameters:
            self.parameters['relative'] = 0

        # The idl rv keyword doesn't seem to work? Do it here instead
        self.rv = 0.0
        if 'rv' in self.parameters:
            self.rv = self.parameters.pop('rv')

    #=========================================================================
    #
    def synplot_input(self):
        """Build the synplot command to IDL/GDL."""

        synplot_command = [key+' = '+str(value)                              \
                           for key, value in self.parameters.iteritems()]

        cmd = "CD, '"+self.spath+"' & synplot, "+ \
                        ', '.join(synplot_command)

        return self.software + ' -e "' + cmd + '"'
    #=========================================================================

    #=========================================================================
    # Run synplot and return the computed spectra
    def run(self):
        """Run synplot and store the computed spectra"""

        # remove old calculated spectrum
        # Stack Overflow #10840533
        try:
            os.remove((self.spath + 'fort.11'))
        except OSError:
            pass

        wrappers.run_command(self.synplot_input(), do_log = True)

        #load synthetized spectra
        try:
            self.spectrum = np.loadtxt(self.spath + 'fort.11')
            self.spectrum[:,0] *= (1.0 + self.rv/3e5)
        except IOError:
            raise IOError('Calculated spectrum is not available. Check if ' +
                'syn(spec|plot) ran correctly.')

    def savetxt(self, file_name, *args):
        """
        Save spectrum fo a file.

        Parameters
        ----------

        file_name: str;
            Name of the file to be saved.

        args:
            Numpy.savetxt arguments.
        """

        self.check_if_run()

        np.savetxt(file_name, self.spectrum, *args)

    # Plot
    def plot(self, ymin = None, ymax = None, windows = None, file_name = None,
             title = None, ident = False):
        """
        Plot the synthetic spectra.
        If the synthetic spectra were not calculated, it will calculate.

        Parameters
        ----------

        ymin : lower limit on y-axis
        ymax : upper limit on y-axis

        file_name: str;
            Name of the file to be saved.
        """

        self.check_if_run()

        # make a copy of array
        spectrum_copy = self.spectrum.copy()
        if hasattr(self, 'observation'):
            observation_copy = self.observation.copy()
            #Apply radial velocity correction if needed.
            if 'rv' in self.parameters:
                observation_copy[:, 0] *= rvcorr(self.parameters['rv'])

        # Apply scale correction needed
        if 'scale' in self.parameters:
            spectrum_copy[:, 1] *= self.parameters['scale']

        # check if figure was already plotted
        if plt.fignum_exists(1):
            fig_exists = True
            plt.clf()
        else:
            fig_exists = False
        # Plot
        fig = plt.figure(num = 1)

        # Set axes.
        # If identification of the linesis requiredm set different size to axes
        if ident:
            ax = fig.add_axes([0.1, 0.1, 0.85, 0.6])
        else:
            ax = fig.gca()

        # If a observation spectra is available, plot it
        if hasattr(self, 'observation'):
            ax.plot(observation_copy[:, 0], observation_copy[:, 1], c='#377eb8',
                    label = 'Observation')

        # Plot synthetuc spectrum
        ax.plot(spectrum_copy[:, 0], spectrum_copy[:, 1], c='#e41a1c',
                label = 'Synthetic')

        # If windows were set, plot it
        if windows is not None:
            plot_windows(windows)

        # set labels
        plt.xlabel(r'Wavelength $(\AA)$')
        if self.parameters['relative'] != 0:
            if ymin is None:
                ymin = 0
            if ymax is None:
                ymax = 1.05
            plt.ylabel('Normalized Flux')
        else:
            plt.ylabel('Flux')
        # Set size of plot
        plt.xlim([self.parameters['wstart'], self.parameters['wend']])
        if ymin is not None:
            plt.ylim(ymin = ymin)
        if ymax is not None:
            plt.ylim(ymax = ymax)
        ####


        # Identify lines, if required
        if ident:
            # Obtain the spectral line wavelength and identification
            line_wave, line_label = self.lineid_select(ident)
            lineid_plot.plot_line_ids(spectrum_copy[:, 0],
                                      spectrum_copy[:, 1],
                                      line_wave, line_label, label1_size = 10,
                                      extend = False, ax = ax,
                                      box_axes_space = 0.15)

        plt.legend(fancybox = True, loc = 'lower right')

        if title is not None:
            plt.title(title, verticalalignment = 'baseline')

        # Plot figure
        if not fig_exists:
            fig.show()
        else:
            fig.canvas.draw()

        # Save file
        if file_name is not None:
            plt.savefig(file_name, dpi = 100)

    #=========================================================================

    #=========================================================================
    # Select lines to line identification
    def lineid_select(self, ident):
        """Identify lines to be plot by lineid_plot"""
        # List of chemical elements
        table = open(self.spath + 'fort.12').read() +                         \
                open(self.spath + 'fort.14').read()

        # Pattern for regex
        ptrn = r'(\d{4}\.\d{3})\s+(\w{1,2}\s+I*V*I*).+(\b\d+\.\d)\s+(\*+)'

        #Find patterns
        regex_table = re.findall(ptrn, table)

        # Parse table
        wavelengths = [float(line[0]) for line in regex_table
                       if float(line[2]) >= ident]
        chem_elements = [line[1] + ' ' + line[0] + '  ' + line[2]
                         for line in regex_table
                         if float(line[2]) >= ident]

        return wavelengths, chem_elements
    #=========================================================================

    #=========================================================================
    #Apply scale
    def apply_scale(self):
        """ Apply scale. """
        self.spectrum[:, 1] *= self.parameters['scale']


    def check_if_run(self):
        """
        Check if spectrum was already calculated. If not, calculate it.
        """

        if not hasattr(self, 'spectrum'):
            self.run()
