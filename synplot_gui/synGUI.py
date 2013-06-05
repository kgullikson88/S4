#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Create and interactive synplot.

I am following the tutorial found at http://zetcode.com/tutorials/pyqt4/
and example from here
http://eli.thegreenplace.net/2009/01/20/matplotlib-with-pyqt-guis/
"""
#==============================================================================

# Import Modules
import os
import sys
from PyQt4 import QtGui, QtCore
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
import json
import s4


#==============================================================================

#==============================================================================
# Global variables
FRAME_WIDTH = 1020
FRAME_HEIGHT = 480

HOME_PATH = os.getenv('HOME') + '/'
CONFIG_FILE = HOME_PATH + '.s4_config.json'


#==============================================================================

#==============================================================================
# 
class Widget(QtGui.QWidget):
    
    def __init__(self):
        super(Widget, self).__init__()

        # set font for tips
        QtGui.QToolTip.setFont(QtGui.QFont('SansSerif', 10))

        # Frame setup
        self.resize(FRAME_WIDTH, FRAME_HEIGHT)
        self.center()
        self.setWindowTitle('synGUI')

        self.create_frame()
       
        # load config
        try:
            self.load_config()
        except IOError:
            # set up a basic configuration
            self.teff = "20000"
            self.logg = "4"
            self.parameters = dict(wstart = "4460", wend = "4480", rv = "0",
                                   vrot = "0", vturb = "0", vmac_rt = "0", 
                                   relative = "1", scale = "1")
        
        # set text boxes              
        self.teff_textbox.setText(self.teff) 
        self.logg_textbox.setText(self.logg) 
        self.wstart_textbox.setText(self.parameters['wstart']) 
        self.wend_textbox.setText(self.parameters['wend']) 
        self.rv_textbox.setText(self.parameters['rv'])        
        self.vrot_textbox.setText(self.parameters['vrot'])
        self.vturb_textbox.setText(self.parameters['vturb'])
        self.vmac_textbox.setText(self.parameters['vmac_rt'])
        self.scale_textbox.setText(self.parameters['scale'])        
        if 'observ' in  self.parameters: 
            self.obs_textbox.setText(self.parameters['observ']) 
        
        # set check boxes status
        if self.parameters["relative"] == "0":
            self.norm_cb.setChecked(False)
        else:                                             
            self.norm_cb.setChecked(True)       
        
        self.synplot()
        
    # About    
    def on_about(self):
        msg = """ A GUI for synplot
        
         * describe using ReST
        """
        QtGui.QMessageBox.about(self, "About synGUI", msg.strip())        

    #=========================================================================
    # Core modules
    # synplotmodule    
    def synplot(self):
        global teff, logg
        # Get and update parameters
        self.teff = self.teff_textbox.text()
        self.logg = float(self.logg_textbox.text())
        self.parameters['wstart'] = self.wstart_textbox.text()
        self.parameters['wend'] = self.wend_textbox.text()
        self.parameters['vrot'] = self.vrot_textbox.text()
        self.parameters['vturb'] = self.vturb_textbox.text()
        self.parameters['vmac_rt'] = self.vmac_textbox.text()
        self.parameters['rv'] = float(self.rv_textbox.text())
        self.parameters['scale'] = float(self.scale_textbox.text())
        if self.obs_textbox.text() != '':
            self.parameters['observ'] = str(self.obs_textbox.text())
        else:
            self.parameters.pop('observ', None)
                
        if self.norm_cb.isChecked():
            self.parameters['relative'] = "1"
        else:
            self.parameters['relative'] = "0"     
        # run synplot
        self.syn = s4.synthesis.Synplot(self.teff, self.logg, 
                                        **self.parameters)
        self.syn.run()
        # make corrections on synthsized spectrum
        self.syn.apply_rvcorr()
        self.syn.apply_scale()
        # draw    
        self.on_draw()
        # save configuration file
        self.save_config()

    #=========================================================================
    
    #=========================================================================
    # GUI             

    def create_frame(self):

        self.main_frame = QtGui.QWidget()
        
        # Create the mpl Figure and FigCanvas objects. 
        # 5x4 inches, 100 dots-per-inch
        #
        self.dpi = 100
        self.fig = Figure((5.0, 4.0), dpi = self.dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.main_frame)  
        
        # Since we have only one plot, we can use add_axes 
        # instead of add_subplot, but then the subplot
        # configuration tool in the navigation toolbar wouldn't
        # work.
        #
        self.axes = self.fig.add_subplot(111)   
        
        # Create the navigation toolbar, tied to the canvas
        #
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame) 
        
        # Other GUI controls
        #
        # teff
        self.teff_label = QtGui.QLabel('teff')
        self.teff_textbox = self.add_text_input('Effective Temperature')
        # logg
        self.logg_label = QtGui.QLabel('logg')
        self.logg_textbox = self.add_text_input('Logarithm of surface ' + \
                                                'gravity')
        # wstart
        self.wstart_label = QtGui.QLabel('wstart')
        self.wstart_textbox = self.add_text_input('Starting Wavelength')
        # wend
        self.wend_label = QtGui.QLabel('wend')
        self.wend_textbox = self.add_text_input('Ending Wavelength')
        # radial velocity
        self.rv_label = QtGui.QLabel('rv')
        self.rv_textbox = self.add_text_input('Radial velocity')        
        # vseni
        self.vrot_label = QtGui.QLabel('vrot')
        self.vrot_textbox = self.add_text_input('Projected rotational ' +\
                                                 'velocity')
        # microturbulent velocity
        self.vturb_label = QtGui.QLabel('vturb')
        self.vturb_textbox = self.add_text_input('Microturbulent velocity')
        # macroturbulent velocity
        self.vmac_label = QtGui.QLabel('vmac_RT')
        self.vmac_textbox = self.add_text_input('Radial-tangential ' + \
                                                 'macroturbulent velocity')
                                                 
        # scale
        self.scale_label = QtGui.QLabel('Scale')
        self.scale_textbox = self.add_text_input('Scale factor')                                                 
        # normalization
        self.norm_cb = QtGui.QCheckBox("Normalization")
        self.norm_cb.setToolTip("If checked, normalize spectrum.")
        
        # observation file
        # label
        self.obs_label = QtGui.QLabel("Observation file")
        # text edit for observation file path
        self.obs_textbox = self.add_text_input()
        self.obs_textbox.setMaximumWidth(500)
        # button to open file dialog
        self.obs_button = QtGui.QPushButton('Open', self)
        self.obs_button.clicked.connect(self.obs_file_dialog)
        self.obs_button.setToolTip('Open observation file')
        #self.obs_button.resize(self.obs_button.sizeHint())
        self.obs_button.setMaximumWidth(60)
 
        # button to run synplot
        self.run_button = QtGui.QPushButton('Run', self)
        self.run_button.clicked.connect(self.synplot)
        self.run_button.setToolTip('Press to run <b>synplot</b>')
        #self.run_button.resize(self.run_button.sizeHint())
        self.run_button.setMaximumWidth(50)
        

        #
        # Layout with box sizers
        #
        # define grid
        grid = QtGui.QGridLayout()
        #grid.setSpacing(10)
        
        #set matplotlib canvas
        grid.addWidget(self.canvas, 0, 0, 11, 1)
        grid.addWidget(self.mpl_toolbar, 12, 0, 1, 1)

        # Define first row
        grid.addWidget(self.teff_label, 0, 1)
        grid.addWidget(self.teff_textbox, 0, 2)
        grid.addWidget(self.logg_label, 0, 3)
        grid.addWidget(self.logg_textbox, 0, 4)
        grid.addWidget(self.wstart_label, 0, 5)
        grid.addWidget(self.wstart_textbox, 0, 6)
        grid.addWidget(self.wend_label, 0, 7)
        grid.addWidget(self.wend_textbox, 0, 8)          
        # Define second row  
        grid.addWidget(self.rv_label, 1, 1)
        grid.addWidget(self.rv_textbox, 1, 2)        
        grid.addWidget(self.vrot_label, 1, 3)
        grid.addWidget(self.vrot_textbox, 1, 4)   
        grid.addWidget(self.vturb_label, 1, 5)
        grid.addWidget(self.vturb_textbox, 1, 6)
        grid.addWidget(self.vmac_label, 1, 7)
        grid.addWidget(self.vmac_textbox, 1, 8)
        # Define third row
        grid.addWidget(self.scale_label, 2, 1)
        grid.addWidget(self.scale_textbox, 2, 2)        
        grid.addWidget(self.norm_cb, 2, 3, 1, 3)
        # Define fourth row
        grid.addWidget(self.obs_label, 3, 1, 1, 3)        
        # Define fifth row
        grid.addWidget(self.obs_textbox, 4, 1, 1, 7)
        grid.addWidget(self.obs_button, 4, 8)
        # Define seventh row        
        grid.addWidget(self.run_button, 12, 8)
        # set grid  
        self.setLayout(grid) 
        
        #self.setCentralWidget(self.main_frame)
       
        self.show()

    # Draw canvas        
    def on_draw(self):
        """ Redraws the figure
        """
        # clear the axes and redraw the plot anew
        #
        self.axes.clear()        
        
        self.axes.plot(
            self.syn.spectrum[:, 0], 
            self.syn.spectrum[:, 1], label = "Synthetic")
        # plot observation, if availabe
        if 'observ' in self.parameters:
            self.axes.plot(
                self.syn.observation[:, 0],     
                self.syn.observation[:, 1], label = "Observation")                
            self.axes.legend(fancybox = True, loc = 'lower right') 
        
        # set x dimension
        xmin = float(self.parameters['wstart'])
        xmax = float(self.parameters['wend'])
        self.axes.set_xlim(xmin, xmax)
        
        self.canvas.draw()        
        
    # Set window to center of desktop
    def center(self):
        
        qr = self.frameGeometry()
        cp = QtGui.QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())        
    #=========================================================================
        
    #=========================================================================
    # Support Modules
    
    # add Label + text input
    def add_text_input(self, tip = None):
        text_input = QtGui.QLineEdit()
        if tip is not None:
            text_input.setToolTip(tip)
        text_input.setMaximumWidth(55)     
        return text_input
        
    def load_config(self):
        """Load the configuration file"""
        self.parameters = json.load(open(CONFIG_FILE))
        self.teff = self.parameters.pop('teff')
        self.logg = self.parameters.pop('logg')
        
    def save_config(self):
        """Save configuration file"""
        with open(CONFIG_FILE, 'w') as f:
            spam = {key : str(value) for key, value in self.parameters.iteritems()}
            spam.update({'teff' : str(self.teff), 'logg' : str(self.logg)})
            json.dump(spam, f, sort_keys = True, indent = 4, 
                      separators=(',', ':'))        
        
    def obs_file_dialog(self):
        """Open a file dialog to open observation file"""
        #self.fileDialog = QtGui.QFileDialog(self)
        #self.fileDialog.show()
        fname = str(QtGui.QFileDialog.getOpenFileName(self,"Open File"))
        self.obs_textbox.setText(fname)
        
                 
                 
    """    
    # ????
    def add_actions(self, target, actions):
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)
        
    # Create (generic) action
    def create_action(self, text, slot = None, shortcut = None, icon = None, 
                      tip = None, checkable = False, signal = "triggered()"):
        action = QtGui.QAction(text, self)
        if icon is not None:
            action.setIcon(QtGui.QIcon(":/%s.png" % icon))
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            self.connect(action, QtCore.SIGNAL(signal), slot)
        if checkable:
            action.setCheckable(True)
        return action
    
    """    
    def test(self):
        print 'this is a test' 
    #=========================================================================    
        
#==============================================================================

#==============================================================================
#
def main():
    
    app = QtGui.QApplication(sys.argv)
    wdg = Widget()
    sys.exit(app.exec_()) 
#==============================================================================

#==============================================================================
# 
if __name__ == '__main__':
    main()
#==============================================================================
