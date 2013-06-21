'''
gdlwrapper is a module containing a funtion to run GDL.

GDL is the open-source alternative to IDl.

This is a workaround. The best way would be to translate all GDL/IDL scripts  
to python.
'''
#=============================================================================
# Modules
from subprocess import Popen, PIPE
#=============================================================================

#=============================================================================
# Run gdl
def run_gdl(inp, do_log = False):
    """Run gdl"""
    
    if do_log:
        with open('gdl.log', 'w') as log:
            gdl = Popen(['nice', '-n0', 'gdl'], stdin = PIPE, \
                           stdout = log, stderr = log)
            gdl.communicate(inp)
                           
    else:
        print '########### gdl ##################'
        gdl = Popen(['nice', '-n0', 'gdl'], stdin = PIPE)
        gdl.communicate(inp)
        print '######## Quitting gdl ############'    
#=============================================================================
