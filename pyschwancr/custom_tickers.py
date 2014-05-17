
import numpy as np
from matplotlib.ticker import ScalarFormatter, Formatter

class SchwancrFormatter(ScalarFormatter):

   def __init__(self,useOffset=True, useMathText=False, precision=None):
      ScalarFormatter.__init__(self,useOffset=useOffset,useMathText=useMathText )
      self.precision = precision

   def _formatSciNotation(self, s):
      # transform 1e+004 into 1e4, for example
      tup = s.split('e')
      try:
         significand = tup[0].rstrip('0').rstrip('.')
         if self.precision:
            significand = '%%.%df'%(self.precision,) %(float(significand),)
         sign = tup[1][0].replace('+', '')
         exponent = tup[1][1:].lstrip('0')
         if self._useMathText or self._usetex:
            if float(significand) == 1:# '1':
               # reformat 1x10^y as 10^y
               significand = ''
            if exponent:
               exponent = '10^{%s%s}'%(sign, exponent)
            if significand and exponent:
               return r'%s{\times}%s'%(significand, exponent)
            else:
               return r'%s%s'%(significand, exponent)
         else:
            s = ('%se%s%s' %(significand, sign, exponent)).rstrip('e')
            return s
      except IndexError, msg:
         return s
   def _set_format(self, vmin, vmax):
      # set the format string to format all the ticklabels
      # The floating point black magic (adding 1e-15 and formatting
      # to 8 digits) may warrant review and cleanup.
      locs = (np.asarray(self.locs)-self.offset) / 10**self.orderOfMagnitude+1e-15
      if self.precision:
         self.format = '%%.%df'%self.precision
      else:
         sigfigs = [len(str('%1.8f'% loc).split('.')[1].rstrip('0')) \
                    for loc in locs]
         sigfigs.sort()
         self.format = '%1.' + str(sigfigs[-1]) + 'f'
      if self._usetex:
          self.format = '$%s$' % self.format
      elif self._useMathText:
          self.format = '$\mathdefault{%s}$' % self.format

