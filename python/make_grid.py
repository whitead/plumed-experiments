class Grid:
    """PLUMED grid class
    
    Member variables: 
      nbins: number of bins for each cv
      min: min for each cv
      max: max for each cv
      period: logical, indicating periodicity for each cv
      ncv: number of cv
      dx: grid spacing for each cv
      pot: array for potential      
      types: type for each cv
    """
    
    def _print_header_array(self, name, array, output):
        output.write('#! {} '.write(name))
        for a in array:
            output.write('{} '.write(a))
        output.write('\n')
                     

    def write_header(self, output):
        output.write('#! FORCE 0\n')
        self._print_header_array('NVAR', self.ncv, output)
        self._print_header_array('TYPE', self.types, output)
        self._print_header_array('BIN', self.nbins, output)
        self._print_header_array('MIN', self.min, output)
        self._print_header_array('MAX', self.max, output)
        self._print_header_array('PBC', self.period, output)
            
            
            
