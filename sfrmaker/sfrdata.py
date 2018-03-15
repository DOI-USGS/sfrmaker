import sys
sys.path.append('/Users/aleaf/Documents/GitHub/flopy3')
import numpy as np
import pandas as pd
import flopy


class sfrdata:

    def __init__(self):

        pass

    @staticmethod
    def from_lines(lines, grid):
        """Create a streamflow routing package dataset from
        a lines object and model grid object.

        Parameters
        ----------
        lines : instance of sfrmaker.lines
        grid : instance of sfrmaker.grid
        """

        if grid.active_area is not None:
            lines.cull(grid.active_area)

        reach_data = lines.intersect(grid)

        j=2


