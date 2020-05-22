"""
Module for creating RIV package input
"""
from sfrmaker.base import DataPackage


class RivData(DataPackage):
    """

    Parameters
    ----------
    stress_period_data : DataFrame
        DataFrame with input information for the RIV package stress period data
        (stress_period_data in FloPy).
        Columns:
        node :
    grid :
    model :
    model_length_units :
    model_time_units :
    package_name :
    kwargs :
    """
    def __init__(self, stress_period_data=None, grid=None,
                 model=None,
                 model_length_units="undefined", model_time_units='d',
                 package_name=None,
                 **kwargs):
        DataPackage.__init__(self, grid=grid, model=model,
                         model_length_units=model_length_units,
                         model_time_units=model_time_units,
                         package_name=package_name)

        self.stress_period_data = stress_period_data

    @classmethod
    def from_lines(cls, lines, grid=None,
                   active_area=None, isfr=None,
                   model=None,
                   model_length_units='undefined',
                   minimum_reach_length=None,
                   cull_flowlines_to_active_area=True,
                   consolidate_conductance=False, one_reach_per_cell=False,
                   model_name=None,
                   **kwargs):
        """
        Create an instance of Riv from an SFRmaker.lines object.


        Parameters
        ----------
        lines : SFRmaker.lines instance :
        grid :
        active_area :
        isfr :
        model :
        model_length_units :
        minimum_reach_length :
        cull_flowlines_to_active_area :
        consolidate_conductance :
        one_reach_per_cell :
        model_name :
        kwargs :

        Returns
        -------
        riv : SFRmaker.RivData instance
        """
        raise NotImplementedError("from_lines not implemented yet.")

