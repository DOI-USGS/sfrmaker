"""
Local tests
"""
import sys
sys.path.append('/Users/aleaf/Documents/GitHub/SFR/')
from preproc import lines

def test_haskell():
    flowlinesshp = 'data/WI_flowlines.shp'
    domainshp = 'data/mfbounds.shp'
    model_grid_shp = 'data/mf.shp'

    lns = lines(flowlinesshp, mf_grid=model_grid_shp, model_domain=domainshp)
    lns.to_sfr()

if __name__ == '__main__':
    test_haskell()