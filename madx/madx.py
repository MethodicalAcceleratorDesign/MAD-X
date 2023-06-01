import cpymad.madx
import re
import numpy as np
import scipy.optimize
from collections import ChainMap

from .tablemixin import TableMixIn


class Table(cpymad.madx.Table, TableMixIn):
    pass


class TableMap(cpymad.madx.TableMap):
    def __getitem__(self, name):
        try:
            t = Table(name, self._libmadx)
            t.set_index("name")
            return t
        except ValueError:
            raise KeyError("Table not found {!r}".format(name)) from None


class Madx(cpymad.madx.Madx):
    def __init__(
        self,
        libmadx=None,
        command_log=None,
        stdout=None,
        history=None,
        prompt=None,
        **Popen_args,
    ):
        super().__init__(libmadx, command_log, stdout, history, prompt, **Popen_args)

        self.table = TableMap(self._libmadx)

    def track_single_4d(mad, coord=np.zeros(4)):
        mad.track(onepass=True)
        mad.start(
            x=coord[0], px=coord[1], y=coord[2], py=coord[3])
        mad.run()
        mad.endtrack()
        t = mad.table["track.obs0001.p0001"]
        out = [t[n][-1] for n in ["x", "px", "y", "py"]]
        return np.array(out)

    def track_single_6d(mad, coord=np.zeros(6)):
        mad.track(onepass=True)
        mad.start(
            x=coord[0], px=coord[1], y=coord[2], py=coord[3], t=coord[4], pt=coord[5])
        mad.run()
        mad.endtrack()
        t = mad.table["track.obs0001.p0001"]
        out = [t[n][-1] for n in ["x", "px", "y", "py", "t", "pt"]]
        return np.array(out)


    def find_closed_orbit_track_4d(mad, x0=None, xtol=1e-10,mode='4d'):
        def ftosolve(x):
            if abs(x).max()>1:
                raise ValueError(f"{x} too large")
            return mad.track_single(x) - x

        res = scipy.optimize.fsolve(ftosolve, x0)
        return res

    def track_single_6d(mad, coord):
        mad.track(onepass=True)
        mad.start(
            x=coord[0], px=coord[1], y=coord[2], py=coord[3], t=coord[4], pt=coord[5]
        )
        mad.run()
        mad.endtrack()
        t = mad.table["track.obs0001.p0001"]
        out = [t[n][-1] for n in ["x", "px", "y", "py", "t", "pt"]]
        return np.array(out)

    def find_closed_orbit_track_4d(mad, x0=np.zeros(6), xtol=1e-10):
        def ftosolve(x):
            if abs(x).max()>1:
                raise ValueError(f"{x} too large")
            return mad.track_single(x) - x

        res = scipy.optimize.fsolve(ftosolve, x0)
        return res

    def find_closed_orbit_twiss(mad):
        t = mad.twiss()
        return np.array([t[n][0] for n in ["x", "px", "y", "py", "t", "pt"]])
