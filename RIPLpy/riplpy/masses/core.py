

# ========================

# RIPLpy
import riplpy.collections as _c
import riplpy.db as _db

# ========================

__all__ = ('MassDatabase', )

class MassDatabase(_db.NuclideDatabase):

    def mass_excess(self, n: _c.Nuclide) -> float:
        """Return the mass excess of the requested nucleus (n). """
        return self.data[n]['Mth']
