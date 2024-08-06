"""RMG IO plugin module."""

from aseconv.pluginbase import AsecIO

class RMGIO(AsecIO):
    """RMG IO Plugin. Currently supports ``write`` only.
    """
    
    def __init__(self):
        pass

    def infos(s):
        return {"help": "RMG file", "typeexts": {"rmg": ".rmg"}}

    def read(s, file, type, **kwargs):
        raise Exception("RMG reading is not supported yet")

    def write(self, args, atom, type, ofile):
        # lstr=['#******** REAL SPACE GRID ********' ]
        # wgrid=(((atom.cell.lengths()/0.2)+1)/2).round()*2
        # lstr.append(f'wavefunction_grid="{vec2str(wgrid,False,"{:.0f}")}"')
        # lstr.append('potential_grid_refinement="2"')
        lstr = []
        lstr.extend(
            [
                "\n#****  LATTICE and ATOMS  ****",
                'bravais_lattice_type = "None"',
                'crds_units           = "Angstrom"',
                'lattice_units        = "Angstrom"',
            ]
        )

        if args.C:
            coord = "Absolute"
            pos = atom.get_positions()
        else:
            coord = "Cell Relative"
            pos = atom.get_scaled_positions(wrap=False)

        lstr.append(f'atomic_coordinate_type = "{coord}"')
        lstr.append('lattice_vector = "')
        for i in atom.cell:
            lstr.append(self.vec2str(i))

        lstr.append('"\natoms = "')  # space must exist before '='
        sym = atom.get_chemical_symbols()
        for i, xyz in enumerate(pos):
            lstr.append(sym[i] + "   " + self.vec2str(xyz))
        lstr.extend('"\n')
        with open(ofile, "wt") as f:
            f.write("\n".join(lstr))
            f.write("\n")
