"""LAMMPS IO plugin module."""

from aseconv.pluginbase import AsecIO

class LammpsIO(AsecIO):
    """Lammps IO plugin.

    Supports atomic(lmp), charge(lmpc), and full(lmpf) type.
    Can read ``Atom Type Labels``
    
    """

    def __init__(s):
        pass

    def infos(s):
        return {
            "help": "LAMMPS file",
            "typeexts": {"lmp": ".lmp", "lmpc": ".lmpc", "lmpf": ".lmpf"},
        }

    def __read_fix(s, args):
        if args.x is None:
            return []
        fixes = glob.glob(args.x)
        arrcon = []
        for i in fixes:
            pf = Path(i)
            fstr = re.sub("[^_]*_", "", pf.stem)
            fix = [False] * 3
            for lst in fstr:
                for p in lst[0]:
                    idx = ord(p.upper()) - ord("X")
                    fix[idx] = True
            print("  - Applying {}, fix {} ".format(i, fstr))
            with open(pf) as f:
                lines = f.readlines()
            alst = []
            lstart = 0
            for ln, l in enumerate(lines, 1):
                if l.startswith("ITEM: ATOMS id"):
                    lstart = ln
                if lstart == 0 or ln <= lstart:
                    continue
                if len(l.strip()) == 0:
                    continue
                alst.append(int(l.split()[0]) - 1)
            arrcon.append([fix, alst])
        # print(arrcon)
        return arrcon

    def _lammps_read(args, pfile):
        from ase.data import (
            atomic_numbers,
        )  # , atomic_names, atomic_masses, covalent_radii

        with open(pfile) as f:
            lines = f.readlines()
            adict = {}
            for il, l in enumerate(lines):
                if "atom types" in l:
                    ntype = int(l.split()[0])
                    continue
                if l.startswith("Atom Type Labels"):
                    for i in range(il + 2, il + 2 + ntype):
                        arr = lines[i].split()
                        adict.update({int(arr[0]): atomic_numbers[arr[1]]})
                    return adict
        return None

    def read(s, file, type, **kwargs):
        format = "lammps-data"
        pfile=Path(file)
        ext=pfile.suffix
        if ext == ".lmp":
            parm = {"style": "atomic"}
        elif ext == ".lmpc":
            parm = {"style": "charge"}
        elif ext == ".lmpf":
            parm = {"style": "full"}
        zdict = _lammps_read(args, pfile)
        if zdict is not None:
            parm.update({"Z_of_type": zdict})
        return ase.io.read(pfile, **kwargs, **parm)

    def write(self, args, atom, type, ofile):
        # global args.Gslabidx
        import numpy as np
        import math

        syms = atom.get_chemical_symbols()
        amass = atom.get_masses()
        pos = atom.get_positions()
        types = []
        masses = []
        lmp_atoms = []
        cstr = ""
        if type == "lmpc":
            cstr = " 0.000"
        for i, s in enumerate(syms):
            if not s in types:
                types.append(s)
                masses.append(amass[i])
            r = pos[i]
            mid = ""  # mid=1
            idx = types.index(s) + 1
            rz = r[2]
            # if (args.Gslabidx > 0): rz=0 # why?
            lmp_atoms.append(
                "{:5} {} {:2}{} {}".format(i + 1, mid, idx, cstr, self.vec2str(r))
            )  # r[0],r[1],rz))

        lstr = []
        hi = atom.cell.lengths()

        lstr.append("#TYPEMAP:%s" % " ".join(types))
        lstr.append("\n{:6} atoms\n{:6} atom types".format(len(syms), len(types)))

        cell = atom.get_cell()
        # vacuum=args.v
        # if (vacuum == None): vacuum=0
        rvac = 0
        if cell.volume != 0:
            lo = np.array([0] * 3)
            cp = cell.cellpar(radians=True)
            a = cp[0]
            b = cp[1]
            c = cp[2]
            A = cp[3]
            B = cp[4]
            C = cp[5]
            lx = a
            xy = b * math.cos(C)
            xz = c * math.cos(B)
            ly = (b**2 - xy**2) ** 0.5
            yz = (b * c * math.cos(A) - xy * xz) / ly
            lz = (c**2 - xz**2 - yz**2) ** 0.5
            hi = [lx, ly, lz]
        else:
            lo = pos.min(axis=0)
            hi = pos.max(axis=0)
            # rvac=vacuum
            xy = 0
            xz = 0
            yz = 0

        avg = pos.mean(axis=0)
        for i in range(3):
            axis = chr(ord("x") + i)
            alo = lo[i]
            ahi = hi[i]
            lstr.append(
                "{0} {1}lo {1}hi".format(
                    self.vec2str([alo - rvac / 2, ahi + rvac / 2]), axis
                    )
            )
        # xy=bcosC, xz=c cos B, yz=(b*c*cosA-xy*xz)/ly

        lstr.append("{} xy xz yz".format(self.vec2str([xy, xz, yz])))

        lstr.append("\nAtom Type Labels\n")
        for i, at in enumerate(types, start=1):
            lstr.append(f" {i:2} {at:3}")

        lstr.append("\nMasses\n")  ## Blank lines are critical
        for i, at in enumerate(masses, 1):
            lstr.append(" {:2} {:8.3f}  # {}".format(i, at, types[i - 1]))

        lstr.append("\nAtoms\n")  ## Blank lines are critical
        lstr.extend(lmp_atoms)

        with open(ofile, "wt") as f:
            f.write("\n".join(lstr))
            f.write("\n")
