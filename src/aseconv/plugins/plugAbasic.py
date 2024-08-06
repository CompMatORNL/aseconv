""" Default geometry plugins.
        
"""
import re
import numpy as np
import ase
import aseconv.geoutil as gu
from aseconv.pluginbase import AsecPlug

def _select_by_xyz(atom, str):
    pos = atom.positions
    x = pos[:, 0]
    y = pos[:, 1]
    z = pos[:, 2]
    arr = eval(f"({str})")
    if arr.sum() == 0:
        raise Exception(f" - No atom is selected by '{sum}' ...")
    catom = eval(f"atom[arr]")
    # print(np.where(arr))
    return catom, np.where(arr)[0]


class APlugNoConst(AsecPlug):
    """Remove all the constraints."""
    
    def __init__(self, asec):
        super().__init__()
        asec.add_argument(self, "--noc", action="store_true", help="Remove all the onstraints.")

    def output_postfix(self, args, opt):
        return f"_NoC"

    def process(self, args, atom, val):
        atom.set_constraint(None)
        return atom


class APlugConst(AsecPlug):
    """Add constraints.

    """

    def __init__(self, asec):
        super().__init__()
        asec.add_argument(
            self,
            "--con",
            metavar='"expr of {x-z}"',
            type=str,
            help="Add additional constraint atoms of the boolean expression. Ex) `(z>0)&(z<10)`.",
        )

    def output_postfix(self, args, opt):
        return f"_C{self.safe_name(opt)}"

    def process(self, args, atom, val):
        from ase.constraints import FixAtoms, FixScaled, FixCartesian
        from collections import OrderedDict

        carr = []
        fixinds = set([])
        fixscaled = {}

        for lst in atom.constraints:
            if isinstance(lst, FixCartesian):  # FHI-aims partial constraint
                fixscaled.update({lst.a: ~np.array(lst.mask)})
            elif isinstance(lst, FixScaled):  # VASP constraint
                ldic = lst.todict()
                fixscaled.update({lst.a: ldic["kwargs"]["mask"]})
            else:
                fixinds = set(lst.get_indices())  # FixAtoms
        __, sidx = _select_by_xyz(atom, val)

        curfix = [[[True] * 3, sidx]]
        for i in curfix:
            co = i[0]
            if co[0] and co[1] and co[2]:
                fixinds = fixinds.union(i[1])
            else:
                for at in i[1]:
                    fixscaled.update({at: co})

        natom = atom.get_global_number_of_atoms()
        ffix = [x for x in sorted(list(fixinds)) if x < natom]
        carr.append(FixAtoms(ffix))
        od = OrderedDict(sorted(fixscaled.items()))
        for k, v in od.items():
            if k in ffix or k >= natom:
                continue  # FixAtoms overrides FixScaled
            carr.append(FixScaled(atom.cell, k, v))
        atom.set_constraint(carr)
        return atom


class APlugStrain(AsecPlug):
    """Add strain."""
    
    def __init__(self, asec):
        super().__init__()
        asec.add_argument(
            self,
            "--strain",
            metavar="[FC]ratio{a-cx-z}",
            dest="strain",
            type=str,
            default=None,
            help="Strain structure. F:Fixed volume. C:add lattice constraints (For FHI-aims xyz axis only).",
        )

    def output_postfix(self, args, opt):
        return f"_S{self.safe_name(opt)}"

    def process(self, args, atom, val):
        from ase.constraints import FixCartesian

        natom = atom.get_global_number_of_atoms()
        const = []
        ore = re.compile("(?P<cmd>[A-Z]+)?(?P<ratio>[0-9.]+)(?P<axis>[a-cx-z]+)")
        orda = ord("a")
        ordx = ord("x")
        for c in val.split(","):
            mtr = np.identity(3)
            straxis = [False] * 3
            reverse = [False] * 3
            mask = [False] * 3
            x = ore.match(c)
            if x == None:
                self.piexception(f"Invalid command '{c}'...")
            cmd = x.group("cmd")
            r = float(x.group("ratio"))
            mstr = x.group("axis")
            for ms in mstr:
                om = ord(ms)
                m = -1
                if om >= orda and om <= orda + 2:
                    m = om - orda
                    reverse[m] = True
                elif om >= ordx and om <= ordx + 2:
                    m = om - ordx
                    mask[m] = True
                else:
                    self.piexception(f"Unknown axis '{ms}'...")
                    break
                mtr[m, m] = r
                straxis[m] = True

            srev = sum(reverse)
            smask = sum(mask)

            if cmd is not None and "F" in cmd:  # Fix volume
                nr = sum(straxis)
                if nr == 3:
                    self.piprint(f"Cannot fix volume due to '{mstr}'...")
                elif nr == 1 or nr == 2:
                    for im, mm in enumerate(straxis):
                        if not mm:
                            mtr[im, im] = r ** (-nr / (3 - nr))
            if srev > 0 and smask > 0:
                self.piexception(f"Mixed axis '{mstr}'...")
            cell = atom.cell
            if srev > 0:
                ncell = np.matmul(mtr, cell)
            else:
                ncell = np.matmul(cell, mtr)
            cell = ase.cell.Cell.new(ncell)

            atom.set_cell(cell, scale_atoms=True, apply_constraint=False)
            if cmd is not None and "C" in cmd:
                if smask == 0:
                    self.piprint(f"'{mstr}' do not add lattice constraints...")
                else:
                    self.piprint(f"Adding lattice constraints to '{mstr}'...")
                    for i in range(natom):
                        const.append(FixCartesian(i, mask))
                    atom.set_constraint(const)
        return atom


class APlugSurface(AsecPlug):
    """Creat a surface."""
    
    def __init__(self, asec):
        super().__init__()
        asec.add_argument(
            self,
            "--surface",
            metavar="h,k,l,Nlayers",
            dest="surface",
            type=str,
            default=None,
            help="Creat a surface of miller indices [h,k,l] with N layers.",
        )

    def output_postfix(self, args, opt):
        return f"_S{self.safe_name(opt)}"

    def process(self, args, atom, val):
        sarr = np.array(val.split(","), int)
        atom = surface(atom, sarr[:3], sarr[3], periodic=True)
        return atom


class APlugScale(AsecPlug):
    """Scale structure."""
    
    def __init__(self, asec):
        super().__init__()
        asec.add_argument(
            self,
            "--scale",
            metavar="F|Fx,Fy,Fz",
            dest="scale",
            type=str,
            default=None,
            help="Scale structure. F:Scale factor.",
        )

    def output_postfix(self, args, opt):
        return f'_SC{self.safe_name(opt,"_")}'

    def process(self, args, atom, val):
        cell = atom.get_cell()
        facts = list(map(float, val.split(",")))
        if len(facts) == 1:
            facts = facts * 3
            if args.Gslabidx > 0:
                facts[args.Gslabidx - 1] = 1
        for i, f in enumerate(facts):
            cell[i] *= float(f)
        atom.set_cell(cell, scale_atoms=True, apply_constraint=False)
        return atom


class APlugRotate(AsecPlug):
    """Rotate around an axis."""
    
    def __init__(self, asec):
        super().__init__()
        asec.add_argument(
            self,
            "--rotate",
            metavar="[C]degree[-]{a-cx-z}[,...]",
            dest="rotate",
            type=str,
            default=None,
            help="Rotate an angle `degree` around an axis `[-]{a-cx-z}`. `[C]`: With cell rotation.",
        )

    def output_postfix(self, args, opt):
        return f"_R{self.safe_name(opt)}"

    def process(self, args, atom, val):
        cell=atom.cell
        ore = re.compile("(?P<cmd>[C])?(?P<degree>[+-]?([0-9]*[.])?[0-9]+)(?P<sign>-?)(?P<axis>[x-za-c])")
        invalid=True
        for c in val.split(","):
            for x in ore.finditer(c):
                invalid=False
                cmd=x.group('cmd')
                rcell=False
                if cmd is not None and 'C' in cmd:
                    rcell=True
                ax=x.group("axis")
                sign=x.group('sign')
                rvec=sign+ax
                deg=x.group("degree")
                if (ax >='a' and ax <='c'):
                    rvec=cell[ord(ax)-ord('a')].copy()*int(sign+'1')

                atom.rotate(float(deg), rvec, rotate_cell=rcell)
        if (invalid):
            self.piexception(f"Invalid rotation expression: '{val}'...")
        return atom


class APlugAlign(AsecPlug):
    """Align a axis."""
    
    def __init__(self, asec):
        super().__init__()
        asec.add_argument(
            self,
            "--align",
            metavar="{a-c}{x-z}[,...]",
            type=str,
            default=None,
            help="Align a source axis`{a-c}` to a target axis`{x-z}`.",
        )

    def output_postfix(self, args, opt):
        return f"_A{self.safe_name(opt)}"

    def process(self, args, atom, val):
        ore = re.compile("(?P<src>[a-c])(?P<tgt>[x-z])")
        for c in val.split(","):
            for x in ore.finditer(c):
                axs = ord(x.group("src")) - ord("a")
                axt = x.group("tgt")
                atom.rotate(atom.cell[axs], axt, rotate_cell=True)
        return atom


class APlugSort(AsecPlug):
    """Sort by elements and/or positions."""
    # TODO this process is called twice
    
    def __init__(self, asec):
        super().__init__()
        asec.add_argument(
            self,
            "--elsort",
            action="store_true",
            dest="elsort",
            help="Sort by elements",
        )
        asec.add_argument(
            self, "--zsort", action="store_true", dest="zsort", help="Sort by Z,Y,X."
        )

    def output_postfix(self, args, opt):
        return f""  # _S{self.safe_name(opt)}'

    def process(self, args, atom, val):
        skey = []
        if args.zsort:
            pos = np.round(atom.positions, 2)
            skey.extend(pos.transpose())
            # skey.append(pos[:,2])
        if args.elsort:
            skey.append(atom.numbers)
        atom = atom[np.lexsort(skey)]  # last key is prime
        return atom


class APlugRepeat(AsecPlug):
    """Repeat the cell."""
    
    def __init__(self, asec):
        super().__init__()
        asec.add_argument(
            self,
            "-r",
            metavar="na,nb,nc[,opt]",
            type=str,
            default=None,
            help="Repeat the cell. n[a-c]:int/float/impr-frac, opt:-1(supercellmode)|tolerance",
        )

    def output_postfix(self, args, opt):
        return f"_{self.safe_name(opt.replace('/','%'),'x')}"

    def process(self, args, atom, val):
        rarr = val.split(",")
        zarr = [eval("/".join(map(str, map(float, self.split("/"))))) for self in rarr]

        catom = atom.get_global_number_of_atoms()
        natom = int(catom * zarr[0] * zarr[1] * zarr[2])
        tols = [0.01]
        if len(zarr) > 3:
            tols = [zarr[3]]
        tols.append(0)
        if tols[0] == -1:
            self.piprint("Supercell creating mode ...")
            atom = ase.build.make_supercell(
                atom, [[zarr[0], 0, 0], [0, zarr[1], 0], [0, 0, zarr[2]]], wrap=False
            )
        else:
            for tol in tols:
                tatom = ase.build.cut(
                    atom,
                    a=(zarr[0], 0, 0),
                    b=(0, zarr[1], 0),
                    c=(0, 0, zarr[2]),
                    tolerance=tol,
                )  # to
                nnatom = tatom.get_global_number_of_atoms()
                if natom != nnatom:
                    self.piprint("tot atom(%d) != expected(%d)'" % (nnatom, natom))
                    if tol != 0:
                        self.piprint(" Trying 0 tolereance...")
                else:
                    atom = tatom
                    break

        return atom


class APlugVacuum(AsecPlug):
    """Add vacuum to c axis."""
    
    def __init__(self, asec):
        super().__init__()
        asec.add_argument(
            self,
            "--vadd",
            metavar="vacuum",
            type=float,
            default=None,
            help="Add vacuum(A) to c axis.",
        )

    def output_postfix(self, args, opt):
        return f"_V{opt}"

    def process(self, args, atom, val):
        if atom.cell.volume == 0:
            pos = atom.get_positions()
            lo = pos.min(axis=0)
            hi = pos.max(axis=0)
            atom.set_cell(hi - lo + val)
            atom.set_positions(pos - lo, apply_constraint=False)
        else:
            ase.build.add_vacuum(atom, float(val))
            # atom.center(axis=(2))
        return atom


class APlugSetVac(AsecPlug):
    """Set actual vacuum along the z axis."""
    
    def __init__(self, asec):
        super().__init__()
        asec.add_argument(
            self,
            "--vset",
            metavar="setvac",
            type=float,
            default=None,
            help="Set actual vacuum(A) in z axis (auto aligned).",
        )

    def output_postfix(self, args, opt):
        return f"_SV{self.safe_name(opt)}"

    def process(self, args, atom, val):
        sidx = gu.identify_slabaxis(atom, True)  # Wrap separted slab
        atom = gu.align_slabaxistoz(atom, sidx)
        pos = atom.get_positions()
        lo = pos.min(axis=0)
        hi = pos.max(axis=0)
        zdiff = (hi - lo)[2]
        zlen = atom.cell.lengths()[2]
        ase.build.add_vacuum(atom, float(val) - zlen + zdiff)
        return atom


class APlugTranslate(AsecPlug):
    """Translate atom positions."""
    
    def __init__(self, asec):
        super().__init__()
        asec.add_argument(
            self,
            "-l",
            metavar="dx,dy,dz",
            type=str,
            default=None,
            help="Translate atoms of d{x-z}(A)",
        )

    def output_postfix(self, args, opt):
        return f"_T{self.safe_name(opt,'_')}"

    def process(self, args, atom, val):
        atom.translate(np.array(val.split(","), float))
        return atom


class APlugSelect(AsecPlug):
    """Select atoms."""
    
    def __init__(self, asec):
        super().__init__()
        asec.add_argument(
            self,
            "--sel",
            metavar='"expr of {x-z}"',
            type=str,
            help="Select atoms of a boolean expression. Ex) `(z>0)&(x>3)`",
        )

    def output_postfix(self, args, opt):
        return f"_S@{self.safe_name(opt)}"

    def process(self, args, atom, val):
        satom, __ = _select_by_xyz(atom, val)
        return satom


class APlugHex2Rec(AsecPlug):
    """Convert current hexagonal cell to a rectangular cell."""
    
    def __init__(self, asec):
        super().__init__()
        asec.add_argument(
            self, "--torec", action="store_true", help="Hexagonal to rectangular cell."
        )

    def output_postfix(self, args, opt):
        return f"_Rect"

    def process(self, args, atom, val):
        shift = -1
        par = atom.cell.cellpar()
        for i in range(2, -1, -1):
            ia = (i + 1) % 3
            ib = (i + 2) % 3
            if abs(par[i + 3] - 60) < 0.1 and abs(par[ia] - par[ib]) < 0.1:
                shift = (i + 1) % 3
                break

        if shift < 0:
            self.piexception(f"No hexagonal plane '{par}' ...")

        vecs = [[1, 1, 0], [-1, 1, 0], [0, 0, 1]]  # when a=b and gamma=60
        rvec = np.roll(vecs, shift, axis=1)
        atom = ase.build.cut(atom, a=rvec[0], b=rvec[1], c=rvec[2], tolerance=0.1)
        # atom.rotate(atom.cell[0],'x',rotate_cell=True)
        return atom


class APlugWrap(AsecPlug):
    """Wrap atom positions into the cell."""
    
    def __init__(self, asec):
        super().__init__()
        asec.add_argument(self, "-w", action="store_true", help="[Wrap](https://wiki.fysik.dtu.dk/ase/ase/atoms.html#ase.Atoms.wrap) atom positions into the cell.")
        asec.add_argument(self, "--wp", action="store_true", help="Wrap with pretty.")

    def output_postfix(self, args, opt):
        return f""

    def process(self, args, atom, val):
        atom.wrap(pretty_translation=args.wp)
        return atom


class APlugWSlab(AsecPlug):
    """Translate and wrap a separated slab into middle."""

    def __init__(self, asec):
        super().__init__()
        asec.add_argument(
            self,
            "--twslab",
            action="store_true",
            dest="twrapslab",
            help="Translate and wrap a separated slab into middle.",
        )

    def output_postfix(self, args, opt):
        return f""  # _S{self.safe_name(opt)}'

    def process(self, args, atom, val):
        slabidx = gu.identify_slabaxis(atom, True)
        return atom


class APlugCTrim(AsecPlug):
    """Remove dangling carbon atoms and/or hydrogenate."""
    
    def __init__(self, asec):
        super().__init__()
        asec.add_argument(
            self,
            "--ctrim",
            metavar="{''|neighcut|hy}[,...]",
            type=str,
            help="Remove dangling carbon atoms [and hydrogenate(hy)]. `neighcut` (def:3.4) ",
        )

    def output_postfix(self, args, opt):
        return f"_CT{self.safe_name(opt)}"

    @staticmethod
    def rot_vec(v, th):
        theta = np.radians(th)
        c, s = np.cos(theta), np.sin(theta)
        R = np.array(((c, -s, 0), (s, c, 0), (0, 0, 1)))
        return np.matmul(R, v.T)

    def process(self, args, atom, val):
        import ase.neighborlist as ann

        hydro = False
        cut = 3.4
        for o in val.split(","):
            if o == "":
                continue
            if o == "hy":
                hydro = True
            else:
                cut = float(o)

        natoms = atom.get_global_number_of_atoms()
        eles = np.array(atom.get_chemical_symbols())
        pos = atom.positions
        cutoffs = [cut / 4] * natoms
        nl = ann.NeighborList(cutoffs, self_interaction=False, bothways=True)
        nl.update(atom)

        hpos = []
        dpos = []
        for aid in range(natoms):
            if eles[aid] != "C":
                continue
            indices, offsets = nl.get_neighbors(aid)
            lind = len(indices)
            bpos = pos[aid]
            if lind == 2:
                vec = 2 * bpos - (pos[indices[0]] + pos[indices[1]])
                vlen = np.linalg.norm(vec)
                newpos = vec * 1.09 / vlen + bpos
                hpos.append(newpos)
            elif lind == 1:
                vec = bpos - pos[indices[0]]
                vlen = np.linalg.norm(vec)
                vec1 = self.rot_vec(vec, 30)
                vec2 = self.rot_vec(vec, -30)
                dpos.append(aid)
                vec3 = pos[indices[0]] + vec * 1.09 / vlen
                hpos.append(vec3)

        if len(dpos) > 0:
            del atom[dpos]
        if len(hpos) > 0 and hydro:
            atom.extend(ase.Atoms(f"H{len(hpos)}", hpos))

        return atom


class APlugRmLayer(AsecPlug):
    """Remove layers along the `slabaxis`."""
    
    def __init__(self, asec):
        super().__init__()
        asec.add_argument(
            self,
            "--rmlayer",
            metavar="l1[,...][,mingap(def:1.9,float)]",
            type=str,
            help="Remove layeres along the `slabaxis`. The layer index(`ln`) starts from the bottom layer, 0.",
        )

    def output_postfix(self, args, opt):
        return f"_RL{self.safe_name(opt)}"

    def process(self, args, atom, val):
        sidx = args.SlabIdx
        if sidx <= 0:
            self.piprint(f"Slab index is not posive({sidx}), please use --sidx")
            return atom
        irm = {}
        mingap = 1.9
        for o in val.split(","):
            if o.isdigit():
                irm.update({int(o): 1})
            else:
                mingap = float(o)
        _, _, ilayers = gu.identify_layers(atom, sidx, mingap)
        # print(ilayers)
        irmatom = []
        nlayers = len(ilayers)
        for i in irm.keys():
            if i >= nlayers:
                self.piprint(f"'{i}' layer doesn't exist... (l<{nlayers})")
                continue
            irmatom.extend(ilayers[i])
        if len(irmatom) > 0:
            del atom[irmatom]
        return atom


# Not added
class NoPlugTemplate(AsecPlug):
    """Plugin class template. Not loaded."""
    def __init__(self, asec):
        super().__init__()
        asec.add_argument(self,'--template', '-l', metavar='temp', type=str)

    def output_postfix(self, args, opt):
        return f"_S{self.safe_name(opt)}"

    def process(self, args, atom, val):
        return atom
