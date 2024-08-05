from __future__ import annotations
import sys
from aseconv.pluginbase import AsecPlug, AsecIO

class APlugKpath(AsecPlug):
    """Kpath generator plugin class.
    
    Supported formats: vasp, aims, rmg only.
    
    Required: spg, ase, seekpath, numpy, scipy
       
    Attributes:
    
    """
    
    def __init__(self, asec):
        #parser=asec.add_subparsers('kpath', help='Kpath generator', description="KPath Generatore")
        super().__init__()
        
        parser = None
        asec.add_argument(
            self,
            "--kp",
            parser=parser,
            metavar="Opt",
            dest="kopt",
            type=str,
            default=None,
            help='Generate kpath. Options: ""(auto)|kfile(.wkp)|kpstring|ss(show)',
        )
        asec.add_argument(
            self,
            "--ksidx",
            parser=parser,
            process=False,
            metavar="KSlabIdx",
            type=int,
            default=-1,
            help="Kpath 2D slab index. -1(auto)|0(bulk)|{1-3}(a-c)",
        )

    def output_postfix(self, args, opt):
        return f"_KP{self.safe_name(opt)}"

    def process(self, args, atom, val):
        ofile = args.pOutFile
        pfile = args.pInFile
        dp = args.pParent
        kfile = dp.joinpath(ofile.stem + ".akp")
        bzfile = dp.joinpath(
            "bz_"
            + ofile.stem
            + "_%s_%s" % (pfile.parent.name, pfile.stem.replace("geometry_", "g"))
        )
        sratom = self._kp_main(args, atom, val, bzfile, kfile)
        return sratom

    def _kp_3d(self, args, atom, nums, syms):
        import seekpath
        from ase import Atoms

        self.piprint(
            """\n\
    ##############################
    ####        WARNING       ####
    ###  Cell is primitivized  ###
    ###  by k-path generator   ###
    ##############################"""
        )

        cell = (atom.cell, atom.get_scaled_positions(), nums)
        kp = seekpath.get_path(
            cell, symprec=0.5, threshold=0.5
        )  # ,with_time_reversal=False)

        pre = "primitive_"
        sr_syms = []
        for i in kp[pre + "types"]:
            sr_syms.append(syms[i - 1])
        sr_atom = Atoms(
            symbols=sr_syms,
            scaled_positions=kp[pre + "positions"],
            cell=kp[pre + "lattice"],
            pbc=True,
        )
        # sr_atom.center()
        return sr_atom, kp

    def _kp_parse(self, lines: list) -> dict:
        """
        Args:
            lines:
                [ "FromK1 K1x K1y K1z ToK2 K2x K2y K2z"
                   , ...
                ]
        
        Returns:
            kp: kpath dictionary
            
        """
        path = []
        points = {}
        for l in lines:
            ar = l.split()
            if len(ar) != 8:
                print("[WARN] '{}' is ignored".format(l))
                continue
            ii = 0
            points.update({ar[ii]: list(map(float, ar[ii + 1 : ii + 4]))})
            ii = 4
            points.update({ar[ii]: list(map(float, ar[ii + 1 : ii + 4]))})
            path.append((ar[0], ar[4]))
        kp = {"path": path, "point_coords": points}
        return kp


    def _kp_auto(self, args, atom):
        syms = []
        nums = []
        for i in atom.get_chemical_symbols():
            if i not in syms:
                syms.append(i)
            nums.append(syms.index(i) + 1)

        rotate = False
        ksidx = args.ksidx
        if ksidx == -1:
            ksidx = args.SlabIdx

        if ksidx > 0:
            import aseconv.plugins._kpath2d as kp2

            sr_atom, kp = kp2.kpath_2d(args, atom, nums, ksidx)
            # Update slab index
            ksidx = 3
            args.SlabIdx = ksidx
        else:
            sr_atom, kp = self._kp_3d(args, atom, nums, syms)

        return sr_atom, kp, ksidx


    def _kp_output(self, args, kp, inst, ksidx):
        dkp = kp["path"]
        cr = kp["point_coords"]
        # cr2d=kp.get('point_2dcoords',{})

        kf = "{:6.3f}"
        tickstr = "#TICK:"
        pes = ""
        nkp = []
        npcr = {}
        for i in dkp:
            b = i[0]
            e = i[1]

            bs = b.replace("GAMMA", "G")
            es = e.replace("GAMMA", "G")
            k1 = AsecIO.vec2str(cr[b], fmt=kf)
            k2 = AsecIO.vec2str(cr[e], fmt=kf)

            npcr.update({b: cr[b], e: cr[e]})
            nkp.append((b, e))
            if pes != "" and bs != pes:
                tickstr += pes + "|"
            inst.addpath(bs, k1, es, k2, pes)
            tickstr += bs + ","
            pes = es
            pk2 = k2

        tickstr += pes

        kp["path"] = nkp
        kp["point_coords"] = npcr

        inst.endpath(args, ksidx, tickstr)


    def _kp_main(self, args, atom, kopt, bzfile, kfile):
        import aseconv.plugins._kpathBZ as bz

        sr_atom = None
        kp = None

        kslabidx = args.ksidx
        if kopt == "" or kopt == "ss":
            sr_atom, kp, kslabidx = self._kp_auto(args, atom)
        else:
            if Path(kopt).is_file():
                lines = distutils.text_file.TextFile(filename=kopt).readlines()
            else:
                lines = kopt.split("@")
            kp = self._kp_parse(lines)
            sr_atom = atom

        dkp = kp["path"]
        if args.t == "vasp":
            inst = KPOvasp(dkp)
        elif args.t == "aims":
            inst = KPOaims(dkp)
        elif args.t == "rmg":
            inst = KPOrmg(dkp)
        else:
            self.piexception(f"Unsupported type '{args.t}' for kpath...")

        self._kp_output(args, kp, inst, kslabidx)
        inst.writepath(kfile)

        try:
            bz.draw_brillouinzone(kopt == "ss", bzfile, sr_atom, kp)
        except:
            # print("X11 server error, please launch X11:", sys.exc_info())
            print(sys.exc_info())
            self.piexception("Draw_brillouinzone error")

        return sr_atom


class KPathOut:
    """Kpath output base class
    
    Attributes:
        lpathstr: List of the path strings.
        nkp: Number of kpoints per segment.
        
    """
    def __init__(self, dkp: dict):
        self.lpathstr: list = []
        self.nkp: int = int(100 / (len(dkp) + 1)) * 4 + 1

    def addpath(self, bs: str, k1: str, es: str , k2: str , pes: str):
        """Add path.

        Args:
                bs: begin k point name.
                k1: begin k point coordinates.
                es: end k point name.
                k2: end k point coordinates.
                pes: previous es.
        """
        pass

    def endpath(self, args: argparse.Namespace, ksidx: int, tickstr: str):
        
        """Endpath

        Args:
                args: Result of ArgumentParser.parse_args().
                ksidx: k path slab index.
                tickstr: tick string.
        """
        pass

    def writepath(self, kfile: str):
        """Write kpath to ``kfile``.

        Args:
                kfile: Output file name
        """
        print(" - Writing '{}'...".format(str(kfile)))
        with open(kfile, "wt") as f:
            f.write("\n".join(self.lpathstr))
            f.write("\n")


class KPOvasp(KPathOut):
    """Kpath output class for vasp/wannier90."""
    
    def __init__(self, dkp):
        super().__init__(dkp)
        self.fmwanr = " {0:>3}  {1}  {2:>3}  {3}"
        self.fmvasp = "{:<20}  ! {}"
        self.lstrwanr = []
        self.lstrvasp = []

    def addpath(self, bs, k1, es, k2, pes):
        self.lstrwanr.append(self.fmwanr.format(bs, k1, es, k2))
        self.lstrvasp.extend(
            [self.fmvasp.format(k1, bs), self.fmvasp.format(k2, es), ""]
        )

    def endpath(self, args, ksidx, tickstr):
        lret = self.lpathstr
        lkmask = [1] * 3
        lamask = [1] * 3

        sidx = args.SlabIdx
        if ksidx > 0:
            lkmask[ksidx - 1] = 0
        if sidx > 0:
            lamask[sidx - 1] = 0
        lret.append(
            "#SLABMASK: " + " ".join([str(x) for x in lamask])
        )  # This is for KPOINTS

        lmall = []
        if sidx == 0:  # Original Bulk, but 2-D Kpath
            if ksidx > 0:
                lm = ["0"] * 3
                lm[ksidx - 1] = "1"
                lmall.append(lm)
            else:
                # TODO: get user input
                lmall.extend([list("001"), list("100"), list("010")])
        else:
            for k in range(3):
                v = lkmask[k]
                if v == 0:
                    lm = ["1"] * 3
                    lm[k] = "0"
                else:
                    lm = ["0"] * 3
                    lm[k] = "1"
                lmall.append(lm)

        for lm in lmall:
            lret.append("#MILLER_{}: ".format("_".join(lm)) + " ".join(lm))

        WAN90_KEY0 = "#<WANNIER90>"
        WAN90_KEY1 = "#<@WANNIER90>"
        aswanr = self.lstrwanr

        lret.extend(
            [
                WAN90_KEY0,
                "begin kpoint_path",
                *aswanr,
                "end kpoint_path",
                WAN90_KEY1,
                "",
            ]
        )
        lret.extend(
            [
                "#<VASP>",
                "k-points along high symmetry lines  " + tickstr,
                str(self.nkp),
                "Line-mode",
                "rec",
            ]
        )
        lret.extend(self.lstrvasp)
        lret.extend(["#<@VASP>"])
        lret.extend(["", tickstr, ""])
        lret.extend(
            [
                "!<WANTOOL>",
                "KPATH_BULK",
                str(len(aswanr)),
                *aswanr,
                "",
                "!<@WANTOOL>",
                "",
            ]
        )


class KPOaims(KPathOut):
    """KPath output class for FHI-aims."""
    
    def __init__(self, dkp):
        super().__init__(dkp)
        self.fmt = "output band  {1} {3} {4:>4} {0:>3}  {2:>3}"  # 21: max for aims?

    def addpath(self, bs, k1, es, k2, pes):
        self.lpathstr.append(self.fmt.format(bs, k1, es, k2, self.nkp))

    def endpath(self, args, ksidx, tickstr):
        self.lpathstr.append(tickstr)


class KPOrmg(KPathOut):
    """KPath output class for RMG."""
    
    def __init__(self, dkp):
        super().__init__(dkp)
        self.lpathstr.append('kpoints_bandstructure = " ')
        self.fmt = "  {1} {2:>4} {0:>3}"
        self.nkp = int(self.nkp * 2 / 3)

    def addpath(self, bs, k1, es, k2, pes):
        nrs = self.nkp
        if pes != "" and bs != pes:
            self.lpathstr.append(self.fmt.format(pes, k2, nrs))
            nrs = 0
        self.lpathstr.append(self.fmt.format(bs, k1, nrs))
        self.es = es
        self.k2 = k2

    def endpath(self, args, ksidx, tickstr):
        self.lpathstr.append(self.fmt.format(self.es, self.k2, self.nkp))
        self.lpathstr.append(' "')

def ___kpath_make_pos(args, kp):
    pre = "primitive_"
    pos = kp[pre + "positions"]
    if sum(sum(pos)) > 0:
        return kp
    cell = kp[pre + "lattice"]
    nkp = kp.copy()
    nkp[pre + "positions"] = -pos
    # nkp[pre+'lattice']=-cell
    print("   [WARN] Kpath negative pos(bug?) were inverted")
    return nkp


def ___kpath_add(args, atom, orgkp):
    # print(orgkp)
    pts = args.kadd.split(",")
    parr = [x.split(":") for x in pts]
    plist = []
    for k in parr:
        plist.append([k[0], np.array([float(x) for x in k[1:]])])
    ptc = orgkp["point_coords"]
    ptp = orgkp["path"]
    cellT = np.array(atom.cell.reciprocal()).T
    validp = []
    for i, bpl in enumerate(plist):
        bpx = bpl[1]
        # remove symmetry points
        over = False
        for k, bpc in ptc.items():
            sqsum = ((np.array(bpc) - bpx) ** 2).sum()
            # print(k,bpc,sqsum)
            if sqsum < 1e-5:
                over = True
                break
        if over:
            continue
        cpx = np.dot(cellT, bpx[1])
        intheline = False
        # remove point in path
        for pi, pv in enumerate(ptp):
            bpa = np.array(ptc[pv[0]])
            bpb = np.array(ptc[pv[1]])
            cvab = cellT.dot(bpb - bpa)
            cvax = cellT.dot(bpx - bpa)
            prox = np.linalg.norm(np.cross(cvab, cvax))
            if prox < 1e-5:  # if colinear
                dotax = cvab.dot(cvax)
                dotab = cvab.dot(cvab)
                if dotax >= 0 and dotax <= dotab:
                    intheline = True
                    break
        if not intheline:
            validp.append(bpl)

    lvalid = len(validp)
    if lvalid == 1:
        # find a closest symmetry point
        valp = validp[0]
        mdist = 9999999
        mpoint = []
        for k, v in ptc.items():
            nv = np.array(v)
            dist = np.linalg.norm(np.dot(cellT, (valp[1] - nv)))
            if dist < mdist:
                mdist = dist
                mpoint = [k, nv]
        # insert path
        pidx = len(ptp) - 1
        for i, p in enumerate(ptp):
            if p[1] == mpoint[0]:
                pidx = i

        if pidx < len(ptp) - 1:  # insert middle
            ptp.insert(pidx + 1, (mpoint[0], valp[0]))
            ptp.insert(pidx + 2, (valp[0], mpoint[0]))
        else:  # add to the end
            ptp.insert(pidx + 1, (mpoint[0], valp[0]))
    elif lvalid == 2:
        ptp.append((validp[0][0], validp[1][0]))

    for p in validp:
        ptc.update({p[0]: tuple(p[1])})

    return orgkp
