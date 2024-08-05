from aseconv.pluginbase import AsecIO

class AimsIO(AsecIO):
    """FHI-aims IO plugin to enhance ``ase``'s write function.
    
    This plugin can convert constraints read from VASP file to FHI-aims.
    
    """
    def __init__(s):
        pass

    def infos(s):
        return {"help": "FHI-aims plugin type", "typeexts": {"aims": ".in"}}

    def _const_FixScaled(arr, kw):
        idx = kw["a"]
        mask = kw["mask"]
        arr[idx] = mask

    def _const_FixAtoms(arr, kw):
        ind = kw["indices"]
        for i in ind:
            if arr[i]:
                print("[WARN] [{}] is fixed by other...".format(i))
            arr[i] = [True] * 3

    def write(s, args, atom, type, ofile):
        iscart = args.C
        cell = atom.get_cell()
        outstr = []
        sort_atom = atom[atom.numbers.argsort()]
        outstr.append("# " + str(sort_atom.symbols))
        conrel = "  constrain_relaxation  "
        if atom.cell.volume != 0:
            # TODO: try not to use args.strain
            # args.strain is not None and
            if len(atom.constraints) > 0:
                mask = atom.constraints[1].mask
                for i in cell:
                    outstr.append("lattice_vector " + s.vec2str(i))
                    for j, m in enumerate(mask):
                        if not m:
                            outstr.append(conrel + chr(ord("x") + j))
                atom.set_constraint(None)
            else:
                for i in cell:
                    outstr.append("lattice_vector " + s.vec2str(i))
        else:
            iscart = True
        # outstr.append('')
        if iscart:
            pos = atom.get_positions()
            atom_str = "atom"
        else:
            pos = atom.get_scaled_positions(wrap=False)
            atom_str = "atom_frac"
        const_arr = [None] * atom.get_global_number_of_atoms()
        err = 0
        for i in atom.constraints:
            d = i.todict()
            name = d["name"]
            func = globals().get("aims_const_" + name, None)
            if not func:
                print(
                    " [ERR] Unhandled constraint exists({}), please report.".format(name)
                )
                err += 1
                continue
            func(const_arr, d["kwargs"])

        if err > 0:
            try:
                ofile.unlink()
            except OSError as e:
                pass
            return

        sym = atom.get_chemical_symbols()

        for i, xyz in enumerate(pos):
            lstr = atom_str + "   " + s.vec2str(xyz) + " " + sym[i]
            outstr.append(lstr)
            con = const_arr[i]
            if not con:
                continue
            if con[0] and con[1] and con[2]:
                lstr = conrel + ".true."
                outstr.append(lstr)
                continue

            for j, xx in enumerate(con):
                if not xx:
                    continue
                lstr = conrel + chr(ord("x") + j)
                outstr.append(lstr)

        # print(" - Writing '{}'...{}".format(str(ofile)," "*20),end=end)
        with open(ofile, "wt") as f:
            f.write("\n".join(outstr))
            f.write("\n")
