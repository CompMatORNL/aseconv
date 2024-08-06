#!/usr/bin/env python3
"""
Aseconv program
"""

# export PYTHONDONTWRITEBYTECODE=1

from aseconv import geoutil as gu
from aseconv.pluginbase import AsecPlug, AsecIO
import aseconv.plugins

import argparse
import sys, os, glob, re, io
from pathlib import Path
import subprocess

import ase.io, ase.build
import ase.io.formats as afmt

import textwrap
import traceback

class AseConv:
    """Main class of the utility."""

    class _AsecHelpFormatter(
        argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter
    ):
        pass

    def __init__(self):
        """Class init
            
            Attributes:
        """
        # argparse.RawTextHelpFormatter
        self.parser = argparse.ArgumentParser(
            description=__doc__,
            formatter_class=self._AsecHelpFormatter,
            allow_abbrev=False,
            exit_on_error=False,
        )
        self._plugins = {}
        self._iopitypeinst = {}
        self._iopitypeext = {}
        self._iopiextinst = {}
        self._instios = None
        self._subparsers = self.parser.add_subparsers(title="commands")
        server_desc = """\
 Simple aseconv server for faster multi processing.
 
 When processing files multiple times, python loading time could be a bottle neck.
 This server runs in background with an idle shutdown timeout 30s (default),
 and receives commands through localhost:17399 (default).
 Commands can be sent using any TCP clinets, but the following example could be used 
 without installing any other software. 
 
 - An example client command in linux.

 `MyCOMMANDS` needs to be replaced properly.
 `--pwd` sets the path of input files, and usually `${PWD}`.
```
    ASEC_SERVER=/dev/tcp/localhost/17399;AFD=793 
    # AFD could be any number (file descriptor).
    exec {AFD}<>${ASEC_SERVER} && echo -n "geo --pwd '${PWD}' MyCOMMANDS [...]" >&${AFD}
    # Open the RW(<>) TCP port and send commands using `echo`.
    ret=1; read -r -u ${AFD} ret; rr=$?; (( rr > 128 )) && ret=1
    # Read a result from server. 0 for success. 
    exec {AFD}<&-   # Port close. Important. 
```
    
"""
        parser_a = self._subparsers.add_parser(
            "service",
            help="asec service",
            description=server_desc,
            formatter_class=self._AsecHelpFormatter,
        )
        parser_a.add_argument("--port", type=int, default=17399, help="Server TCP port.")
        parser_a.add_argument(
            "--timeout", type=int, default=30, help="Idle shutdown timeout(s)."
        )
        parser_a.add_argument(
            "--stop", action="store_true", help="Send stop signal to the server."
        )
        parser_a.set_defaults(func=self._server)

        self.pparser = argparse.ArgumentParser(add_help=False)
        self.pparser.add_argument(
            "-t",
            metavar="OutputFormat",
            type=str,
            required=True,
            help="Output Format.",
        )
        self.pparser.add_argument(
            "ins", metavar="File/Dir", type=str, nargs="+", help="File/Dir of sources."
        )
        self.pparser.add_argument(
            "-i", metavar="InputFormat", type=str, default=None, help="Input Format."
        )
        self.pparser.add_argument(
            "-C", action="store_true", help="Save in Cartesian coordinates (def:fractional)."
        )
        self.pparser.add_argument(
            "-f", "--force", dest="f", action="store_true", help="Force write."
        )
        self.pparser.add_argument(
            "--frame",
            metavar="start:[stop[:step]]",
            type=str,
            default=":",
            help="Specify frame(s) using python array slicing for multiframe format file such as `.xyz`. ",
        )
        self.pparser.add_argument(
            "-o", metavar="Outputfile", type=str, default=None, help="Output file name."
        )
        self.pparser.add_argument(
            "--pwd",
            metavar="PWD",
            dest="pwd",
            type=str,
            default="",
            help="Working path when using `service`.",
        )
        self.pparser.add_argument(
            "--sidx",
            metavar="SlabIdx",
            type=int,
            default=-1,
            help="-1(auto), 0(non-slab), {1-3}(a-c axis).",
        )
#        self.pparser.add_argument("--aslab", action="store_true", help="Temporary")
   
        self.parser_geo = self.add_subparsers(
            "geo",
            help="geometry modificatin",
            description="Geometry modifier",
        )
    def add_subparsers(self, cmd: str, **kwargs) -> argparse.ArgumentParser:
        """Add subparsers for plugins. 
        
        Args:
            cmd: Command for subparser. 
            kwargs: Keywoard args for ``argparse.add_parser``
            
        Returns:
             Added subparser.
        """
        parser = self._subparsers.add_parser(
            cmd,
            parents=[self.pparser],
            formatter_class=self._AsecHelpFormatter,
            allow_abbrev=False,
            **kwargs,
        )
        parser.set_defaults(func=self._main_handler, argparser=parser)
        return parser

    def add_argument(self, clplug: AsecPlug, *args, parser: argparse.ArgumentParser=None, process: bool=True, **kwargs):
        """Add argument when registering a plugin class.
        
        Args:
            clplug: Plugin class. 
            parser: A parser for the arguments.
            process: Whther to be used for the atom image processing. True for yes and false for no.
            kwargs: Keyword args for ``argparse.add_argument``
            
        Returns:
             Added subparser.
        """
       
        if parser is None:
            parser = self.parser_geo
        if process:
            prog = parser.prog
            if prog not in self._plugins:
                self._plugins.update({prog: {}})
            dplugs = self._plugins[prog]
            for v in args:
                dplugs.update({v.lstrip("-"): clplug})
        parser.add_argument(*args, **kwargs)

    def _ordered_loop(self, sargs, args, mode, ret=""):
        dplugs = self._plugins[args.argparser.prog]
        pitems = dplugs.keys()
        nargs = len(sargs)
        for i, v in enumerate(sargs):
            if not v.startswith("-"):
                continue
            vv = v.lstrip("-")
            if not vv in pitems:
                continue
            cls = dplugs[vv]
            if i < nargs - 1:
                vopt = sargs[i + 1]
            else:
                vopt = ""
            if mode == "pfix":
                ret += cls.output_postfix(args, vopt)
            elif mode == "process":
                ret = cls.process(args, ret, vopt)

        return ret

    def _onefile(self, sargs, args, inp):
        din = Path(inp)
        outset = None

        # searching for plugin type
        oext = ""  # without dot
        for type, ext in self._iopitypeext.items():
            if args.t == type:
                oext = ext
                if oext.startswith('.'):
                    oext=oext[1:]
                break
        if oext == "":
            fmt = afmt.get_ioformat(args.t)
            if len(fmt.extensions) > 0:
                oext = fmt.extensions[0] # No dot
            else:
                oext = fmt.name

        pfix = self._ordered_loop(sargs, args, "pfix")
        if din.is_dir():
            if args.o != None:
                dp = Path(args.o)
            else:
                dp = din.parent.joinpath("0conv_" + args.t + "_" + din.name)
                
            files = [
                x
                for x in glob.glob(str(din.joinpath("*")))
                if Path(x).name != "desktop.ini"
            ]
        else:
            if args.o != None:
                dp = Path(args.o).parent
                outset = Path(args.o)
            else:
                dp = din.parent
            files = [str(inp)]

        dp.mkdir(parents=True, exist_ok=True)
        args.pParent = dp

        for f in files:
            pfile = Path(f)
            if outset != None:
                ofile = outset
            else:
                ofile = dp.joinpath(pfile.stem + pfix + '.' + oext)

            isdev = str(ofile).startswith("/dev")
            # TODO folder type plugin
            if (
                not args.f and ofile.exists() and not isdev
            ):  
                print("[INFO] '{}' exists...".format(str(ofile)))
                continue
            if not pfile.exists():
                print("[ERR] No input file({}) exists...".format(str(pfile)))
                continue

            args.pOutFile = ofile
            args.pInFile = pfile

            allatom = self._read(args, pfile)
            if ofile != pfile and not isdev:
                ofile.unlink(missing_ok=True)
            # global args.Gslabidx
            outatoms = []
            for atom in allatom:
                slabidx = args.sidx
                if slabidx < 0:
                    slabidx = gu.identify_slabaxis(atom)
                args.SlabIdx = slabidx
                atom = self._ordered_loop(sargs, args, "process", atom)
                self._write(args, atom, ofile)

    def _read(self, args, pfile):
        ext = pfile.suffix
        type = args.i
        defkwargs = {"index": args.frame, "do_not_split_by_at_sign": True}
        if type is not None:
            if type in self._iopitypeinst:
                cls=self._iopitypeinst[type]
                return cls.read(pfile,type,**defkwargs)
        else:
            if ext in self._iopiextinst:
                cls = self._iopiextinst[ext]
                return cls.read(pfile,type,**defkwargs)

        return ase.io.read(pfile, format=type, **defkwargs)

    def _write(self, args, atom, ofile, log=True):
        type = args.t
        if log:
            print(
                " - Writing [{}] '{}'...".format(
                    atom.get_chemical_formula(), str(ofile)
                )
            )
        # if (args.noconst):
        # 	atom.set_constraint(None)
        if type in self._iopitypeinst:
            cls = self._iopitypeinst[type]
            return cls.write(args, atom, type, ofile)

        parm = {}
        if type == "vasp":
            parm = {"direct": not args.C, "wrap": False}
        ase.io.write(ofile, atom, format=type, **parm)

    def _main_handler(self, sargs, args):
        # lmpfix=self.read_fix(args)
        for i in args.ins:
            print(f">>> Globbing '{i}' ...")
            f = glob.glob(i)
            if len(f) < 1:
                print(f" - No file in '{i}/' ...")
                return 1

            for j in f:
                if len(f) > 1:
                    print(" ====")
                print(" > Processing '{}'...".format(j))
                self._onefile(sargs, args, j)  # ,lmpfix)
        return 0

    def _plug_update_warn(self, type, ext, inst):
        if type in self._iopitypeinst:
            print(f"[WARN] '{type}' is already registered to {self._iopitypeinst[type]}...")
        self._iopitypeinst.update({type: inst})
        self._iopiextinst.update({ext: inst})
        self._iopitypeext.update({type: ext})
        

    def _server(self, sargv, args):
        import socketserver, socket, queue, threading, shlex, traceback

        class _TCPHandler(socketserver.BaseRequestHandler):
            def __init__(self, mainself):
                self.mainself = mainself

            def __call__(self, *args, **kwargs):
                super().__init__(*args, **kwargs)

            def handle(self):
                req = self.request
                msg = str(req.recv(4096), "ascii")
                req.settimeout(0)
                if msg != "":
                    self.mainself.qmsg.put(msg)
                    ret = self.mainself.qret.get(msg)
                    response = bytes(ret + "\n", "ascii")
                    req.sendall(response)
                # req.close()

        port = args.port
        deftimeout = args.timeout

        if args.stop:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as soc:
                try:
                    print(f" - ## Sending stop to server ::{port}")
                    soc.connect(("localhost", port))
                    soc.sendall(b"stop")
                    response = soc.recv(4096)
                except:
                    # print(traceback.format_exc())
                    print(" - ## No server is running...")
            sys.exit(1)

        print(f" - ## Starting aseconv server ::{port}, timeout ({deftimeout}s)")
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as soc:
            try:
                soc.bind(("localhost", port))
            except:
                print(" - ## binding error...")
                print(traceback.format_exc())
                sys.exit(1)

        self.qmsg = queue.Queue()
        self.qret = queue.Queue()

        # server.daemon_threads=True
        socketserver.TCPServer.allow_reuse_address = True
        server = socketserver.ThreadingTCPServer(("localhost", port), _TCPHandler(self))
        # server = socketserver.ThreadingTCPServer(('localhost',port), TCPHandler(self),bind_and_activate=False)
        with server:
            # server.allow_reuse_address = True
            # server.daemon_threads = True
            # server.server_bind()
            # server.server_activate()

            server_thread = threading.Thread(target=server.serve_forever)
            server_thread.daemon = True
            server_thread.start()
            timeout = deftimeout
            while True:
                try:
                    msg = self.qmsg.get(timeout=1)
                    timeout = deftimeout
                    spm = shlex.split(msg)
                    print(f" - ## Requested args: {msg}")
                    if spm[0] == "stop":
                        break
                    hargs = self.parser.parse_args(spm)
                    if hargs.pwd == "":
                        os.chdir(hargs.pwd)
                    ret = hargs.func(spm, hargs)
                    self.qret.put(str(ret))
                except queue.Empty:
                    timeout -= 1
                    if timeout <= 0:
                        break
                    continue
                except KeyboardInterrupt:
                    break
                except argparse.ArgumentError:
                    self.qret.put(str(1))
                    self.parser.print_help()
                except:
                    # exit_on_error not working (SysExit from argparse)
                    self.qret.put(str(1))
                    print("Unexpected error:", sys.exc_info()[0])
                    # traceback.print_exc()

            # print(dir(server))
            # server.shutdown()
            # print(" - ## aseconv server was closed...")
            self.qret.put(str(0))
        print(" - ## aseconv server was closed...")
        sys.exit(0)

    def init_parser(self, iaio: AsecIO):
        """Initializer parser.
        
            Args:
                iaio: Instance of ioplugins. 
        """
        self._instios = iaio._instances
        for i in self._instios:
            for type, ext in i.infos().get("typeexts", {}).items():
                self._plug_update_warn(type, ext, i)

        lstio=[]
        aiodict=ase.io.formats.ioformats
        for x in aiodict:
            if (x not in self._iopitypeext):
                lstio.append(x)
        
        #print(os.get_terminal_size())
        #OSError: [Errno 25] Inappropriate ioctl for device
        #w=os.get_terminal_size().columns
        # Description for supported formats.
        defw=80
        atype=textwrap.fill(f"**ASE** supported formats: {', '.join(lstio)}", width=defw, subsequent_indent="        ")
        #atype=', '.join(lstio)
        ptype=textwrap.fill(f"**aseconv** supported formats: {', '.join(self._iopitypeext)}", width=defw, subsequent_indent="        ")
        geo_desc="""Geometry modification tool using [ASE](https://wiki.fysik.dtu.dk/ase/index.html) library.
    
{}
    
{}

    """.format(atype, ptype)
        self.parser_geo.description=geo_desc
        
        return self.parser

    def start(self):
        """Main start function."""
        
        try:
            args = self.parser.parse_args()
            if len(sys.argv) == 1:
                self.parser.print_help()
                sys.exit(1)
                
            args.func(sys.argv, args)
        except argparse.ArgumentError:
            self.parser.print_help()
        #except:
        #   traceback.print_exc()
        #  sys.exit(1)

ascparser=None
def acmparser():
    """Init function for sphinx-argparse"""
    
    global asp, ascparser

    if (ascparser is not None):
        return ascparser
    
    asp = AseConv()
    class NoPlugMain(AsecPlug):
        def process(self, atom):
            pass
    
    class NoPlugIO(AsecIO):
        def infos(self):
            return {}
            
    inst = NoPlugMain()
    inst.init_plugins(asp)
    iaio = NoPlugIO()
    iaio.init_plugins()
    ascparser=asp.init_parser(iaio)
    return ascparser
    
def main():
    """Main function."""
    
    global asp
    acmparser()
    asp.start()


if __name__ == "__main__":
    main()
