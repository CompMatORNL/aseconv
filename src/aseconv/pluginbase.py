"""
Base plugin classes for aseconv
"""

from __future__ import annotations
from typing import TYPE_CHECKING
import ase.io
import sys
from ase import Atoms
from abc import ABC, abstractmethod

if TYPE_CHECKING:  # Only imports the below statements during type checking
    from .main import AseConv

import traceback
from collections import OrderedDict


class _AsecBase(ABC):
    """Plugin base class.

    Child class should have ``_plugins`` instance. Currently two differnt plugin class cannot have the same name.

    """

    # _plugins: OrderedDict = OrderedDict()
    # _plugins = []

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)

        cname = cls.__name__
        if (
            not cname.startswith("NoPlug")
            and not cname.startswith("Asec")
            and not cname.startswith("_Asec")
        ):
            if cname not in cls._plugins:
                cls._plugins.update({cname: cls})
            else:
                print(f" - [_AsecBase]: Warn '{cname}' is duplicated...")


class AsecPlug(_AsecBase):
    """Plugin base class for ``atom`` modification."""

    _plugins: OrderedDict = OrderedDict({})
    _instances: list = []
    _piinit = False

    def __init__(self):
        """Initializing the class.

        Attributes:
            clspre: Class name prefix for output.

        """
        self.clspre: str = self.__class__.__name__.replace("APlug", "").lower()

    def init_plugins(self, asp: aseconv.main.AseConv):
        """Loading all registered _plugins."""
        if self._piinit:
            return

        for key, p in AsecPlug._plugins.items():
            self._instances.append(p(asp))

        self._piinit = True

    def piprint(self, *args, **kwargs):
        """Plugin print with ``clspre`` tag."""
        print(f" - [{self.clspre}] " + args[0], *args[1:], **kwargs)

    def piexception(self, *args, **kwargs):
        """Plugin exeception with ``clspre`` tag."""
        raise Exception(f" - [{self.clspre}] " + args[0], *args[1:], **kwargs)

    def safe_name(self, src: str, rstr: str = "") -> str:
        """Returns safe file name.

        Args:
            src: A source string.
            rstr: A replacement string for invalid characters.

        Returns:
            Corrected file name.
        """

        chars = "/,:<> ()+*{}[]!$&=|"
        for c in chars:
            if c in src:
                src = src.replace(c, rstr)
        return src

    @abstractmethod
    def process(self, atom: ase.Atoms) -> ase.Atoms:
        """Actual process method to modify ``atom``.

        Args:
            atom: An ``atom`` image.

        Returns:
            The modified ``atom`` image.

        """
        pass

    def output_postfix(self, args: argparse.Namespace, opt: str) -> str:
        """Postfix to be added to default output file name.

        Args:
            args: Argument list.
            opt: The option string from an argument.

        Return:
            Postfix to be appended.

        """
        return ""


class AsecIO(_AsecBase):
    """Base class for type IO plugins for specific extensions.

    If duplicated, this plugin will be used instead of ``ASE`` library's.

    """

    _plugins: OrderedDict = OrderedDict({})
    _instances = []
    _piinit = False

    def init_plugins(self):
        """Initialize registered _plugins."""
        if self._piinit:
            return

        for key, p in self._plugins.items():
            self._instances.append(p())

        self._piinit = True

    @abstractmethod
    def infos(self) -> dict:
        """Get IO plugin info. Must be implemented in the subclasses.

        The ``typeexts`` can include multiple dicts.
        The ``type`` is a name for the geometry format such as 'vasp', 'aims', or etc.
        The ``extwithdot`` is the corresponding extenstion including '.'.

        Returns:
            ``dict`` containing

            - 'help': 'Help message'
            - 'typeexts': { 'type1': 'ext1withdot' , ... }
        """
        return {}

    def read(self, file: str, type: str, **kwargs) -> ase.Atoms:
        """Read function for specfic extensions and returns an ``ase.Atoms`` object.

           If not implemented in a plugin class, ``ase.io.read`` is used.

        Args:
            file: Input file name.
            type: Input file type.
            kwargs: Arbitrary keyword arguments.

        Returns:
            The opened atom image.

        """
        return ase.io.read(file, format=type, **kwargs)

    def write(self, args: argparse.Namespace, atom: ase.Atoms, type: str, file: str):
        """Write function of an ``atom`` image for the plugin.

        Args:
            args: Processed arguments from ``parse_args``.
            atom: An atom image.
            type: Output file type.
            file: Output file name.

        """
        raise Exception(" - [AsecIO] write is not implemented...")

    @staticmethod
    def vec2str(vec: list, fmt: str = "{:>20.16f}") -> str:
        """Converts a vector of arbitary size into a formmated string ``str``.

        Args:
            vec: A vector of arbitary size.
            fmt: A format string for a single ``float``.

        Returns:
            The formatted string joined by a single space ' '.

        """

        return " ".join([str.format(fmt, x) for x in vec])
