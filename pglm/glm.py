"""
GLM library interface - runs GLM in-process via libglm.so.
"""

import os
import ctypes
from pathlib import Path
from typing import Optional


def _find_libglm() -> str:
    """Find libglm.so or libglm.dylib."""
    if os.name == "nt":
        names = ["libglm.dll", "glm.dll"]
    elif os.uname().sysname == "Darwin":
        names = ["libglm.dylib", "libglm.so"]
    else:
        names = ["libglm.so"]

    search_dirs = []
    if "GLM_LIB" in os.environ:
        search_dirs.append(os.environ["GLM_LIB"])
    # glm.py lives at <glm_repo>/pglm/glm.py; the Makefile builds libglm.so
    # into the repo root (<glm_repo>), which is parent.parent from here.
    glm_root = Path(__file__).resolve().parent.parent
    search_dirs.extend([
        glm_root,
        Path.cwd(),
    ])
    search_dirs.extend(os.environ.get("LD_LIBRARY_PATH", "").split(":"))

    for d in search_dirs:
        if not d:
            continue
        d = Path(d) if isinstance(d, str) else d
        if not d.exists():
            continue
        for name in names:
            p = d / name
            if p.exists():
                return str(p)

    for name in names:
        try:
            ctypes.CDLL(name)
            return name
        except OSError:
            pass

    raise FileNotFoundError(
        f"libglm not found. Build with: make libglm.so\n"
        f"Searched: {search_dirs}"
    )


class GLM:
    """
    GLM (General Lake Model) as a Python library.

    Runs GLM in-process - no subprocess. Use like the glm executable:
    set working directory, then run(nml_path).

    Example:
        glm = GLM()
        glm.run(workspace="model_ws", nml_path="glm3.nml")
    """

    def __init__(self, lib_path: Optional[str] = None):
        path = lib_path or _find_libglm()
        self._lib = ctypes.CDLL(path)
        self._lib.glm_run_with_nml.argtypes = [ctypes.c_char_p]
        self._lib.glm_run_with_nml.restype = None

    def run(
        self,
        nml_path: str = "glm3.nml",
        workspace: Optional[str] = None,
    ) -> None:
        """Run GLM simulation."""
        orig = None
        if workspace:
            orig = os.getcwd()
            os.chdir(workspace)
        try:
            nml = os.path.abspath(nml_path) if os.path.isabs(nml_path) else nml_path
            self._lib.glm_run_with_nml(nml.encode("utf-8"))
        finally:
            if orig:
                os.chdir(orig)
