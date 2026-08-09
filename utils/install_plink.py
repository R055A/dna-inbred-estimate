from urllib.request import Request, urlopen
from stat import S_IXUSR, S_IXGRP, S_IXOTH
from tempfile import TemporaryDirectory
from platform import system, machine
from zipfile import ZipFile, ZipInfo
from urllib.error import URLError
from shutil import copyfileobj
from pathlib import Path
from enum import Enum
from os import name


class ValidOpSys(Enum):
    LINUX = "linux"
    DARWIN = "darwin"
    WINDOWS = "windows"


class ValidArch(Enum):
    BIT_32 = "32"
    BIT_64 = "64"


PLINK_STABLE_BUILD: str = "20250819"
PLINK_DOWNLOAD_URLS = {
    (ValidOpSys.LINUX, ValidArch.BIT_64): (
        f"https://s3.amazonaws.com/plink1-assets/"
        f"plink_linux_x86_64_{PLINK_STABLE_BUILD}.zip"
    ),
    (ValidOpSys.LINUX, ValidArch.BIT_32): (
        f"https://s3.amazonaws.com/plink1-assets/"
        f"plink_linux_i686_{PLINK_STABLE_BUILD}.zip"
    ),
    (ValidOpSys.WINDOWS, ValidArch.BIT_64): (
        f"https://s3.amazonaws.com/plink1-assets/"
        f"plink_win64_{PLINK_STABLE_BUILD}.zip"
    ),
    (ValidOpSys.WINDOWS, ValidArch.BIT_32): (
        f"https://s3.amazonaws.com/plink1-assets/"
        f"plink_win32_{PLINK_STABLE_BUILD}.zip"
    ),
    (ValidOpSys.DARWIN, ValidArch.BIT_64): (
        f"https://s3.amazonaws.com/plink1-assets/" f"plink_mac_{PLINK_STABLE_BUILD}.zip"
    ),
}


def get_arch(device: str) -> ValidArch:
    arch: dict[str, ValidArch] = {
        "x86_64": ValidArch.BIT_64,
        "amd64": ValidArch.BIT_64,
        "x64": ValidArch.BIT_64,
        "i386": ValidArch.BIT_32,
        "i486": ValidArch.BIT_32,
        "i586": ValidArch.BIT_32,
        "i686": ValidArch.BIT_32,
        "x86": ValidArch.BIT_32,
    }

    device = device.lower()
    try:
        return arch[device]
    except KeyError as exc:
        raise SystemExit(f"Unsupported device: {device}") from exc


def get_platform() -> tuple[str, tuple[ValidOpSys, ValidArch]]:
    op_sys: str = system().lower()
    if op_sys not in {os_.value for os_ in ValidOpSys}:
        raise SystemExit(f"Unsupported OS: {op_sys}")

    device: str = machine().lower()
    arch: ValidArch = get_arch(device=device)

    platform = (ValidOpSys(op_sys), arch)
    if platform not in PLINK_DOWNLOAD_URLS:
        raise SystemExit(f"Unsupported platform: {op_sys} {device}")
    return op_sys, platform


def get_exec_filepath(os_sys: str, dest_dir: str = "./plink") -> tuple[Path, str]:
    dest_dir: Path = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    exec_filename: str = "plink.exe" if os_sys == ValidOpSys.WINDOWS.value else "plink"
    return dest_dir, exec_filename


def install_plink() -> Path | None:
    os_sys, platform = get_platform()
    dest_dir, exec_filename = get_exec_filepath(os_sys=os_sys)
    plink_exec_path = Path(dest_dir / exec_filename)
    if plink_exec_path.exists():
        return plink_exec_path

    print("Installing PLINK...")
    with TemporaryDirectory[str](prefix="plink-tmp") as tmp_dir:
        url: str = PLINK_DOWNLOAD_URLS[platform]
        arch_path: Path = Path(tmp_dir) / "plink.zip"
        req: Request = Request(url=url)

        try:
            with urlopen(url=req, timeout=120) as res, arch_path.open(mode="wb") as f:
                copyfileobj(fsrc=res, fdst=f)
        except URLError as exc:
            raise SystemExit(f"Installing {exc} failed.") from exc

        with ZipFile(file=arch_path) as arch:
            is_member = [
                i for i in arch.namelist() if Path(i).name.lower() == exec_filename
            ]
            if is_member:
                member: ZipInfo = arch.getinfo(name=is_member[0])
                with arch.open(name=member) as src, plink_exec_path.open(
                    mode="wb"
                ) as f:
                    copyfileobj(fsrc=src, fdst=f)

                if name != "nt":
                    plink_exec_path.chmod(
                        mode=plink_exec_path.stat().st_mode
                        | S_IXUSR
                        | S_IXGRP
                        | S_IXOTH
                    )
                return plink_exec_path
    return None


if __name__ == "__main__":
    install_plink()
