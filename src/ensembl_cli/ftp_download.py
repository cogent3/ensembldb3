import os
import pathlib
import re

from ftplib import FTP
from typing import Callable, Iterable

from rich.progress import track
from unsync import unsync

from ensembl_cli.util import (
    atomic_write,
    checksum,
    load_ensembl_checksum,
    load_ensembl_md5sum,
    md5sum,
)


_sig_load_funcs = dict(CHECKSUMS=load_ensembl_checksum, MD5SUM=load_ensembl_md5sum)
_sig_calc_funcs = dict(CHECKSUMS=checksum, MD5SUM=md5sum)
_dont_write = re.compile("(CHECKSUMS|MD5SUM|README)")
_sig_file = re.compile("(CHECKSUMS|MD5SUM)")


def is_signature(path: os.PathLike) -> bool:
    return _sig_file.search((path.name)) is not None


def get_signature_data(path: os.PathLike) -> Callable:
    return _sig_load_funcs[path.name](path)


def configured_ftp(host: str = "ftp.ensembl.org") -> FTP:
    ftp = FTP(host)
    ftp.login()
    return ftp


def listdir(host: str, path: str, pattern: Callable = None):
    """returns directory listing"""
    pattern = pattern or (lambda x: True)
    ftp = configured_ftp(host=host)
    ftp.cwd(path)
    for fn in ftp.nlst():
        if pattern(fn):
            yield f"{path}/{fn}"
    ftp.close()


@unsync
def unsynced_copy_to_local(
    host: str, src: os.PathLike, dest: os.PathLike
) -> os.PathLike:
    if dest.exists():
        return dest
    ftp = configured_ftp(host=host)
    # pass in checksum and keep going until it's correct?
    with atomic_write(dest, mode="wb") as outfile:
        ftp.retrbinary(f"RETR {src}", outfile.write)

    ftp.close()
    return dest


def download_data(
    *,
    host: str,
    local_dest: os.PathLike,
    remote_paths: Iterable[os.PathLike],
    description,
) -> bool:
    tasks = [
        unsynced_copy_to_local(host, path, local_dest / pathlib.Path(path).name)
        for path in remote_paths
    ]
    saved_paths = [
        task.result() for task in track(tasks, description=description, transient=True)
    ]
    downloaded_chksums = {}
    checksums = {}
    calc_sig = None
    for path in saved_paths:
        if _dont_write.search(path.name):
            if is_signature(path):
                checksums = get_signature_data(path)
                calc_sig = _sig_calc_funcs[path.name]
            continue

        signature = calc_sig(path.read_bytes(), path.stat().st_size)
        downloaded_chksums[path.name] = signature

    for fn in downloaded_chksums:
        assert checksums[fn] == downloaded_chksums[fn], fn
    return True
