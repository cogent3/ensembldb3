import os
import pathlib
import re

from ftplib import FTP
from typing import Iterable

from rich.progress import track
from unsync import unsync

from ensembl_cli.util import checksum, load_ensembl_checksum


dont_write = re.compile("(CHECKSUMS|README)")


def configured_ftp(host: str = "ftp.ensembl.org") -> FTP:
    ftp = FTP(host)
    ftp.login()
    return ftp


def listdir(host: str, path: str, debug=True):
    """returns directory listing"""
    ftp = configured_ftp(host=host)
    ftp.cwd(path)
    yield from ftp.nlst()
    ftp.close()


@unsync
def unsynced_copy_to_local(host, src, dest):
    #  TODO check if path exists and satisfies chksum
    # return when both conditions satisfied
    ftp = configured_ftp(host=host)
    # pass in checksum and keep going until it's correct?
    with open(dest, mode="wb") as outfile:
        ftp.retrbinary(f"RETR {src}", outfile.write)

    ftp.close()
    return dest


def download_data(
    host: str,
    local_dest: os.PathLike,
    remote_paths: Iterable[os.PathLike],
    description,
    checkpoint_file,
) -> bool:
    tasks = [
        unsynced_copy_to_local(host, path, local_dest / pathlib.Path(path).name)
        for path in remote_paths
    ]
    saved_paths = [
        task.result() for task in track(tasks, description=description, transient=True)
    ]
    downloaded_chksums = {}
    for path in saved_paths:
        if dont_write.search(path.name):
            if path.name == "CHECKSUMS":
                checksums = load_ensembl_checksum(path)
            continue
        summed, blocks = checksum(path)
        checkpoint_file.write(f"{summed}\t{blocks}\t{path}\n")
        downloaded_chksums[path.name] = summed, blocks
    for fn in checksums:
        assert checksums[fn] == downloaded_chksums[fn], fn
    return True
