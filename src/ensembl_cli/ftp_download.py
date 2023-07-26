import os
import pathlib
import re
import shutil
import uuid

from ftplib import FTP
from tempfile import mkdtemp
from typing import IO, Callable, Iterable

from rich.progress import track
from unsync import unsync

from ensembl_cli.util import checksum, load_ensembl_checksum


dont_write = re.compile("(CHECKSUMS|README)")


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


class atomic_write:
    """performs atomic write operations, cleans up if fails"""

    def __init__(self, path: os.PathLike, tmpdir=None, mode="wb", encoding=None):
        """

        Parameters
        ----------
        path
            path to file
        tmpdir
            directory where temporary file will be created
        mode
            file writing mode
        encoding
            text encoding
        """
        path = pathlib.Path(path).expanduser()

        self._path = path
        self._mode = mode
        self._file = None
        self._encoding = encoding
        self._tmppath = self._make_tmppath(tmpdir)

        self.succeeded = None
        self._close_func = self._close_rename_standard

    def _make_tmppath(self, tmpdir):
        """returns path of temporary file

        Parameters
        ----------
        tmpdir: Path
            to directory

        Returns
        -------
        full path to a temporary file

        Notes
        -----
        Uses a random uuid as the file name, adds suffixes from path
        """
        suffixes = "".join(self._path.suffixes)
        parent = self._path.parent
        name = f"{uuid.uuid4()}{suffixes}"
        tmpdir = (
            pathlib.Path(mkdtemp(dir=parent))
            if tmpdir is None
            else pathlib.Path(tmpdir)
        )

        if not tmpdir.exists():
            raise FileNotFoundError(f"{tmpdir} directory does not exist")

        return tmpdir / name

    def _get_fileobj(self):
        """returns file to be written to"""
        if self._file is None:
            self._file = open(self._tmppath, self._mode)

        return self._file

    def __enter__(self) -> IO:
        return self._get_fileobj()

    def _close_rename_standard(self, src):
        dest = pathlib.Path(self._path)
        try:
            dest.unlink()
        except FileNotFoundError:
            pass
        finally:
            src.rename(dest)

        shutil.rmtree(src.parent)

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._file.close()
        if exc_type is None:
            self._close_func(self._tmppath)
            self.succeeded = True
        else:
            self.succeeded = False
            shutil.rmtree(self._tmppath.parent)

    def write(self, text):
        """writes text to file"""
        fileobj = self._get_fileobj()
        fileobj.write(text)

    def close(self):
        """closes file"""
        self.__exit__(None, None, None)


@unsync
def unsynced_copy_to_local(host, src, dest):
    #  TODO check if path exists and satisfies chksum
    # return when both conditions satisfied
    ftp = configured_ftp(host=host)
    # pass in checksum and keep going until it's correct?
    with atomic_write(dest, mode="wb") as outfile:
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

        summed, blocks = checksum(path.read_bytes(), path.stat().st_size)
        checkpoint_file.write(f"{summed}\t{blocks}\t{path}\n")
        downloaded_chksums[path.name] = summed, blocks

    for fn in downloaded_chksums:
        assert checksums[fn] == downloaded_chksums[fn], fn
    return True
