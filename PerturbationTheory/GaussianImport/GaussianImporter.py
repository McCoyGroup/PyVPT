"""Implements an importer from a Gaussian log file

"""

import numpy as np
from mmap import mmap

class FileCheckPoint:
    def __init__(self, parent):
        self._parent = parent
        self._chk = parent.tell
    def __enter__(self, ):
        return self
    def __exit__(self, exc_type, exc_val, exc_tb):
        self._parent.seek(self._chk)

class GaussianLogReader:
    """Implements a stream based reader for a Gaussian file... a bit messy

    """
    def __init__(self, file, mode="r", encoding="utf-8", **kw):
        self._file = file
        self._mode = mode
        self._encoding = encoding
        self._kw = kw
        self._stream = None

    def __enter__(self):
        self._fstream = open(self._file, mode=self._mode, encoding=self._encoding, **self._kw)
        self._stream = mmap(self._fstream.fileno(), 0)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._stream.close()
        self._fstream.close()

    def __iter__(self):
        return iter(self._fstream)

    def __getattr__(self, item):
        return getattr(self._stream, item)

    def readline(self):
        return self._fstream.readline()
    def read(self, n=1):
        return self._fstream.read(n)

    def seek_tag(self, tag):
        """

        :param header: a string specifying a header to look for
        :type header: str
        :return: if header was found
        :rtype: bool
        """
        pos = self._stream.find(tag.encode(self._encoding))
        if pos > 0:
            self._stream.seek(pos)
        return pos
    def get_tagged_block(self, tag_start, tag_end, block_size = 500):
        """Pulls the string between tag_start and tag_end

        :param tag_start:
        :type tag_start:
        :param tag_end:
        :type tag_end:
        :return:
        :rtype:
        """
        start = self.seek_tag(tag_start)
        if start > 0:
            with FileCheckPoint(self):
                end = self.seek_tag(tag_end)
            if end > 0:
                return self._stream.read(start-end).decode(self._encoding)

        # implict None return if no block found




