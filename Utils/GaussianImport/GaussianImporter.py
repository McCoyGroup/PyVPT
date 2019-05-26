"""Implements an importer from a Gaussian log file

"""

import numpy as np, re, math, io
from mmap import mmap
from .GaussianLogComponents import GaussianLogComponents, GaussianLogDefaults, GaussianLogOrdering
from .GaussianFChkComponents import FormattedCheckpointComponents, FormattedCheckpointCommonNames

########################################################################################################################
#
#                                           FileStreamReader
#
class FileCheckPoint:
    def __init__(self, parent):
        self._parent = parent
        self._chk = parent.tell()
    def __enter__(self, ):
        return self
    def __exit__(self, exc_type, exc_val, exc_tb):
        self._parent.seek(self._chk)

class FileStreamReaderException(IOError):
    pass

class FileStreamReader:
    def __init__(self, file, mode="r", encoding="utf-8", **kw):
        self._file = file
        self._mode = mode.strip("+b")+"+b"
        self._encoding = encoding
        self._kw = kw
        self._stream = None

    def __enter__(self):
        self._fstream = open(self._file, mode=self._mode, **self._kw)
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
        return self._stream.readline()
    def read(self, n=1):
        return self._stream.read(n)
    def seek(self, *args, **kwargs):
        # self._fstream.seek(*args, **kwargs)
        return self._stream.seek(*args, **kwargs)
    def tell(self):
        return self._stream.tell()

    def find_tag(self, tag, skip_tag = True, seek = True):
        """

        :param header: a string specifying a header to look for
        :type header: str
        :return: if header was found
        :rtype: bool
        """
        enc_tag = tag.encode(self._encoding)
        pos = self._stream.find(enc_tag)
        if seek and pos > 0:
            if skip_tag:
                pos = pos + len(enc_tag)
            self._stream.seek(pos)
        return pos
    def get_tagged_block(self, tag_start, tag_end, block_size = 500):
        """Pulls the string between tag_start and tag_end

        :param tag_start:
        :type tag_start: str or None
        :param tag_end:
        :type tag_end: str
        :return:
        :rtype:
        """
        if tag_start is not None:
            start = self.find_tag(tag_start)
            if start > 0:
                with FileCheckPoint(self):
                    end = self.find_tag(tag_end, seek=False)
                if end > 0:
                    return self._stream.read(end-start).decode(self._encoding)
        else:
            start = self.tell()
            with FileCheckPoint(self):
                end = self.find_tag(tag_end, seek=False)
            if end > 0:
                return self._stream.read(end-start).decode(self._encoding)

        # implict None return if no block found

    def parse_key_block(self, tag_start=None, tag_end=None, mode=None, parser = None, num = None, **ignore):
        """

        :param key: registered key pattern to pull from Gaussian
        :type key: str
        :return:
        :rtype:
        """
        if tag_start is None:
            raise FileStreamReaderException("{}.{}: needs a '{}' argument".format(
                type(self).__name__,
                "parse_key_block",
                "tag_start"
            ))
        if tag_end is None:
            raise FileStreamReaderException("{}.{}: needs a '{}' argument".format(
                type(self).__name__,
                "parse_key_block",
                "tag_end"
            ))
        if mode == "List":
            with FileCheckPoint(self):
                # we do this in a checkpointed fashion only for list-type tokens
                # for all other tokens we introduce an ordering to apply when checking
                # does it need to be done like this... probably not?
                # I suppose we could be a significantly more efficient by returning a
                # generator statement in these multi-block cases
                # and then introducing a sorting order across these multi-blocks
                # that tells us which to check first, second, etc.
                # but this is probably good enough
                if isinstance(num, int):
                    blocks = [None]*num
                    for i in range(num):
                        block = self.get_tagged_block(tag_start, tag_end)
                        if block is None:
                            break
                        blocks[i] = block
                    if parser is None:
                        parser = lambda a:a
                    blocks = parser(blocks[:i+1])
                else:
                    blocks = []
                    block = self.get_tagged_block(tag_start, tag_end)
                    while block is not None:
                        blocks.append(block)
                        block = self.get_tagged_block(tag_start, tag_end)
                    if parser is None:
                        parser = lambda a:a
                    blocks = parser(blocks)
                return blocks
        else:
            block = self.get_tagged_block(tag_start, tag_end)
            if parser is not None:
                block = parser(block)
            return block

########################################################################################################################
#
#                                           GaussianLogReader
#
class GaussianLogReader(FileStreamReader):
    """Implements a stream based reader for a Gaussian .log file... a bit messy

    """

    registered_components = GaussianLogComponents
    default_keys = GaussianLogDefaults
    default_ordering = GaussianLogOrdering

    def parse(self, keys = None, num = None):
        """The main function we'll actually use. Parses bits out of a .log file.

        :param keys: the keys we'd like to read from the log file
        :type keys: str or list(str)
        :param num: for keys with multiple entries, the number of entries to pull
        :type num: int or None
        :return:
        :rtype:
        """
        if keys is None:
            keys = self.default_keys
        # important for ensuring correctness of what we pull
        if isinstance(keys, str):
            keys = (keys,)
        keys = sorted(keys, key = lambda k:self.default_ordering[k])
        res = { k:self.parse_key_block(**self.registered_components[k], num=num) for k in keys }
        return res

########################################################################################################################
#
#                                           GaussianFChkReader
#
class GaussianFChkReaderException(FileStreamReaderException):
    pass

class GaussianFChkReader(FileStreamReader):
    """Implements a stream based reader for a Gaussian .fchk file. Pretty generall I think. Should be robust-ish.
    One place to change things up is convenient parsers for specific commonly pulled parts of the fchk

    """

    registered_components = FormattedCheckpointComponents
    common_names = {to_:from_ for from_, to_ in FormattedCheckpointCommonNames.items()}
    to_common_name = FormattedCheckpointCommonNames
    def read_header(self):
        """Reads the header and skips the stream to where we want to be

        :return: the header
        :rtype: str
        """
        return self.get_tagged_block(None, "Number of atoms")

    fchk_re_pattern = r"^(.+?)\s+(I|R)\s+(N=)?\s+(.+)\s+" # matches name, type, num (if there), and val
    fchk_re = re.compile(fchk_re_pattern)
    def get_next_block_params(self):
        """Pulls the tag of the next block, the type, the number of bytes it'll be,
        and if it's a single-line block it'll also spit back the block itself

        :return:
        :rtype: dict
        """
        with FileCheckPoint(self):
            tag_line = self.readline()
            if tag_line is b'':
                return None
            tag_line = tag_line.decode(self._encoding)
            match = re.match(self.fchk_re, tag_line)
            if match is None:
                raise FileStreamReaderException("{}.{}: line '{}' couldn't be read as a tag line".format(
                    type(self).__name__,
                    "get_next_block_params",
                    tag_line
                ))
            jump = self.tell()
        self.seek(jump)

        name, btype, numQ, val = gg = match.groups()
        # print(gg)
        if numQ:
            byte_count = 0
            shits = int(val)
            # hard coded these block formats since they're documented by Gaussian and thus unlikely to change
            if btype == "I":
                byte_per_shit = 12
                shits_per_line = 6
                btype = int
            elif btype == "R":
                byte_per_shit = 16
                shits_per_line = 5
                btype = float
            elif btype == "C":
                byte_per_shit = 12
                shits_per_line = 5
                btype = str
            elif btype == "L":
                byte_per_shit = 12
                shits_per_line = 5
                btype = bool
            byte_count = shits * byte_per_shit + math.ceil( shits / shits_per_line ) # each newline needs a byte
        else:
            byte_count = None
            if btype == "I":
                val = int(val)
            elif btype == "R":
                val = float(val)
            elif btype == "L":
                val = bool(int(val))

        return {
            "name": name,
            "dtype": btype,
            "byte_count": byte_count,
            "value": val
        }

    def get_block(self, name = None, dtype = None, byte_count = None, value = None):
        """Pulls the next block by first pulling the block tag

        :return:
        :rtype:
        """

        if byte_count is not None:

            block_str = self.read(byte_count)
            if dtype == int:
                value = np.fromstring(block_str, np.int8)
            elif dtype == float:
                block_str = io.StringIO(block_str.decode(self._encoding).replace("\n", "")) # flatten it out
                value = np.loadtxt(block_str)
            else:
                value = block_str

        try:
            parser = self.registered_components[name]
            value = parser(value)
        except KeyError:
            pass

        return value

    def skip_block(self, name = None, dtype = None, byte_count = None, value = None):
        """Skips the next block

        :return:
        :rtype:
        """

        if byte_count is not None:
            self.seek(self.tell() + byte_count)

    def parse(self, keys = None):
        if keys is None:
            keys_to_go = None
        else:
            if isinstance(keys, str):
                keys = (keys,)
            keys_original = set(keys)
            keys_to_go = { (self.common_names[k] if k in self.common_names else k) for k in keys }

        parse_results = {}
        header = self.read_header()

        parse_results['Header'] = header
        if keys_to_go is None:
            while True: # I'll just break once I've exhausted everything
                next_block = self.get_next_block_params()
                if next_block is None:
                    break
                tag = next_block["name"]
                parse_results[tag] = self.get_block(**next_block)
        else:
            while len(keys_to_go)>0:
                next_block = self.get_next_block_params()
                if next_block is None:
                    raise GaussianFChkReaderException("{}.{}: couldn't find keys {}".format(
                        type(self).__name__,
                        "parse",
                        keys_to_go
                        )
                    )
                tag = next_block["name"]
                if tag in keys_to_go:
                    keys_to_go.remove(tag)
                    if tag not in keys_original:
                        tag = self.to_common_name[tag]
                    parse_results[tag] = self.get_block(**next_block)
                else:
                    self.skip_block(**next_block)

        return parse_results


