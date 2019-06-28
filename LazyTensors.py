import numpy as np
from abc import *

########################################################################################################################
#
#                                           Tensor
#
class Tensor:
    """A semi-symbolic representation of a tensor. Allows for lazy processing of tensor operations."""
    def __init__(self, a, shape = None):
        self._a = a
        if shape is None:
            shape = a.shape
        self.shape = shape
    @classmethod
    def from_array(cls, a, shape = None):
        if isinstance(a, Tensor):
            return a
        else:
            cls(a, shape = shape)
    @property
    def array(self):
        if isinstance(self._a, np.ndarray):
            return self._a
        else:
            return self._a.array
    def get_dim(self):
        return len(self.shape)
    @property
    def dim(self):
        return self.get_dim()
    def __mul__(self, other):
        if isinstance(other, float):
            return TensorMul(self, other)
        else:
            return TensorDot(self, other)
    def __rmul__(self, other):
        if isinstance(other, float):
            return TensorMul(self, other)
        else:
            return TensorDot(self, other)
    def __add__(self, other):
        return TensorPlus(self, other)
    def __pow__(self, power, modulo=None):
        return TensorPow(self, power)

    def handle_missing_indices(self, missing, extant):
        # the default assumption is basically to use the diagonal element -- i.e. to assume that only the extant indices matter
        return self._a[extant]
    def pull_index(self, *idx):
        dim = len(idx)
        sdim = self.dim
        if dim > sdim:
            raise ValueError("{} has only {} dimensions but {} were requested".format(self, sdim, dim))
        else:
            real_shape = self._a.shape
            # here we add our hacks
            real_dim = len(real_shape)
            if real_dim < dim:
                a = self.handle_missing_indices(idx[:dim-real_dim], idx[dim-real_dim:])
            else:
                a = self._a[idx]
            try:
                new = type(self)(a)
            except ValueError:
                new = Tensor(a)
            return new

    def __getitem__(self, item):
        self.pull_index(*item)

    def __repr__(self):
        return "{}(<{}>)".format(type(self).__name__, ", ".join(*(str(x) for x in self.shape)))

########################################################################################################################
#
#                                           TensorOps
#
class TensorOp(metaclass=ABCMeta, Tensor):
    """A lazy representation of tensor operations to save memory"""
    def __init__(self, a, b):
        # if isinstance(a, Tensor):
        #     a = a._a
        # if isinstance(b, Tensor):
        #     a = b._a
        self._a = a
        self._b = b
    @abstractmethod
    def op(self, a, b):
        pass
    @abstractmethod
    def get_shape(self, a, b):
        return a.shape
    @property
    def shape(self):
        return self.get_shape(self._a, self._b)
    @property
    def array(self):
        """Ought to always compile down to a proper ndarray

        :return:
        :rtype: np.ndarray
        """
        return self.op(self._a, self._b)
    def __getitem__(self, i):
        return self.op(self._a[i], self._b)
class TensorPlus(TensorOp):
    """Represents an addition of two tensors"""
    def op(self, a, b):
        if isinstance(a, Tensor):
            a = a.array # with addition we can't really avoid it...
        if isinstance(b, Tensor):
            b = b.array
        return a + b
    def get_shape(self, a, b):
        return super().get_shape(a, b)
    def __getitem__(self, i):
        return type(self)(self._a[i], self._b[i])
class TensorPow(TensorOp):
    """Represents an addition of two tensors"""
    def op(self, a, b):
        if isinstance(a, Tensor):
            a = a.array
        return np.linalg.matrix_power(a, b)
    def get_shape(self, a, b):
        return super().get_shape(a, b)
    def __getitem__(self, i):
        return type(self)(self._a[i], self._b)
class TensorMul(TensorOp):
    """Represents a multiplication of a tensor and a scalar"""
    def op(self, a, b):
        if isinstance(a, Tensor):
            a = a.array # with scalar mult we can't really avoid it...
        return a * b
    def get_shape(self, a, b):
        return super().get_shape(a, b)
    def __getitem__(self, i):
        return type(self)(self._a[i], self._b)
class TensorDot(TensorOp):
    def get_shape(self, a, b):
        return a.shape[:-1] + b.shape[1:]
    def op(self, a, b):
        # we'll allow this to cast to dense as we'll assume people are only ever calling this on low dimensional tensors...
        if isinstance(a, Tensor):
            a = a.array
        if isinstance(b, Tensor):
            b = b.array
        return np.tensordot(a, b)
    def __getitem__(self, i):
        return type(self)(self._a[i], self._b)