import copy

from dislib.data.array import Array


class Preprocessing:
    def __init__(self, X=None, y=None, method=None, data_analytics=None):
        if X is not None:
            if not isinstance(X, Array):
                raise TypeError("X data should be of class ds-array")
            else:
                self.X = X
        if y is not None:
            if not isinstance(y, Array):
                raise TypeError("X data should be of class ds-array")
            else:
                self.y = y
        if method is not None:
            if isinstance(method, list):
                self.methods = method
            else:
                self.methods = [method]
        else:
            self.methods = []
        self.methods_x = None
        self.methods_y = None
        if data_analytics is not None:
            if isinstance(data_analytics, list):
                self.data_analytics = data_analytics
            else:
                self.data_analytics = [data_analytics]
        else:
            self.data_analytics = []

    def add_method(self, estimator):
        module_tree = getattr(estimator, '__module__', None)
        parent = module_tree.split('.')[0] if module_tree else None
        if parent == sklearn.__name__ or parent == dislib.__name__:
            self.methods.append(estimator)
        else:
            pass#TODO identificar si tiene fit y transform metodos o dar una excepción

    def fit(self, X, y=None):
        if len(self.methods) == 0:
            raise ValueError("There are no methods to fit with the data")
        if y is not None:
            self.methods_x = copy.deepcopy(self.methods)
            for method in self.methods_x:
                method.fit(X)
            self.methods_y = copy.deepcopy(self.methods)
            for method in self.methods_y:
                method.fit(y)
        else:
            self.methods_x = copy.deepcopy(self.methods)
            for method in self.methods_x:
                method.fit(X)

    def fit_transform(self, X, y=None):
        self.fit(X, y=y)
        return self.transform(X, y=y)

    def transform(self, X, y=None):
        if len(self.methods) == 0:
            raise ValueError("There are no methods to apply on data")
        if self.methods_x is not None:
            for method in self.methods_x:
                X = method.transform(X)
        else:
            raise ValueError("There are not fitted methods to transform the data")
        if y is not None:
            if self.methods_y is not None:
                for method in self.methods_y:
                    y = method.transform(y)
            else:
                raise ValueError("There are not fitted methods to transform the target data")
        return X, y

    def save_model(self, key=None):
        if key is not None:
            pass
        else:
            for method in self.methods:
                pass

    def load_model(self, path):
        self.models.append(path)#TODO CAMBIAR ESTO POR EL CORRESPONDIENTE LOAD

    def read(self, path):
        pass

    def write(self, path):
        pass

    def set_data(self, X, y=None):
        if not isinstance(X, Array) and not isinstance(y, Array):#Todo Falta el tensor
            raise ValueError("Data should be contained on a distributed stored object like:"
                             "ds-array or ds-Tensor")
        self.X = X
        self.y = y
