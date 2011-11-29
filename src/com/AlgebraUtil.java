package com;
import java.util.Arrays;

import javax.vecmath.GMatrix;

public class AlgebraUtil {
    
    public static double sigmoid(double z) {
        return 1. / (1. + Math.exp(-z));
    }

    public static double sumAllSquared(GMatrix mtx) {
        double sum = 0;
        for (int i = 0; i < mtx.getNumRow(); i++)
            for (int j = 0; j < mtx.getNumCol(); j++)
                sum += Math.pow(mtx.getElement(i, j), 2);
        return sum;
    }

    public static GMatrix transpose(GMatrix mtx) {
        GMatrix t = new GMatrix(mtx.getNumCol(), mtx.getNumRow());
        t.transpose(mtx);
        return t;
    }

    public static GMatrix scalarProduct(GMatrix mtx1, double value) {
        int n = mtx1.getNumCol();
        int m = mtx1.getNumRow();

        GMatrix z = new GMatrix(m, n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                z.setElement(i, j, mtx1.getElement(i, j) * value);
            }
        }
        return z;
    }

    public static GMatrix scalarDivide(GMatrix mtx1, double value) {
        int n = mtx1.getNumCol();
        int m = mtx1.getNumRow();

        GMatrix z = new GMatrix(m, n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                z.setElement(i, j, mtx1.getElement(i, j) / value);
            }
        }
        return z;
    }

    public static GMatrix product(GMatrix mtx1, GMatrix mtx2) {
        assert mtx1.getNumCol() == mtx2.getNumRow();
        GMatrix z = new GMatrix(mtx1.getNumRow(), mtx2.getNumCol());
        z.mul(mtx1, mtx2);
        return z;
    }

    public static GMatrix omitFirstColumn(GMatrix mtx) {
        int n = mtx.getNumCol() - 1;
        int m = mtx.getNumRow();
        GMatrix omit = new GMatrix(m, n);
        mtx.copySubMatrix(0, 1, m, n, 0, 0, omit);
        return omit;
    }

    public static GMatrix prependColumn(GMatrix z, double value) {

        int rows = z.getNumRow();

        GMatrix a = new GMatrix(rows, z.getNumCol() + 1);
        double[] column = new double[rows];
        Arrays.fill(column, value);
        a.setColumn(0, column);
        for (int i = 1; i < a.getNumCol(); i++) {
            z.getColumn(i - 1, column);
            a.setColumn(i, column);
        }
        return a;
    }

    interface ElementwiseFunction {
        public double computeElement(double el);
    }

    public static GMatrix apply(GMatrix z, ElementwiseFunction f) {
        GMatrix a = new GMatrix(z.getNumRow(), z.getNumCol());
        for (int i = 0; i < z.getNumRow(); i++) {
            for (int j = 0; j < z.getNumCol(); j++) {
                a.setElement(i, j, f.computeElement(z.getElement(i, j)));
            }
        }
        return a;
    }

    public static GMatrix dotProduct(GMatrix x, GMatrix y) {
        int n = x.getNumCol();
        int m = x.getNumRow();
        assert n == y.getNumCol();
        assert m == y.getNumRow();

        GMatrix z = new GMatrix(m, n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                z.setElement(i, j, x.getElement(i, j) * y.getElement(i, j));
            }
        }
        return z;
    }

    public static GMatrix add(GMatrix x, GMatrix y) {
        int n = x.getNumCol();
        int m = x.getNumRow();
        assert n == y.getNumCol();
        assert m == y.getNumRow();

        GMatrix z = new GMatrix(m, n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                z.setElement(i, j, x.getElement(i, j) + y.getElement(i, j));
            }
        }
        return z;
    }

    public static GMatrix substract(GMatrix x, GMatrix y) {
        int n = x.getNumCol();
        int m = x.getNumRow();
        assert n == y.getNumCol();
        assert m == y.getNumRow();

        GMatrix z = new GMatrix(m, n);
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                z.setElement(i, j, x.getElement(i, j) - y.getElement(i, j));
            }
        }
        return z;
    }

    public static GMatrix sumRows(GMatrix mtx) {
        int m = mtx.getNumRow();
        GMatrix z = new GMatrix(m, 1);
        for (int i = 0; i < m; i++) {
            double sum = 0;
            for (int j = 0; j < mtx.getNumCol(); j++) {
                sum += mtx.getElement(i, j);
            }
            z.setElement(i, 0, sum);
        }
        return z;
    }

    public static final ElementwiseFunction mtxSigmoid = new ElementwiseFunction() {
        @Override
        public double computeElement(double el) {
            return sigmoid(el);
        }
    };
    public static final ElementwiseFunction mtxLog = new ElementwiseFunction() {
        @Override
        public double computeElement(double el) {
            return Math.log(el);
        }
    };
    public static final ElementwiseFunction mtxNeg = new ElementwiseFunction() {
        @Override
        public double computeElement(double el) {
            return el * -1;
        }
    };
    
    public static final ElementwiseFunction mtxMinusOne = new ElementwiseFunction() {
        @Override
        public double computeElement(double el) {
            return el - 1;
        }
    };


    public static final ElementwiseFunction mtxOneMinus = new ElementwiseFunction() {
        @Override
        public double computeElement(double el) {
            return 1 - el;
        }
    };

    public static final ElementwiseFunction mtxSigmoidGradient = new ElementwiseFunction() {
        @Override
        public double computeElement(double el) {
            return sigmoid(el) * (1 - sigmoid(el));
        }
    };

}
