package com.smallnn;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Arrays;

import javax.vecmath.GMatrix;

import org.junit.Test;


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
        t.setZero();
        t.transpose(mtx);
        return t;
    }

    public static GMatrix scalarProduct(GMatrix mtx1, double value) {
        int n = mtx1.getNumCol();
        int m = mtx1.getNumRow();

        GMatrix z = new GMatrix(m, n); z.setZero();
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

        GMatrix z = new GMatrix(m, n); z.setZero();
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                z.setElement(i, j, mtx1.getElement(i, j) / value);
            }
        }
        return z;
    }

    public static GMatrix scalarVectorDivide(GMatrix mtx, GMatrix vec){
        int n = mtx.getNumCol();
        int m = mtx.getNumRow();
        GMatrix sub = new GMatrix(m,n); sub.setZero();
        for(int i=0;i<m;i++)
            for (int j=0;j<n;j++) {
                double element = vec.getElement(0, j);
                double element2 = mtx.getElement(i, j);
                if (element ==0 && element2==0){
                    sub.setElement(i, j, 0.);
                } else {
                    sub.setElement(i, j, element2/element);
                }
            }
        return sub;
    }
    
    public static GMatrix product(GMatrix mtx1, GMatrix mtx2) {
        assert mtx1.getNumCol() == mtx2.getNumRow();
        GMatrix z = new GMatrix(mtx1.getNumRow(), mtx2.getNumCol()); z.setZero();
        z.mul(mtx1, mtx2);
        return z;
    }

    public static GMatrix omitFirstColumn(GMatrix mtx) {
        int n = mtx.getNumCol() - 1;
        int m = mtx.getNumRow();
        GMatrix omit = new GMatrix(m, n); omit.setZero();
        mtx.copySubMatrix(0, 1, m, n, 0, 0, omit);
        return omit;
    }

    public static GMatrix prependColumn(GMatrix z, double value) {

        int rows = z.getNumRow();

        GMatrix a = new GMatrix(rows, z.getNumCol() + 1); a.setZero();
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
        GMatrix a = new GMatrix(z.getNumRow(), z.getNumCol()); a.setZero();
        for (int i = 0; i < z.getNumRow(); i++) {
            for (int j = 0; j < z.getNumCol(); j++) {
                a.setElement(i, j, f.computeElement(z.getElement(i, j)));
            }
        }
        return a;
    }
    
    public static GMatrix vectorize(GMatrix mtx){
        int m = mtx.getNumRow();
        int n = mtx.getNumCol();
        GMatrix z = new GMatrix(m, n); z.setZero();
        for (int i=0; i<m; i++){
            double[] row = new double[n];
            mtx.getRow(i, row);
            int maxidx = 0; double maxval = row[maxidx];
            for (int j=0;j<n;j++){
                if (row[j]>maxval){
                    maxval = row[j];
                    maxidx = j;
                }
            }
            double[] vectorized = new double[n];
            Arrays.fill(vectorized, 0.);
            vectorized[maxidx] = 1;
            z.setRow(i, vectorized);
        }
        return z;
    }

    public static GMatrix dotProduct(GMatrix x, GMatrix y) {
        int n = x.getNumCol();
        int m = x.getNumRow();
        assert n == y.getNumCol();
        assert m == y.getNumRow();

        GMatrix z = new GMatrix(m, n); z.setZero();
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

        GMatrix z = new GMatrix(m, n); z.setZero();
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

        GMatrix z = new GMatrix(m, n); z.setZero();
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                z.setElement(i, j, x.getElement(i, j) - y.getElement(i, j));
            }
        }
        return z;
    }

    public static GMatrix sumRows(GMatrix mtx) {
        int m = mtx.getNumRow();
        GMatrix z = new GMatrix(m, 1); z.setZero();
        for (int i = 0; i < m; i++) {
            double sum = 0;
            for (int j = 0; j < mtx.getNumCol(); j++) {
                sum += mtx.getElement(i, j);
            }
            z.setElement(i, 0, sum);
        }
        return z;
    }
    
    public static GMatrix sumColumns(GMatrix mtx){
        int n = mtx.getNumCol();
        GMatrix z = new GMatrix(1, n);  z.setZero();
        
        for (int i =0; i<mtx.getNumRow();i++){
            double[] row = new double[n];
            mtx.getRow(i, row);
            GMatrix rowMtx = new GMatrix(1,n,row);
            z.add(rowMtx);
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
    
    public static final ElementwiseFunction mtxSqrt = new ElementwiseFunction() {
        @Override
        public double computeElement(double el) {
            return Math.sqrt(el);
        }
    };
    
    @Test
    public void testMean(){
        GMatrix mtx = new GMatrix(3,2, new double[]{1,2,3,4,5,6});
        GMatrix mean = mean(mtx);
        assertTrue(3. == mean.getElement(0, 0));
        assertTrue(4. == mean.getElement(0, 1));
        
        GMatrix mtx2 = new GMatrix(2,2, new double[]{0,0,0,0});
        GMatrix mean2 = mean(mtx2);
        assertEquals(0., mean2.getElement(0, 0), 0.00001);
        assertEquals(0., mean2.getElement(0, 1), 0.00001);
    }
    
    public static GMatrix mean(GMatrix mtx){
        GMatrix colSums = sumColumns(mtx);
        return scalarDivide(colSums, mtx.getNumRow());
    }
    
    @Test
    public void testStd(){
        GMatrix mtx = new GMatrix(3,2, new double[]{1,2,3,4,5,6});
        GMatrix std = std(mtx);
        
        assertTrue(2. == std.getElement(0, 0));
        assertTrue(2. == std.getElement(0, 1));
        
        GMatrix mtx2 = new GMatrix(2,2, new double[]{0,0,0,0});
        GMatrix std2 = std(mtx2);
        
        assertEquals(0., std2.getElement(0, 0), 1e-9);
        assertEquals(0., std2.getElement(0, 1), 1e-9);
    }
    
    public static GMatrix std(GMatrix mtx){
        //sqrt(sumsq(x-mean(x))/(length(x)-1))
        int m = mtx.getNumRow();
        GMatrix sub = substractVector(mtx, mean(mtx));
        GMatrix sumSQ = sumSQColumns(sub);
        GMatrix divided = scalarDivide(sumSQ, m-1.);
        return apply(divided, mtxSqrt);
    }

    public static GMatrix substractVector(GMatrix mtx, GMatrix mean) {
        int n = mtx.getNumCol();
        int m = mtx.getNumRow();
        GMatrix sub = new GMatrix(m,n); sub.setZero();
        for(int i=0;i<m;i++)
            for (int j=0;j<n;j++)
                sub.setElement(i, j, mtx.getElement(i, j)-mean.getElement(0, j));
        return sub;
    }
    
    private static GMatrix sumSQColumns(GMatrix mtx){
        int n = mtx.getNumCol();
        GMatrix z = new GMatrix(1, n); z.setZero();
        for(int i=0;i<mtx.getNumRow();i++){
            double[] rowSQ = new double[n];
            mtx.getRow(i, rowSQ);
            for (int j=0;j<n;j++){
                rowSQ[j] = Math.pow(rowSQ[j], 2);
            }
            z.add(new GMatrix(1,n,rowSQ));
        }
        return z;
    }
    
    public static GMatrix[] computeNormalizationParams(GMatrix X) {
        /* mu = mean(X); */
        GMatrix mu = mean(X);
        /* X_norm = X - kron(mu, ones(m, 1)); */
        GMatrix X_mean = substractVector(X, mu);

        /* sigma = std(X_norm) = sqrt(sumsq(x-mean(x))/(length(x)-1)); */
        GMatrix sigma = std(X_mean);
        /* X_norm = X_norm ./ kron(sigma, ones(m, 1)); */
        GMatrix X_norm = scalarVectorDivide(X_mean, sigma);
        return new GMatrix[]{mu, sigma, X_norm};
    }
    
    public static GMatrix normalize(GMatrix mu, GMatrix sigma, GMatrix x){
        GMatrix xMinusMean = substractVector(x, mu);
        return scalarVectorDivide(xMinusMean, sigma);
    }
}
