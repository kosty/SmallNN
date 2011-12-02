package com;

import static com.AlgebraUtil.apply;
import static com.AlgebraUtil.dotProduct;
import static com.AlgebraUtil.mean;
import static com.AlgebraUtil.mtxLog;
import static com.AlgebraUtil.mtxNeg;
import static com.AlgebraUtil.mtxOneMinus;
import static com.AlgebraUtil.mtxSigmoid;
import static com.AlgebraUtil.mtxSigmoidGradient;
import static com.AlgebraUtil.omitFirstColumn;
import static com.AlgebraUtil.prependColumn;
import static com.AlgebraUtil.product;
import static com.AlgebraUtil.scalarProduct;
import static com.AlgebraUtil.scalarVectorDivide;
import static com.AlgebraUtil.std;
import static com.AlgebraUtil.substract;
import static com.AlgebraUtil.substractVector;
import static com.AlgebraUtil.sumAllSquared;
import static com.AlgebraUtil.sumRows;
import static com.AlgebraUtil.transpose;
import static com.AlgebraUtil.vectorize;
import static com.ChartUtil.showChart;
import static com.FileUtil.filenames;
import static com.FileUtil.testcases;
import static com.ImageUtil.HEIGHT;
import static com.ImageUtil.SMALL_HEIGHT;
import static com.ImageUtil.SMALL_WIDTH;
import static com.ImageUtil.WIDTH;
import static com.ImageUtil.imageToGrayscale;
import static com.ImageUtil.readImage;
import static com.ImageUtil.resize;
import static com.TrainDataUtil.TEST_DIR;
import static com.TrainDataUtil.getMixedData;
import static com.TrainDataUtil.listSetFiles;
import static com.TrainDataUtil.readData;

import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import javax.imageio.ImageIO;
import javax.vecmath.GMatrix;

import org.junit.Before;

import com.TrainDataUtil.Data;

public class TrainNN {
    private static final int NUM_TRIES = 1000;
    Random r = new Random();

    int features = 60;
    int classes = 2;

    public static final int IMAGE_SIZE = WIDTH * HEIGHT;
    public static final int SMALL_IMAGE_SIZE = SMALL_WIDTH * SMALL_HEIGHT;

    int m = filenames.length; 
    int inputLayerSize = IMAGE_SIZE;
    double[] x = new double[m * inputLayerSize];
    double[] y = new double[m * classes];

    int m_test = testcases.length;

    public void randomInit() {
        m = 50;
        inputLayerSize = 2;
        features = 600;
        x = new double[m * inputLayerSize];
        y = new double[m * classes];
        for (int i = 0; i < m; i++) {
            double acc = 0.;
            for (int j = 0; j < inputLayerSize; j++) {
                double val = -1. + (Math.random() * 2.);
                x[i * inputLayerSize + j] = val;
                acc += val;
            }
            if (acc <= 0) {
                y[i * classes] = 1.;
                y[i * classes + 1] = 0.;
            } else {
                y[i * classes] = 0.;
                y[i * classes + 1] = 1.;
            }
        }

    }
    
    public void blackWhiteInit(){
        double[] white = new double[]{16777215., 16777215., 16777215., 16777215., 16777215., 16777215., 16777215., 16777215.}; 
        double[] black = new double[]{0., 0., 0., 0., 0., 0., 0., 0.}; 
        m = 50;
        inputLayerSize = white.length;
        features = 600;
        x = new double[m * inputLayerSize];
        y = new double[m * classes];
        for(int i=0;i<m;i++){
            double[] tmparr = i%2 == 0 ? white : black; 
            for (int j=0;j<inputLayerSize;j++){
                x[i*inputLayerSize+j] = tmparr[j];
            }
            y[i*classes]   = i%2 == 0 ? 1. : 0.;
            y[i*classes+1] = i%2 == 0 ? 0. : 1.;
        }
        
    }
    
    public static double[] longToDoubleArray(long[] arr){
        double[] result = new double[arr.length];
        for(int i=0;i<arr.length;i++){
            result[i] = (double)arr[i];
        }
        return result;
    }
    
    public void controlledImageSetInit() throws Exception{
        double[] bars0 = longToDoubleArray(imageToGrayscale(resize(readImage(filenames[0]), SMALL_HEIGHT, SMALL_WIDTH))); //bars
        double[] nonBars1 =longToDoubleArray(imageToGrayscale(resize(readImage(filenames[1]), SMALL_HEIGHT, SMALL_WIDTH))); //non-bars
        
        double[] bars2 = longToDoubleArray(imageToGrayscale(resize(readImage(filenames[2]), SMALL_HEIGHT, SMALL_WIDTH))); //bars
        double[] nonBars3 =longToDoubleArray(imageToGrayscale(resize(readImage(filenames[3]), SMALL_HEIGHT, SMALL_WIDTH))); //non-bars
        
        double[] bars4 = longToDoubleArray(imageToGrayscale(resize(readImage(filenames[4]), SMALL_HEIGHT, SMALL_WIDTH))); //bars
        double[] nonBars5 =longToDoubleArray(imageToGrayscale(resize(readImage(filenames[5]), SMALL_HEIGHT, SMALL_WIDTH))); //non-bars
        
        double[] bars6 = longToDoubleArray(imageToGrayscale(resize(readImage(filenames[6]), SMALL_HEIGHT, SMALL_WIDTH))); //bars
        double[] nonBars7 =longToDoubleArray(imageToGrayscale(resize(readImage(filenames[7]), SMALL_HEIGHT, SMALL_WIDTH))); //non-bars
        
        double[] bars8 = longToDoubleArray(imageToGrayscale(resize(readImage(filenames[8]), SMALL_HEIGHT, SMALL_WIDTH))); //bars
        double[] nonBars9 =longToDoubleArray(imageToGrayscale(resize(readImage(filenames[9]), SMALL_HEIGHT, SMALL_WIDTH))); //non-bars
        
        double[][] imagesArray = new double[][]{bars0, nonBars1, bars2, nonBars3, bars4, nonBars5, bars6, nonBars7, bars8, nonBars9};
        inputLayerSize = bars0.length;
        m=imagesArray.length*40;
        x = new double[m*inputLayerSize];
        y = new double[m*classes];
        features = 600;
        x = new double[m * inputLayerSize];
        y = new double[m * classes];
        for(int i=0;i<m;i++){
            double[] tmparr = imagesArray[i%imagesArray.length]; 
            for (int j=0;j<inputLayerSize;j++){
                x[i*inputLayerSize+j] = tmparr[j];
            }
            y[i*classes]   = i%2 == 0 ? 1. : 0.;
            y[i*classes+1] = i%2 == 0 ? 0. : 1.;
        }
    }

    @Before
    public void smallImageInit() throws Exception {
        inputLayerSize = SMALL_IMAGE_SIZE;
        x = new double[m * inputLayerSize];
        y = new double[m * classes];

        for (int l = 0; l < m; l++) {
            BufferedImage img = resize(ImageIO.read(new File(filenames[l])), SMALL_WIDTH, SMALL_HEIGHT);
            // System.out.println(filenames[l]);
            int height = img.getHeight();
            assert height == 5;
            for (int i = 0; i < height; i++) {
                int width = img.getWidth();
                assert width == 6;
                for (int j = 0; j < width; j++) {
                    int pixelInx = l * inputLayerSize + i * height + j;
                    x[pixelInx] = img.getRGB(j, i);
                }
            }
            // System.out.println();
            int idx = l * classes;
            if (filenames[l].contains("non-bars")) {
                y[idx] = 0;
                y[idx + 1] = 1;
            } else {
                y[idx] = 1;
                y[idx + 1] = 0;
            }
        }
    }

    public void fullImageInit() throws Exception {
        for (int l = 0; l < filenames.length; l++) {
            StringBuilder sb = new StringBuilder();
            int cnt = 0;

            BufferedImage img = ImageIO.read(new File(filenames[l]));
            sb.append(filenames[l]);
            int height = img.getHeight();
            assert height == HEIGHT;
            for (int i = 0; i < height; i++) {
                int width = img.getWidth();
                assert width == WIDTH;
                for (int j = 0; j < width; j++) {
                    int pixelInx = l * IMAGE_SIZE + i * height + j;
                    long grayscale = (long) (img.getRGB(j, i) & 0xffffffffl);
                    x[pixelInx] = grayscale;
                    if (cnt < 10) {
                        sb.append(grayscale).append("[").append(pixelInx).append("]");
                        cnt++;
                    }
                }
            }
            for (int i = 0; i < y.length; i++) {
                if (filenames[l].contains("non-bars")) {
                    y[i] = 0;
                    y[i + 1] = 1;
                } else {
                    y[i] = 1;
                    y[i + 1] = 0;
                }
            }
            // System.out.println(sb.toString());
        }

    }

    public static class NeuralNetworkConfig {
        final GMatrix theta1;
        final GMatrix theta2;
        final double[] trainingCosts;
        final double lambda;
        final Normalized nrm;

        public NeuralNetworkConfig(GMatrix theta1, GMatrix theta2, double[] costs, double lambda, Normalized nrm) {
            this.theta1 = theta1;
            this.theta2 = theta2;
            this.trainingCosts = costs;
            this.lambda = lambda;
            this.nrm = nrm;
        }
    }

    public int testset(NeuralNetworkConfig config, GMatrix x_test_mtx, GMatrix y_test_mtx) {
        TrainingStep test_s = doTrainingStep(x_test_mtx, y_test_mtx, config.theta1, config.theta2, config.lambda);
//        System.out.println("x");
//        System.out.print(x_test_mtx);
//        System.out.println("y");
//        System.out.print(y_test_mtx);
//        System.out.print(test_s.h);
        GMatrix vectorize = vectorize(test_s.h);
//        System.out.print(vectorize);
//        System.out.println("-----------------------------");
        return compare(y_test_mtx, vectorize);
    }

    private int compare(GMatrix expected, GMatrix actual) {
        if (expected.getNumCol() != actual.getNumCol())
            return 0;
        if (expected.getNumRow() != actual.getNumRow())
            return 0;
        for(int i=0;i<actual.getNumRow();i++)
            for(int j=0;j<actual.getNumCol();j++)
                if (Math.round(expected.getElement(i, j)) != Math.round(actual.getElement(i, j)))
                    return 0;
        return 1;
    }

    public NeuralNetworkConfig train(GMatrix x_mtx, GMatrix y_mtx) throws Exception {
        GMatrix theta1 = initTheta(features, x_mtx.getNumCol() + 1);
        GMatrix theta2 = initTheta(y_mtx.getNumCol(), features + 1);
        // System.out.println(theta2);

        long timing = 0;
        double alpha = 0.001;
        double previousCost = Double.MAX_VALUE;
        double totalJdrop = 0;
        boolean numericGradients = false;

        double lambda = 1.;
        double[] stepCosts = new double[NUM_TRIES];

         Normalized nrm = featureNormalize(x_mtx);
//         System.out.println("X mtx: " + size(x_mtx));
//         System.out.println(x_mtx);
         System.out.println("X nrm:"+size(nrm.values));
         x_mtx = nrm.values;
//         printMatrixByRows(nrm.values);

        for (int k = 0; k < NUM_TRIES; k++) {
            long currentTimeMillis = System.currentTimeMillis();
            TrainingStep step = doTrainingStep(x_mtx, y_mtx, theta1, theta2, lambda);
            long n = System.currentTimeMillis() - currentTimeMillis;
            timing += n;
            double diff = previousCost - step.cost;
            totalJdrop += diff;
//            if (diff < 0)
//                System.out.println("! OLOLO <!" + " " + diff);
            previousCost = step.cost;

            if (numericGradients)
                compareWithNumericGradients(step.grad1, step.grad2, x_mtx, y_mtx, theta1, theta2, 1);

            stepCosts[k] = step.cost;

            theta1.sub(scalarProduct(step.grad1, alpha));
            theta2.sub(scalarProduct(step.grad2, alpha));
        }

        return new NeuralNetworkConfig(theta1, theta2, stepCosts, lambda, nrm);
    }

    private void compareWithNumericGradients(GMatrix grad1, GMatrix grad2, GMatrix x_mtx, GMatrix y_mtx,
            GMatrix theta1, GMatrix theta2, int lambda) {
        double epsilon = 0.0001;
        GMatrix theta1_eps = new GMatrix(grad1.getNumRow(), grad1.getNumCol());
        theta1_eps.set(theta1);
        List<Double> gradPlus = new ArrayList<Double>();
        List<Double> gradMinus = new ArrayList<Double>();

        int cnt = 0;
        for (int i = 0; i < grad1.getNumRow(); i++)
            for (int j = 1; j < grad1.getNumCol(); j++) {
                double theta_i_j = theta1_eps.getElement(i, j);
                theta1_eps.setElement(i, j, theta_i_j + epsilon);
                TrainingStep plus = doTrainingStep(x_mtx, y_mtx, theta1_eps, theta2, lambda);
                theta1_eps.setElement(i, j, theta_i_j - epsilon);
                TrainingStep minus = doTrainingStep(x_mtx, y_mtx, theta1_eps, theta2, lambda);
                double numericcGrad = (plus.cost - minus.cost) / (2 * epsilon);
                double analyticGrad = grad1.getElement(i, j);
                gradMinus.add(numericcGrad - analyticGrad);
                gradPlus.add(numericcGrad + analyticGrad);
                // if (cnt++ < 2) {
                // System.out.format("%.9f %.9f\n", numericcGrad, analyticGrad);
                // }
                theta1_eps.setElement(i, j, theta_i_j);
            }

        GMatrix theta2_eps = new GMatrix(grad2.getNumRow(), grad2.getNumCol());
        theta2_eps.set(theta2);
        for (int i = 0; i < grad2.getNumRow(); i++)
            for (int j = 1; j < grad2.getNumCol(); j++) {
                double theta_i_j = theta2_eps.getElement(i, j);
                theta2_eps.setElement(i, j, theta_i_j + epsilon);
                TrainingStep plus = doTrainingStep(x_mtx, y_mtx, theta1, theta2_eps, lambda);
                theta2_eps.setElement(i, j, theta_i_j - epsilon);
                TrainingStep minus = doTrainingStep(x_mtx, y_mtx, theta1, theta2_eps, lambda);
                double numericcGrad = (plus.cost - minus.cost) / (2 * epsilon);
                double analyticGrad = grad2.getElement(i, j);
                gradMinus.add(numericcGrad - analyticGrad);
                gradPlus.add(numericcGrad + analyticGrad);
                // if (cnt++ < 4) {
                // System.out.format("%.9f %.9f\n", numericcGrad, analyticGrad);
                // }
                theta2_eps.setElement(i, j, theta_i_j);
            }

        double numerator = 0;
        double denominator = 0;
        for (int i = 0; i < gradPlus.size(); i++) {
            numerator += Math.pow(gradMinus.get(i), 2);
            denominator += Math.pow(gradPlus.get(i), 2);
        }
        double diff = Math.sqrt(numerator) / Math.sqrt(denominator);
        System.out.format("diff == %.14f\n", diff);
    }

    private GMatrix initTheta(int outputLayerSize, int inputLayerSize) {
        double init_epsilon = 1;// Math.sqrt(6.) / Math.sqrt(outputLayerSize +
                                // inputLayerSize);
        GMatrix mtx = new GMatrix(outputLayerSize, inputLayerSize);
        for (int i = 0; i < mtx.getNumRow(); i++)
            for (int j = 0; j < mtx.getNumCol(); j++) {
                double val = -init_epsilon + (Math.random() * ((2 * init_epsilon - init_epsilon) + 1.));

                mtx.setElement(i, j, val);
            }

//        System.out.println(mtx);
        return mtx;
    }

    public static class TrainingStep {
        final double cost;
        final GMatrix h;
        final GMatrix grad1;
        final GMatrix grad2;

        public TrainingStep(double j, GMatrix h, GMatrix grad1, GMatrix grad2) {
            this.cost = j;
            this.h = h;
            this.grad1 = grad1;
            this.grad2 = grad2;
        }
    }

    public static String size(GMatrix mtx) {
        return mtx.getNumRow() + "x" + mtx.getNumCol();
    }

    private TrainingStep doTrainingStep(GMatrix x, GMatrix y, GMatrix theta1, GMatrix theta2, double lambda) {
        // double lambda = 1;
        // size(Theta1) = [25 401]
        // size(Theta2) = [10 26]
        // GMatrix Theta1 = new GMatrix(features, features+1); // 10x(10+1)
        // GMatrix Theta2 = new GMatrix(classes, features+1); // 10x(10+1)
        GMatrix theta1_t = transpose(theta1);
        GMatrix theta2_t = transpose(theta2);

        // GMatrix X = new GMatrix(m, features);
        // TODO: Y = eye(num_labels)(y,:);
        // GMatrix Y = new GMatrix(m, classes); // 50000 x 1
        // X1 = [ones(m, 1), X];
        GMatrix X1 = prependColumn(x, 1); // 50000 x (2+1)
        // GMatrix Z2 = new GMatrix(m, features); // 50000 x 2
        // Z2 = X1 * Theta1';
        // Z2.mul(X1, theta1_t);
        GMatrix Z2 = product(X1, theta1_t);
        

        // A2 = [ones(m, 1), sigmoid(Z2)];
        GMatrix A2 = prependColumn(apply(Z2, mtxSigmoid), 1); // 50000x(2+1)

        // GMatrix Z3 = new GMatrix(m, classes); // 50000 x1

        // Z3 = A2 * Theta2';
        // Z3.mul(A2, theta2_t);
        GMatrix Z3 = product(A2, theta2_t);
        // A3 = sigmoid(Z3);
        GMatrix A3 = apply(Z3, mtxSigmoid);
        // System.out.println(size(A3));
        // A3 = vectorize(A3);
        // System.out.println(A3);

        // case1 = -Y .* log(A3);
        // [50000, 10] = size(case1);
        GMatrix case1 = dotProduct(apply(y, mtxNeg), apply(A3, mtxLog));

        // case2 = (1 - Y) .* log(1-A3);
        // [50000, 10] = size(case2);
        GMatrix oneMinusY = apply(y, mtxOneMinus);
        GMatrix case2 = dotProduct(oneMinusY, apply(apply(A3, mtxOneMinus), mtxLog));
        // System.out.println("\n++++++++++++++++++++");
        // System.out.print(A3);
        // System.out.println("++++++++++++++++++++");
        // System.out.print(oneMinusY);
        // System.out.println("--------------------");
        // System.out.print(case1);
        // System.out.println("--------------------");
        // System.out.print(case2);
        // System.out.println("====================\n");

        // J = sum(sum(case1 - case2, 2))/m;

        GMatrix rowSums = sumRows(substract(case1, case2));
        double sum = 0;
        for (int i = 0; i < rowSums.getNumRow(); i++)
            sum += rowSums.getElement(i, 0);
        double J = sum / m;

        // sum_Theta1 = sum(sum(Theta1(:,2:end).^2,2));
        double sumTheta1 = sumAllSquared(theta1);
        // sum_Theta2 = sum(sum(Theta2(:,2:end).^2,2));
        double sumTheta2 = sumAllSquared(theta2);

        // regularization = lambda * (sum_Theta1 + sum_Theta2)/(2*m);
        // J = J + regularization;
        J += lambda * (sumTheta1 + sumTheta2) / (2 * m);

        // delta_3 = A3 - Y; % 5000x10
        GMatrix delta3 = substract(A3, y);

        // delta_2 = delta_3 * Theta2(:,2:end) .* sigmoidGradient(Z2); % 5000x10
        // * 10x25 .* 5000x25 = 5000x25
        GMatrix delta2 = dotProduct(product(delta3, omitFirstColumn(theta2)), apply(Z2, mtxSigmoidGradient));

        // Theta1_with_0 = [zeros(hidden_layer_size,1), Theta1(:,2:end)]; %
        // replace first colum with 0s;
        GMatrix Theta1_with_0 = prependColumn(omitFirstColumn(theta1), 0);
        // Theta2_with_0 = [zeros(num_labels,1), Theta2(:,2:end)]; % replace
        // first column with 0s
        GMatrix Theta2_with_0 = prependColumn(omitFirstColumn(theta2), 0);

        // Theta1_grad = (delta_2' * X1)./m; % ((25x5000 * 5000x401))./5000 =
        // (25x401)./5000
        GMatrix Theta1_grad = product(transpose(delta2), X1);

        // Theta1_grad = ((delta_2' * X1) + lambda*Theta1_with_0)./m; %
        // ((25x5000 * 5000x401) + lambda*(25x401))./5000 = (25x401)./5000
        Theta1_grad.add(scalarProduct(Theta1_with_0, lambda));
        GMatrix grad1 = scalarProduct(Theta1_grad, (double) 1. / m);

        // Theta2_grad = (delta_3' * A2)./m; % ((10x5000 * 5000x26))./5000 =
        // (10x26)./5000
        GMatrix Theta2_grad = product(transpose(delta3), A2);

        // Theta2_grad = ((delta_3' * A2) + lambda*Theta2_with_0)./m; %
        // ((10x5000 * 5000x26) + lambda*(10x26))./5000 = (10x26)./5000
        Theta2_grad.add(scalarProduct(Theta2_with_0, lambda));
        GMatrix grad2 = scalarProduct(Theta2_grad, (double) 1. / m);

        return new TrainingStep(J, A3, grad1, grad2);
    }

    private void printMatrixByRows(GMatrix mtx) {
        for (int i = 0; i < mtx.getNumRow(); i++) {
            int n = mtx.getNumCol();
            double[] row = new double[n];
            mtx.getRow(i, row);
            for (int j = 0; j < n; j++)
                System.out.print(row[j] + " ");
            System.out.println();
        }
    }
    
    public static GMatrix featureNormalize(GMatrix x_orig, Normalized norm) {
        GMatrix x_mean = substractVector(x_orig, norm.mu);
        
        return scalarVectorDivide(x_mean, norm.sigma);
    }

    public static Normalized featureNormalize(GMatrix X) {
        // m = size(X, 1);
        //
        // mu = mean(X);
        // X_norm = X - kron(mu, ones(m, 1));
        //
        // sigma = std(X_norm) = sqrt(sumsq(x-mean(x))/(length(x)-1));
        // X_norm = X_norm ./ kron(sigma, ones(m, 1));
        GMatrix mu = mean(X);
        GMatrix X_mean = substractVector(X, mu);
        GMatrix sigma = std(X_mean);
        GMatrix X_norm = scalarVectorDivide(X_mean, sigma);
        return new Normalized(X_norm, mu, sigma);
    }

    public static class Normalized {
        GMatrix values;
        GMatrix mu;
        GMatrix sigma;

        public Normalized(GMatrix vals, GMatrix mu, GMatrix sigma) {
            this.values = vals;
            this.mu = mu;
            this.sigma = sigma;
        }
    }

    public static void main(String[] args) throws Exception {
        TrainNN nn = new TrainNN();
//        nn.controlledImageSetInit();
//        GMatrix x = new GMatrix(nn.m, nn.inputLayerSize, nn.x);
//        GMatrix y = new GMatrix(nn.m, nn.classes, nn.y);
        Data data = getMixedData();
        NeuralNetworkConfig cfg = nn.train(data.x, data.y);
        int accuracy = 0;
        
//        double[] bars = longToDoubleArray(imageToGrayscale(resize(readImage(filenames[0]), SMALL_HEIGHT, SMALL_WIDTH))); //bars
//        double[] nonBars =longToDoubleArray(imageToGrayscale(resize(readImage(filenames[1]), SMALL_HEIGHT, SMALL_WIDTH))); //non-bars
        
        
        List<File> testSet = listSetFiles(TEST_DIR);
        for (File f : testSet) {
            
            double[] y_test = new double[nn.classes];
            y_test[0]   = f.getAbsolutePath().contains("non-bar") ? 0. : 1.;
            y_test[1]   = f.getAbsolutePath().contains("non-bar") ? 1. : 0.;
            
            GMatrix x_norm = featureNormalize(readData(f), cfg.nrm);
            GMatrix y = new GMatrix(1, nn.classes, y_test);
            int val = nn.testset(cfg, x_norm, y);
            if (val == 0){
                System.err.println(f.getAbsolutePath());
                System.err.print(y);
                System.err.println("-------------------------\n");
            }
            accuracy += val;
            
        }
        
        System.out.println(((((float)accuracy)/testSet.size())*100)+"% accuracy");
        showChart(cfg.trainingCosts);
    }
}
