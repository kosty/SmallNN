import static com.AlgebraUtil.apply;
import static com.AlgebraUtil.dotProduct;
import static com.AlgebraUtil.mtxLog;
import static com.AlgebraUtil.mtxNeg;
import static com.AlgebraUtil.mtxOneMinus;
import static com.AlgebraUtil.mtxSigmoid;
import static com.AlgebraUtil.mtxSigmoidGradient;
import static com.AlgebraUtil.omitFirstColumn;
import static com.AlgebraUtil.prependColumn;
import static com.AlgebraUtil.product;
import static com.AlgebraUtil.scalarProduct;
import static com.AlgebraUtil.substract;
import static com.AlgebraUtil.sumAllSquared;
import static com.AlgebraUtil.sumRows;
import static com.AlgebraUtil.transpose;

import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import javax.imageio.ImageIO;
import javax.swing.JFrame;
import javax.vecmath.GMatrix;

import org.junit.Before;
import org.junit.Test;

public class TrainNN {
    private static final int NUM_TRIES = 2;
    Random r = new Random();

    int features = 20;
    int classes = 2;

    private static final int HEIGHT = 60;
    private static final int WIDTH = 72;

    private static final int SMALL_HEIGHT = 5;
    private static final int SMALL_WIDTH = 6;

    public static final int IMAGE_SIZE = WIDTH * HEIGHT;
    public static final int SMALL_IMAGE_SIZE = SMALL_WIDTH * SMALL_HEIGHT;

    int m = filenames.length; // 50000;
    int inputLayerSize = IMAGE_SIZE;
    double[] x = new double[m * inputLayerSize];
    double[] y = new double[m * classes];

    int m_test = testcases.length;

    public void randomInit() {
        m = 2;
        inputLayerSize = 2;
        features = 2;
        x = new double[m * inputLayerSize];
        y = new double[m * classes];
        for (int i = 0; i < x.length; i++) {
            x[i] = r.nextDouble();
        }
        for (int i = 0; i < y.length; i++) {
            if (r.nextBoolean()) {
                y[i] = 0;
                y[i + 1] = 1;
            } else {
                y[i] = 1;
                y[i + 1] = 0;
            }
        }
    }

    @Before
    public void smallImageInit() throws Exception {
        inputLayerSize = SMALL_IMAGE_SIZE;
        x = new double[m * inputLayerSize];
        y = new double[m * classes];

        for (int l = 0; l < m; l++) {
            int cnt = 0;

            BufferedImage img = resize(ImageIO.read(new File(filenames[l])), SMALL_WIDTH, SMALL_HEIGHT);
            // System.out.println(filenames[l]);
            int height = img.getHeight();
            assert height == 5;
            for (int i = 0; i < height; i++) {
                int width = img.getWidth();
                assert width == 6;
                for (int j = 0; j < width; j++) {
                    int pixelInx = l * inputLayerSize + i * height + j;
                    int rgb = img.getRGB(j, i);
                    double grayscale = (double) rgb / (double) Integer.MAX_VALUE;
                    x[pixelInx] = grayscale;
                    // System.out.format("%.4f [%d] ", grayscale, pixelInx);
                    cnt++;
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
                    int rgb = img.getRGB(j, i);
                    double grayscale = (double) rgb / (double) Integer.MAX_VALUE;
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

    private static BufferedImage resize(BufferedImage image, int width, int height) {
        BufferedImage resizedImage = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
        Graphics2D g = resizedImage.createGraphics();
        g.setRenderingHint(RenderingHints.KEY_INTERPOLATION, RenderingHints.VALUE_INTERPOLATION_BILINEAR);

        g.setRenderingHint(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);

        g.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g.drawImage(image, 0, 0, width, height, null);
        g.dispose();
        return resizedImage;
    }

    @Test
    public double[] test() throws Exception {
        GMatrix x_mtx = new GMatrix(m, inputLayerSize, x);
        GMatrix y_mtx = new GMatrix(m, classes, y);

        GMatrix theta1 = initTheta(features, inputLayerSize + 1);
        GMatrix theta2 = initTheta(classes, features + 1);
        // System.out.println(theta2);

        long timing = 0;
        double alpha = 0.01;
        double previousCost = Double.MAX_VALUE;
        double totalJdrop = 0;
        boolean numericGradients = true;

        double lambda = 1.;
        double[] stepCosts = new double[NUM_TRIES];

        for (int k = 0; k < NUM_TRIES; k++) {
            long currentTimeMillis = System.currentTimeMillis();
            SimulationStep step = doTrainingStep(x_mtx, y_mtx, theta1, theta2, lambda);
            long n = System.currentTimeMillis() - currentTimeMillis;
            timing += n;
            double diff = previousCost - step.cost;
            totalJdrop += diff;
            if (diff < 0)
                System.out.println("! OLOLO <!" + " " + diff);
            previousCost = step.cost;

            if (numericGradients)
                compareWithNumericGradients(step.grad1, step.grad2, x_mtx, y_mtx, theta1, theta2, 1);

            stepCosts[k] = step.cost;

            theta1.sub(scalarProduct(step.grad1, alpha));
            theta2.sub(scalarProduct(step.grad2, alpha));
        }

        for (int l = 0; l < testcases.length; l++) {

            BufferedImage img = resize(ImageIO.read(new File(testcases[l])), SMALL_WIDTH, SMALL_HEIGHT);
            double[] x_test = new double[inputLayerSize];

            int height = img.getHeight();
            assert height == HEIGHT;
            for (int i = 0; i < height; i++) {
                int width = img.getWidth();
                assert width == WIDTH;
                for (int j = 0; j < width; j++) {
                    int pixelInx = i * height + j;
                    int rgb = img.getRGB(j, i);
                    x_test[pixelInx] = (double) rgb / (double) Integer.MAX_VALUE;
                    // System.out.print(x_test[pixelInx]+" ");
                }
            }
            double[] y_test = testcases[l].contains("non-bars") ? new double[] { 1., 0. } : new double[] { 0., 1. };
            GMatrix x_test_mtx = new GMatrix(1, inputLayerSize, x_test);
            GMatrix y_test_mtx = new GMatrix(1, classes, y_test);
            SimulationStep test_s = doTrainingStep(x_test_mtx, y_test_mtx, theta1, theta2, lambda);
            System.out.format("%.1f : %s\n", y_test[0], testcases[l]);
            System.out.println(test_s.h);
        }
        return stepCosts;
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
                SimulationStep plus = doTrainingStep(x_mtx, y_mtx, theta1_eps, theta2, lambda);
                theta1_eps.setElement(i, j, theta_i_j - epsilon);
                SimulationStep minus = doTrainingStep(x_mtx, y_mtx, theta1_eps, theta2, lambda);
                double numericcGrad = (plus.cost - minus.cost) / (2 * epsilon);
                double analyticGrad = grad1.getElement(i, j);
                gradMinus.add(numericcGrad - analyticGrad);
                gradPlus.add(numericcGrad + analyticGrad);
                if (cnt++ < 2) {
                    System.out.format("%.9f %.9f\n", numericcGrad, analyticGrad);
                }
                theta1_eps.setElement(i, j, theta_i_j);
            }

        GMatrix theta2_eps = new GMatrix(grad2.getNumRow(), grad2.getNumCol());
        theta2_eps.set(theta2);
        for (int i = 0; i < grad2.getNumRow(); i++)
            for (int j = 1; j < grad2.getNumCol(); j++) {
                double theta_i_j = theta2_eps.getElement(i, j);
                theta2_eps.setElement(i, j, theta_i_j + epsilon);
                SimulationStep plus = doTrainingStep(x_mtx, y_mtx, theta1, theta2_eps, lambda);
                theta2_eps.setElement(i, j, theta_i_j - epsilon);
                SimulationStep minus = doTrainingStep(x_mtx, y_mtx, theta1, theta2_eps, lambda);
                double numericcGrad = (plus.cost - minus.cost) / (2 * epsilon);
                double analyticGrad = grad2.getElement(i, j);
                gradMinus.add(numericcGrad - analyticGrad);
                gradPlus.add(numericcGrad + analyticGrad);
                if (cnt++ < 4) {
                    System.out.format("%.9f %.9f\n", numericcGrad, analyticGrad);
                }
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
        double init_epsilon = Math.sqrt(6.) / Math.sqrt(outputLayerSize + inputLayerSize);
        GMatrix mtx = new GMatrix(outputLayerSize, inputLayerSize);
        for (int i = 0; i < mtx.getNumRow(); i++)
            for (int j = 0; j < mtx.getNumCol(); j++) {
                double val = r.nextGaussian() * (2 * init_epsilon) - init_epsilon;
                mtx.setElement(i, j, val);
            }

        System.out.println(mtx);
        return mtx;
    }

    public static class SimulationStep {
        final double cost;
        final GMatrix h;
        final GMatrix grad1;
        final GMatrix grad2;

        public SimulationStep(double j, GMatrix h, GMatrix grad1, GMatrix grad2) {
            this.cost = j;
            this.h = h;
            this.grad1 = grad1;
            this.grad2 = grad2;
        }
    }

    public static String size(GMatrix mtx) {
        return mtx.getNumRow() + "x" + mtx.getNumCol();
    }

    private SimulationStep doTrainingStep(GMatrix x, GMatrix y, GMatrix theta1, GMatrix theta2, double lambda) {
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

        return new SimulationStep(J, A3, grad1, grad2);
    }

    static final String[] filenames = {
            "/Users/arsen/References/misc/bars/72x60/61.scaled.72x60.png", // [1
                                                                           // 0]
            "/Users/arsen/References/misc/non-bars/72x60/bar-041.2.72x60.png",// [0
                                                                              // 1]
            "/Users/arsen/References/misc/bars/72x60/61.40x40.to.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-041.2.72x60.png",// [0
                                                                              // 1]

            "/Users/arsen/References/misc/bars/72x60/61.ffmpeg.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/61.rescaled.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/badteacher_2011_rated_hd_16x9_185_2398_english_6207_JPEG2000.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/burlesque_2010_hd_16x9_240_2398_english_6146_JPEG2000.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/fogthe_unrated_2005_hd_16x9_235_2398_english_74_JPEG2000.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/groove_2000_hd_16x9_178_2398_english_4814_JPEG2000.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/harttohart_sd_63_las_no7241_test_JPEG2000.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/jeopardy_5301_test_JPEG2000.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/replacementkillers_2005_ec_CC_16x9_240_english_2136_JPEG2000.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/socialnetworkthe_2010_hd_16x9_240_2398_english_5135_JPEG2000.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/steve_dont_get_nun_1361_tv_2727.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/swat_1975_16_sd_4x3_133_25_english_1180_JPEG2000.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-001.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-002.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-003.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-004.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-005.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-006.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-007.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-008.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-009.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-010.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-011.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-012.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-013.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-014.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-015.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-016.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-017.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-018.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-019.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-020.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-021.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-022.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-023.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-024.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-025.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-026.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-027.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-028.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-029.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-030.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-031.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-032.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-033.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-034.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-035.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-036.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-037.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-038.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-039.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-040.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-041.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-042.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-043.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-044.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-045.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-046.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-047.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-048.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-049.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-050.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-051.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-052.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-053.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-054.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-055.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-056.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-057.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-058.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-059.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-060.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-061.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-062.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-063.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-064.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-065.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-066.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-067.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-068.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-069.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-070.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-071.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-072.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-073.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-074.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-075.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-076.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-077.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-078.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-079.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-080.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-042.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-043.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-044.2.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-081.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-082.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-045.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-046.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-047.2.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-083.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-085.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-040.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-048.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-049.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-050.2.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-088.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-089.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-002.2.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-086.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-087.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-003.2.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-084.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-004.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-005.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-006.2.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-090.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-091.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-092.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-093.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-094.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-095.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-096.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-100.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-001.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-007.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-008.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-009.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-010.2.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-097.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-011.2.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-099.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-012.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-013.2.72x60.png",
            "/Users/arsen/References/misc/bars/72x60/tvproxy-bar-098.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-014.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-015.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-016.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-017.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-018.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-019.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-020.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-021.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-022.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-023.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-024.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-025.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-026.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-027.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-028.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-029.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-030.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-031.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-032.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-033.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-034.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-035.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-036.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-037.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-038.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-039.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-051.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-052.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-053.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-054.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-055.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-056.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-057.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-058.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-059.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-060.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-061.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-062.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-063.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-064.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-065.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-066.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-067.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-068.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-069.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-070.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-071.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-072.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-073.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-074.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-075.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-076.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-077.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-078.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-079.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-080.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-081.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-082.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-083.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-084.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-085.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-086.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-087.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-088.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-089.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-090.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-091.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-092.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-093.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-094.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-095.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-096.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-097.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-098.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-099.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/bar-100.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/non-bar-001.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/non-bar-001.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/non-bar-002.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/non-bar-002.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/non-bar-003.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/non-bar-003.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/non-bar-004.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/non-bar-004.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/non-bar-005.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/non-bar-005.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/non-bar-006.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/non-bar-006.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/non-bar-007.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/non-bar-007.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/non-bar-008.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/non-bar-008.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/non-bar-009.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/non-bar-009.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/non-bar-010.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/non-bar-010.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/non-bar-011.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/non-bar-011.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/non-bar-012.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/non-bar-012.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/non-bar-013.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/non-bar-013.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/non-bar-014.2.72x60.png",
            "/Users/arsen/References/misc/non-bars/72x60/non-bar-014.72x60.png"

    };

    public static final String[] testcases = { "/Users/arsen/References/misc/bars/bars-test/72x60/bar-001.72x60.png",
            "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-100.72x60.png",
    /*
     * "/Users/arsen/References/misc/bars/bars-test/72x60/bar-002.72x60.png" ,
     * "/Users/arsen/References/misc/bars/bars-test/72x60/bar-003.72x60.png" ,
     * "/Users/arsen/References/misc/bars/bars-test/72x60/bar-004.72x60.png" ,
     * "/Users/arsen/References/misc/bars/bars-test/72x60/bar-005.72x60.png" ,
     * "/Users/arsen/References/misc/bars/bars-test/72x60/bar-006.72x60.png" ,
     * "/Users/arsen/References/misc/bars/bars-test/72x60/bar-007.72x60.png" ,
     * "/Users/arsen/References/misc/bars/bars-test/72x60/bar-008.72x60.png" ,
     * "/Users/arsen/References/misc/bars/bars-test/72x60/bar-009.72x60.png" ,
     * "/Users/arsen/References/misc/bars/bars-test/72x60/bar-010.72x60.png" ,
     * "/Users/arsen/References/misc/bars/bars-test/72x60/bar-011.72x60.png" ,
     * "/Users/arsen/References/misc/bars/bars-test/72x60/bar-012.72x60.png" ,
     * "/Users/arsen/References/misc/bars/bars-test/72x60/bar-013.72x60.png" ,
     * "/Users/arsen/References/misc/bars/bars-test/72x60/bar-014.72x60.png" ,
     * "/Users/arsen/References/misc/bars/bars-test/72x60/bar-015.72x60.png" ,
     * "/Users/arsen/References/misc/bars/bars-test/72x60/bar-016.72x60.png" ,
     * "/Users/arsen/References/misc/bars/bars-test/72x60/bar-017.72x60.png" ,
     * "/Users/arsen/References/misc/bars/bars-test/72x60/bar-018.72x60.png" ,
     * "/Users/arsen/References/misc/bars/bars-test/72x60/bar-019.72x60.png" ,
     * "/Users/arsen/References/misc/bars/bars-test/72x60/bar-020.72x60.png" ,
     * "/Users/arsen/References/misc/bars/bars-test/72x60/bar-021.72x60.png" ,
     * "/Users/arsen/References/misc/bars/bars-test/72x60/bar-022.72x60.png" ,
     * "/Users/arsen/References/misc/bars/bars-test/72x60/bar-023.72x60.png" ,
     * "/Users/arsen/References/misc/bars/bars-test/72x60/bar-024.72x60.png" ,
     * "/Users/arsen/References/misc/bars/bars-test/72x60/bar-025.72x60.png" ,
     * "/Users/arsen/References/misc/bars/bars-test/72x60/bar-026.72x60.png" ,
     * "/Users/arsen/References/misc/bars/bars-test/72x60/bar-027.72x60.png" ,
     * "/Users/arsen/References/misc/bars/bars-test/72x60/bar-028.72x60.png" ,
     * "/Users/arsen/References/misc/bars/bars-test/72x60/bar-029.72x60.png" ,
     * "/Users/arsen/References/misc/bars/bars-test/72x60/bar-030.72x60.png" ,
     * "/Users/arsen/References/misc/bars/bars-test/72x60/bar-031.72x60.png" ,
     * "/Users/arsen/References/misc/bars/bars-test/72x60/bar-032.72x60.png" ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-001.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-002.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-003.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-004.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-005.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-006.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-007.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-008.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-009.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-010.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-011.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-012.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-013.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-014.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-015.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-016.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-017.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-018.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-019.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-020.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-021.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-022.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-023.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-024.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-025.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-026.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-027.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-028.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-029.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-030.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-031.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-032.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-033.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-034.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-035.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-036.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-037.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-038.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-039.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-040.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-041.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-042.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-043.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-044.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-045.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-046.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-047.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-048.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-049.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-050.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-051.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-052.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-053.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-054.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-055.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-056.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-057.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-058.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-059.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-060.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-061.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-062.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-063.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-064.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-065.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-066.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-067.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-068.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-069.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-070.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-071.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-072.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-073.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-074.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-075.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-076.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-077.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-078.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-079.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-080.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-081.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-082.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-083.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-084.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-085.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-086.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-087.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-088.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-089.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-090.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-091.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-092.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-093.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-094.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-095.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-096.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-097.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-098.72x60.png"
     * ,
     * "/Users/arsen/References/misc/non-bars/non-bars-test/72x60/bar-099.72x60.png"
     * ,
     */
    };

    public static void main(String[] args) throws Exception {
        TrainNN nn = new TrainNN();
        nn.smallImageInit();
        double[] data = nn.test();
        JFrame f = new JFrame();
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        f.getContentPane().add(new PlotTest(data));
        f.setSize(400, 400);
        f.setLocation(200, 200);
        f.setVisible(true);
    }
}
