package com.smallnn;

import static com.smallnn.AlgebraUtil.apply;
import static com.smallnn.AlgebraUtil.dotProduct;
import static com.smallnn.AlgebraUtil.mtxLog;
import static com.smallnn.AlgebraUtil.mtxNeg;
import static com.smallnn.AlgebraUtil.mtxOneMinus;
import static com.smallnn.AlgebraUtil.mtxSigmoid;
import static com.smallnn.AlgebraUtil.mtxSigmoidGradient;
import static com.smallnn.AlgebraUtil.omitFirstColumn;
import static com.smallnn.AlgebraUtil.prependColumn;
import static com.smallnn.AlgebraUtil.product;
import static com.smallnn.AlgebraUtil.scalarProduct;
import static com.smallnn.AlgebraUtil.substract;
import static com.smallnn.AlgebraUtil.sumAllSquared;
import static com.smallnn.AlgebraUtil.sumRows;
import static com.smallnn.AlgebraUtil.transpose;
import static com.smallnn.input.ImageUtil.IMAGE_SIZE;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import javax.vecmath.GMatrix;

/**
 * Inline code comments give better insight on what's happening. In order to make them short and
 * insightful at once Matlab/Octave syntax is used.
 */
public class SingleLayerNetwork {
//    public static final int NUM_TRIES = 2;
    Random r = new Random();

    int features = 40;
    int classes = 2;

//    public static final int IMAGE_SIZE = WIDTH * HEIGHT;
//    public static final int SMALL_IMAGE_SIZE = SMALL_WIDTH * SMALL_HEIGHT;
    
    final GMatrix theta1;
    final GMatrix theta2;
    
    GMatrix X1;
    GMatrix Z2;
    GMatrix A2;
    GMatrix Z3;
    GMatrix A3;

//    int m = 0; 
    int inputSize = IMAGE_SIZE;
    
    /* Learning rate for gradient descent */
    double alpha = 0.001;
    boolean numericGradients = false;
    
    public SingleLayerNetwork(){
        theta1 = initTheta(this.features, this.inputSize + 1);
        theta2 = initTheta(this.classes, this.features + 1);
    }
    
    public SingleLayerNetwork(int features, int classes, int inputSize){
        this.features = features;
        this.classes  = classes;
        this.inputSize = inputSize;
        theta1 = initTheta(this.features, this.inputSize + 1);
        theta2 = initTheta(this.classes, this.features + 1);
    }
    
    public SingleLayerNetwork(GMatrix theta1, GMatrix theta2){
        assert theta1.getNumCol() > 1; // There is at least one input, apart from bias unit
        assert theta2.getNumCol() > 1; // There is at least one hidden unit, apart from bias unit;
        assert theta1.getNumRow() +1 == theta2.getNumCol(); // number of hidden units match in both theta matrices
        this.theta1 = theta1;
        this.theta2 = theta2;
        this.features = this.theta1.getNumRow();
        this.classes  = this.theta2.getNumRow();
        this.inputSize = this.theta1.getNumCol()-1;
    }

    public double[] train(GMatrix x, GMatrix y, double lambda, int numberOfSteps) throws Exception {
        assert this.inputSize == x.getNumCol();
        assert this.classes == y.getNumCol();
        double[] stepCosts = new double[numberOfSteps];

//        double previousCost = Double.MAX_VALUE;
        for (int k = 0; k < numberOfSteps; k++) {
            activate(x);
            double J = computeCost(x, y, lambda);
            GMatrix[] gradients = computeGradients(x, y, lambda);
//            double diff = previousCost - step.cost;
//            if (diff < 0)
//                System.out.println("! OLOLO <!" + " " + diff);
//            previousCost = step.cost;

            if (numericGradients)
                compareWithNumericGradients(gradients[0], gradients[1], x, y, theta1, theta2, 1);

            stepCosts[k] = J;

            theta1.sub(scalarProduct(gradients[0], alpha));
            theta2.sub(scalarProduct(gradients[1], alpha));
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

        for (int i = 0; i < grad1.getNumRow(); i++)
            for (int j = 1; j < grad1.getNumCol(); j++) {
                double theta_i_j = theta1_eps.getElement(i, j);
                theta1_eps.setElement(i, j, theta_i_j + epsilon);
                double plusCost = new SingleLayerNetwork(theta1_eps, theta2).computeCost(x_mtx, y_mtx, lambda);
                theta1_eps.setElement(i, j, theta_i_j - epsilon);
                double minusCost = new SingleLayerNetwork(theta1_eps, theta2).computeCost(x_mtx, y_mtx, lambda);
                double numericcGrad = (plusCost - minusCost) / (2 * epsilon);
                double analyticGrad = grad1.getElement(i, j);
                
                gradMinus.add(numericcGrad - analyticGrad);
                gradPlus.add(numericcGrad + analyticGrad);
                theta1_eps.setElement(i, j, theta_i_j);
            }

        GMatrix theta2_eps = new GMatrix(grad2.getNumRow(), grad2.getNumCol());
        theta2_eps.set(theta2);
        for (int i = 0; i < grad2.getNumRow(); i++)
            for (int j = 1; j < grad2.getNumCol(); j++) {
                double theta_i_j = theta2_eps.getElement(i, j);
                
                theta2_eps.setElement(i, j, theta_i_j + epsilon);
                double plusCost = new SingleLayerNetwork(theta1, theta2_eps).computeCost(x_mtx, y_mtx, lambda);
                
                theta2_eps.setElement(i, j, theta_i_j - epsilon);
                double minusCost = new SingleLayerNetwork(theta1, theta2_eps).computeCost(x_mtx, y_mtx, lambda);
                double numericcGrad = (plusCost - minusCost) / (2 * epsilon);
                double analyticGrad = grad2.getElement(i, j);
                
                gradMinus.add(numericcGrad - analyticGrad);
                gradPlus.add(numericcGrad + analyticGrad);
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

    private static GMatrix initTheta(int outputLayerSize, int inputLayerSize) {
        double init_epsilon = 1;
        /* Math.sqrt(6.) / Math.sqrt(outputLayerSize + inputLayerSize); */
        GMatrix mtx = new GMatrix(outputLayerSize, inputLayerSize);
        for (int i = 0; i < mtx.getNumRow(); i++)
            for (int j = 0; j < mtx.getNumCol(); j++) {
                double val = -init_epsilon + (Math.random() * ((2 * init_epsilon - init_epsilon) + 1.));

                mtx.setElement(i, j, val);
            }

        return mtx;
    }
    
    private double computeCost(GMatrix x, GMatrix y, double lambda){
        int m = x.getNumRow();
        /* Start computing the value function J, case1 corresponds values marked as 1 in training data
         * case1 = -Y .* log(A3);
         */
        GMatrix case1 = dotProduct(apply(y, mtxNeg), apply(A3, mtxLog));

        /* case2 corresponds to values marked as 0 in training data
         * case2 = (1 - Y) .* log(1-A3); 
         */
        GMatrix oneMinusY = apply(y, mtxOneMinus);
        GMatrix case2 = dotProduct(oneMinusY, apply(apply(A3, mtxOneMinus), mtxLog));

        /* J = sum(sum(case1 - case2, 2))/m */
        GMatrix rowSums = sumRows(substract(case1, case2));
        double sum = 0;
        for (int i = 0; i < rowSums.getNumRow(); i++)
            sum += rowSums.getElement(i, 0);
        
        double J = sum / m;

        /* Start computing the regularization values, in order to avoid overfitting,
         * bias unit is not regularized, hence Theta1(:,2:end)
         * TODO: sumAllSquared(theta1) - why all the thata1 values are summed?
         */
        // sum_Theta1 = sum(sum(Theta1(:,2:end).^2,2));
        double sumTheta1 = sumAllSquared(theta1);
        // sum_Theta2 = sum(sum(Theta2(:,2:end).^2,2));
        double sumTheta2 = sumAllSquared(theta2);

        // J = J + lambda * (sum_Theta1 + sum_Theta2)/(2*m);
        /* Add regularization to cost function J */
        J += lambda * (sumTheta1 + sumTheta2) / (2 * m);
        
        return J;
    }

    private GMatrix[] computeGradients(GMatrix x, GMatrix y, double lambda) {
        int m = x.getNumRow();

        /* Start computing the theta gradients. Gradients for output layer
         * delta_3 = A3 - Y;
         */ 
        GMatrix delta3 = substract(A3, y);

        /* Gradients for hiddent layer:
         * delta_2 = delta_3 * Theta2(:,2:end) .* sigmoidGradient(Z2); 
         */
        GMatrix delta2 = dotProduct(product(delta3, omitFirstColumn(theta2)), apply(Z2, mtxSigmoidGradient));

        /* Bias units are not regularized, thus gradient computation must take 
         * this into account.
         * One simple way to acheave this is add 0's where regularization values for
         * bias unit should be
         * Theta1_with_0 = [zeros(hidden_layer_size,1), Theta1(:,2:end)]; 
         */
        GMatrix Theta1_with_0 = prependColumn(omitFirstColumn(theta1), 0);
        /* Theta2_with_0 = [zeros(num_labels,1), Theta2(:,2:end)]; 
         * replace first column with 0s
         */
        GMatrix Theta2_with_0 = prependColumn(omitFirstColumn(theta2), 0);

        /* Actually compute gradients:
         *  Theta1_grad = (delta_2' * X1)./m; 
         */
        GMatrix Theta1_grad = product(transpose(delta2), X1);

        /* Theta1_grad = ((delta_2' * X1) + lambda*Theta1_with_0)./m; */ 
        Theta1_grad.add(scalarProduct(Theta1_with_0, lambda));
        GMatrix grad1 = scalarProduct(Theta1_grad, (double) 1. / m);

        /* Theta2_grad = (delta_3' * A2)./m; */ 
        GMatrix Theta2_grad = product(transpose(delta3), A2);

        /* Theta2_grad = ((delta_3' * A2) + lambda*Theta2_with_0)./m; */ 
        Theta2_grad.add(scalarProduct(Theta2_with_0, lambda));
        GMatrix grad2 = scalarProduct(Theta2_grad, (double) 1. / m);

        return new GMatrix[]{grad1, grad2};
    }

    public void activate(GMatrix x) {
        GMatrix theta1_t = transpose(theta1);
        GMatrix theta2_t = transpose(theta2);

        /* Prepend bias unit
         * X1 = [ones(m, 1), X]; 
         */
        X1 = prependColumn(x, 1);

        /* Compute the hiddent layer - non logistic values
         * Z2 = X1 * Theta1'; 
         */
        Z2 = product(X1, theta1_t);
        
        /* Compute the logistic values of hidden layer and prepend the bias unit
         * A2 = [ones(m, 1), sigmoid(Z2)]; 
         */
        A2 = prependColumn(apply(Z2, mtxSigmoid), 1);

        /* Compute the output layer - non logistic values
         * Z3 = A2 * Theta2';
         */
        Z3 = product(A2, theta2_t);

        /* Compute the output layer - logistic values
         * A3 = sigmoid(Z3);
         */
        A3 = apply(Z3, mtxSigmoid);
    }
}
