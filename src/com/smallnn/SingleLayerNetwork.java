/**
 *  This file is part of SmallNN, a small neural network implementation
 *  Copyright (C) 2011, 2012 Arsen Kostenko <arsen.kostenko@gmail.com>
 *     
 *  SmallNN is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  SmallNN is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with SmallNN.  If not, see <http://www.gnu.org/licenses/>.
 */
package com.smallnn;

import static com.smallnn.AlgebraUtil.apply;
import static com.smallnn.AlgebraUtil.computeNormalizationParams;
import static com.smallnn.AlgebraUtil.dotProduct;
import static com.smallnn.AlgebraUtil.mtxLog;
import static com.smallnn.AlgebraUtil.mtxNeg;
import static com.smallnn.AlgebraUtil.mtxOneMinus;
import static com.smallnn.AlgebraUtil.mtxSigmoid;
import static com.smallnn.AlgebraUtil.mtxSigmoidGradient;
import static com.smallnn.AlgebraUtil.normalize;
import static com.smallnn.AlgebraUtil.omitFirstColumn;
import static com.smallnn.AlgebraUtil.prependColumn;
import static com.smallnn.AlgebraUtil.product;
import static com.smallnn.AlgebraUtil.scalarProduct;
import static com.smallnn.AlgebraUtil.substract;
import static com.smallnn.AlgebraUtil.sumAllSquared;
import static com.smallnn.AlgebraUtil.sumRows;
import static com.smallnn.AlgebraUtil.transpose;
import static java.util.Arrays.binarySearch;
import static java.util.Arrays.sort;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import javax.vecmath.GMatrix;

/**
 * Inline code comments give better insight on what's happening. In order to make them short and
 * insightful at once Matlab/Octave syntax is used.
 */
public class SingleLayerNetwork implements NeuralNetwork {
    
    private static final int PRECISION_QUEUE_SIZE = 5;

    private static final int MAX_TRAIN_ITERATIONS = 5000;

    private static final double PRECISSION = 1e-9;

    Random r = new Random();

    int inputSize = 72*40;
    int classes = 2;
    int features = 40;

    final GMatrix theta1;
    final GMatrix theta2;
    
    GMatrix X1;
    GMatrix Z2;
    GMatrix A2;
    GMatrix Z3;
    GMatrix A3;
    
    GMatrix mu;
    GMatrix sigma;

    boolean numericGradients = false;
    
    public SingleLayerNetwork(int features, int classes, int inputSize){
        this.features = features;
        this.classes  = classes;
        this.inputSize = inputSize;
        theta1 = initTheta(this.features, this.inputSize + 1);
        theta2 = initTheta(this.classes, this.features + 1);
    }
    
    public SingleLayerNetwork(GMatrix mu, GMatrix sigma, GMatrix theta1, GMatrix theta2){
        assert theta1.getNumCol() > 1; // There is at least one input, apart from bias unit
        assert theta2.getNumCol() > 1; // There is at least one hidden unit, apart from bias unit;
        assert theta1.getNumRow() +1 == theta2.getNumCol(); // number of hidden units match in both theta matrices
        assert mu.getNumCol() == theta1.getNumCol()-1;
        assert sigma.getNumCol() == theta1.getNumCol()-1;
        this.theta1 = theta1;
        this.theta2 = theta2;
        this.mu = mu;
        this.sigma = sigma;
        this.features = this.theta1.getNumRow();
        this.classes  = this.theta2.getNumRow();
        this.inputSize = this.theta1.getNumCol()-1;
    }
    
    @Override
    public Double[] train(GMatrix x, GMatrix y, double lambda, double alpha) throws Exception {
        assert this.inputSize == x.getNumCol();
        assert this.classes == y.getNumCol();
        LearningRates learningRates = new LearningRates(alpha);
        GMatrix[] muSigmaX = computeNormalizationParams(x);
        this.mu = muSigmaX[0];
        this.sigma = muSigmaX[1];
        x = muSigmaX[2];
        List<Double> stepCosts = new ArrayList<Double>();
        PushOutQueue recentCosts = new PushOutQueue(PRECISION_QUEUE_SIZE, Double.MAX_VALUE);
        
        for (int k=0; k < MAX_TRAIN_ITERATIONS; k++){
            recentCosts.add(computeCost(x, y, lambda));
            if (isPreciseEnough(recentCosts)){
                break;
            }
            
            if (isPecissionGoingUp(recentCosts)){
                System.out.println("Lowering learning rate with neighboring costs: "+ recentCosts.peek(PRECISION_QUEUE_SIZE-2)+", "+ recentCosts.peek(PRECISION_QUEUE_SIZE-1));
                alpha = learningRates.getLowerRate();
            }
            

            GMatrix[] gradients = computeGradients(x, y, lambda);
            if (numericGradients)
                compareWithNumericGradients(gradients[0], gradients[1], x, y, mu, sigma, theta1, theta2, 1);

            stepCosts.add(recentCosts.peek(PRECISION_QUEUE_SIZE-1));

            theta1.sub(scalarProduct(gradients[0], alpha));
            theta2.sub(scalarProduct(gradients[1], alpha));
        }
        Double[] result = new Double[stepCosts.size()];
        stepCosts.toArray(result);
        return result;
    }

    private boolean isPecissionGoingUp(PushOutQueue recentCosts) {
        return recentCosts.peek(PRECISION_QUEUE_SIZE-2) - recentCosts.peek(PRECISION_QUEUE_SIZE-1) < 0;
    }

    private boolean isPreciseEnough(PushOutQueue recentCosts) {
        for(int i=0;i<PRECISION_QUEUE_SIZE-1;i++)
            if (recentCosts.peek(i)-recentCosts.peek(i+1) > PRECISSION)
                return false;
        return true;
    }

    private static void compareWithNumericGradients(GMatrix grad1, GMatrix grad2, GMatrix x, GMatrix y,
            GMatrix mu, GMatrix sigma, GMatrix theta1, GMatrix theta2, int lambda) {
        double epsilon = 0.0001;
        GMatrix theta1_eps = new GMatrix(grad1.getNumRow(), grad1.getNumCol());
        theta1_eps.set(theta1);
        List<Double> gradPlus = new ArrayList<Double>();
        List<Double> gradMinus = new ArrayList<Double>();

        for (int i = 0; i < grad1.getNumRow(); i++)
            for (int j = 1; j < grad1.getNumCol(); j++) {
                double theta_i_j = theta1_eps.getElement(i, j);
                theta1_eps.setElement(i, j, theta_i_j + epsilon);
                double plusCost = new SingleLayerNetwork(mu, sigma, theta1_eps, theta2).computeCost(x, y, lambda);
                theta1_eps.setElement(i, j, theta_i_j - epsilon);
                double minusCost = new SingleLayerNetwork(mu, sigma, theta1_eps, theta2).computeCost(x, y, lambda);
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
                double plusCost = new SingleLayerNetwork(mu, sigma, theta1, theta2_eps).computeCost(x, y, lambda);
                
                theta2_eps.setElement(i, j, theta_i_j - epsilon);
                double minusCost = new SingleLayerNetwork(mu, sigma, theta1, theta2_eps).computeCost(x, y, lambda);
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

    public static GMatrix initTheta(int outputLayerSize, int inputLayerSize) {
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
        doActivate(x);
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


    @Override
    public GMatrix activate(GMatrix x) {
        return doActivate(normalize(this.mu, this.sigma, x));
    }
    
    public GMatrix doActivate(GMatrix x) {
        GMatrix theta1_t = transpose(theta1);
        GMatrix theta2_t = transpose(theta2);

        /* Prepend bias unit
         * X1 = [ones(m, 1), X]; 
         */
        X1 = prependColumn(x, 1);

        /* Compute the hidden layer - non logistic values
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
        
        return A3;
    }
    
    @Override
    public GMatrix[] getConfig(){
        return new GMatrix[]{mu, sigma, theta1, theta2};
    }
    
    public static class LearningRates {
        //double[] rates = {1.3, 0.9, 0.6, 0.3, 0.09, 0.06, 0.03, 0.009, 0.006, 0.003, 0.0009, 0.0006, 0.0003};
        double[] rates = {1.3, 0.9, 0.6, 0.3, 0.09, 0.06, 0.03};
        int idx=0;
        
        public LearningRates(double r){
            sort(rates);
            idx = binarySearch(rates, r);
        }
        
        public double getLowerRate(){
            if (idx > 0)
                idx--;
            return rates[idx];
        }
    }

}
