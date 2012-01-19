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
import static com.smallnn.AlgebraUtil.normalize;
import static com.smallnn.AlgebraUtil.omitFirstColumn;
import static com.smallnn.AlgebraUtil.prependColumn;
import static com.smallnn.AlgebraUtil.product;
import static com.smallnn.AlgebraUtil.scalarProduct;
import static com.smallnn.AlgebraUtil.substract;
import static com.smallnn.AlgebraUtil.sumAllSquared;
import static com.smallnn.AlgebraUtil.sumRows;
import static com.smallnn.AlgebraUtil.transpose;
import static com.smallnn.SingleLayerNetwork.initTheta;

import java.util.ArrayList;
import java.util.List;

import javax.vecmath.GMatrix;

public class NoHiddenLayerNetwork implements NeuralNetwork {
    
    private static final double COST_PRECISSION = 1e-30;

    private GMatrix theta1;
    private GMatrix X1;
    private GMatrix Z2;
    private GMatrix A2;
    private int classes;
    private int inputSize;
    
    GMatrix mu;
    GMatrix sigma;
    
    /* Learning rate for gradient descent */
    double alpha = 0.001;
    boolean numericGradients = false;
    
    
    
    public NoHiddenLayerNetwork(int classes, int inputSize){
        this.classes  = classes;
        this.inputSize = inputSize;
        theta1 = initTheta(this.classes, this.inputSize + 1);
    }
    
    public NoHiddenLayerNetwork(GMatrix mu, GMatrix sigma, GMatrix theta1){
        assert theta1.getNumCol() > 1; // There is at least one input, apart from bias unit
        assert mu.getNumCol() == theta1.getNumCol()-1;
        assert sigma.getNumCol() == theta1.getNumCol()-1;
        this.theta1 = theta1;
        this.mu = mu;
        this.sigma = sigma;
        this.classes  = this.theta1.getNumRow();
        this.inputSize = this.theta1.getNumCol()-1;
    }

    @Override
    public Double[] train(GMatrix x, GMatrix y, double lambda, double alpha) throws Exception {
        assert this.inputSize == x.getNumCol();
        assert this.classes == y.getNumCol();
        GMatrix[] muSigmaX = computeNormalizationParams(x);
        this.mu = muSigmaX[0];
        this.sigma = muSigmaX[1];
        x = muSigmaX[2];
        List<Double> stepCosts = new ArrayList<Double>();

        double previousCost = Double.MAX_VALUE;
        double diff = Double.MAX_VALUE;
        while (diff > COST_PRECISSION){
            double cost = computeCost(x, y, lambda);
            GMatrix[] gradients = computeGradients(x, y, lambda);
            diff = previousCost - cost;
            if (diff < COST_PRECISSION){
                continue;
            }
            previousCost = cost;

            if (numericGradients)
                compareWithNumericGradients(gradients[0], x, y, mu, sigma, theta1, 1);

            stepCosts.add(cost);

            theta1.sub(scalarProduct(gradients[0], alpha));
        }
        Double[] result = new Double[stepCosts.size()];
        stepCosts.toArray(result);
        return result;
    }

    @Override
    public GMatrix activate(GMatrix x) {
        return doActivate(normalize(this.mu, this.sigma, x));
    }
    
    private GMatrix doActivate(GMatrix x){
        GMatrix theta1_t = transpose(theta1);

        /* Prepend bias unit
         * X1 = [ones(m, 1), X]; 
         */
        X1 = prependColumn(x, 1);

        /* Compute the output layer - non logistic values
         * Z2 = X1 * Theta1'; 
         */
        Z2 = product(X1, theta1_t);
        
        /* Compute the output layer - logistic values
         * A3 = sigmoid(Z2);
         */
        A2 = apply(Z2, mtxSigmoid);
        
        return A2;
    }
    
    @Override
    public GMatrix[] getConfig(){
        return new GMatrix[]{mu, sigma, theta1};
    }
    
    private double computeCost(GMatrix x, GMatrix y, double lambda){
        doActivate(x);
        int m = x.getNumRow();
        /* Start computing the value function J, case1 corresponds values marked as 1 in training data
         * case1 = -Y .* log(A3);
         */
        GMatrix case1 = dotProduct(apply(y, mtxNeg), apply(A2, mtxLog));

        /* case2 corresponds to values marked as 0 in training data
         * case2 = (1 - Y) .* log(1-A3); 
         */
        GMatrix oneMinusY = apply(y, mtxOneMinus);
        GMatrix case2 = dotProduct(oneMinusY, apply(apply(A2, mtxOneMinus), mtxLog));

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

        // J = J + lambda * (sum_Theta1)/(2*m);
        /* Add regularization to cost function J */
        J += lambda * (sumTheta1) / (2 * m);
        
        return J;
    }
    
    private GMatrix[] computeGradients(GMatrix x, GMatrix y, double lambda) {
        int m = x.getNumRow();

        /* Start computing the theta gradients. Gradients for output layer
         * delta_3 = A3 - Y;
         */ 
        GMatrix delta3 = substract(A2, y);

        /* Bias units are not regularized, thus gradient computation must take 
         * this into account.
         * One simple way to acheave this is add 0's where regularization values for
         * bias unit should be
         * Theta1_with_0 = [zeros(hidden_layer_size,1), Theta1(:,2:end)]; 
         */
        GMatrix Theta1_with_0 = prependColumn(omitFirstColumn(theta1), 0);

        /* Theta2_grad = (delta_3' * A2)./m; */ 
        GMatrix Theta1_grad = product(transpose(delta3), X1);

        /* Theta2_grad = ((delta_3' * A2) + lambda*Theta2_with_0)./m; */ 
        Theta1_grad.add(scalarProduct(Theta1_with_0, lambda));
        GMatrix grad1 = scalarProduct(Theta1_grad, (double) 1. / m);

        return new GMatrix[]{grad1};
    }
    
    private void compareWithNumericGradients(GMatrix grad1, GMatrix x_mtx, GMatrix y_mtx, GMatrix mu, GMatrix sigma, GMatrix theta1, int lambda) {
        double epsilon = 0.0001;
        GMatrix theta1_eps = new GMatrix(grad1.getNumRow(), grad1.getNumCol());
        theta1_eps.set(theta1);
        List<Double> gradPlus = new ArrayList<Double>();
        List<Double> gradMinus = new ArrayList<Double>();

        for (int i = 0; i < grad1.getNumRow(); i++)
            for (int j = 1; j < grad1.getNumCol(); j++) {
                double theta_i_j = theta1_eps.getElement(i, j);
                theta1_eps.setElement(i, j, theta_i_j + epsilon);
                double plusCost = new NoHiddenLayerNetwork(mu, sigma, theta1_eps).computeCost(x_mtx, y_mtx, lambda);
                theta1_eps.setElement(i, j, theta_i_j - epsilon);
                double minusCost = new NoHiddenLayerNetwork(mu, sigma, theta1_eps).computeCost(x_mtx, y_mtx, lambda);
                double numericcGrad = (plusCost - minusCost) / (2 * epsilon);
                double analyticGrad = grad1.getElement(i, j);
                
                gradMinus.add(numericcGrad - analyticGrad);
                gradPlus.add(numericcGrad + analyticGrad);
                theta1_eps.setElement(i, j, theta_i_j);
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

}
