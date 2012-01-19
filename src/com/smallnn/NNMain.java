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

import static com.smallnn.AlgebraUtil.vectorize;
import static com.smallnn.input.TrainDataUtil.getData;
import static com.smallnn.input.TrainDataUtil.getTestData;

import javax.vecmath.GMatrix;

import com.smallnn.input.ImageUtil.Resolution;
import com.smallnn.input.TrainDataUtil.Data;
import com.smallnn.input.TrainDataUtil.TestData;
import com.smallnn.output.ChartUtil;

public class NNMain {

    private static final int BOOST_LEVEL = 1;
    private static final int MANY_REPETITIONS = 50000;

    public static void main(String[] args) throws Exception {
        String kit = args[0];
        Resolution res = new Resolution(71, 40);
        
        System.out.println("Preparing test data...");
        Data data = getData(kit, res, res, BOOST_LEVEL);
        TestData testData = getTestData(kit, res, res);

        NeuralNetwork nn = new SingleLayerNetwork(30, data.y.getNumCol(), data.x.getNumCol());
        
        
        Double[] trainingCosts = measureTrainingTime(data, nn);

        measurePrecisionAndRecall(testData, nn);
        ChartUtil.showChart(trainingCosts);
        
        measureBulkActivationTime(testData, nn);
    }

    private static Double[] measureTrainingTime(Data data, NeuralNetwork nn) throws Exception {
        
        /* Learning rate for gradient descent */
        double alpha = 0.6;
        long startTime = System.currentTimeMillis();
        

        /* Regularization cost: 1 */
        Double[] trainingCosts = nn.train(data.x, data.y, 1., alpha);
        double traintime = (System.currentTimeMillis() - startTime) / 1000.;
        System.out.println("Training time: " + traintime + " sec; " + (traintime / trainingCosts.length) + " sec per train, total train iterations: "+trainingCosts.length);
        return trainingCosts;
    }

    private static void measurePrecisionAndRecall(TestData testData, NeuralNetwork nn) {
        int truePositives = 0, falsePositives = 0, falseNegatives = 0, trueNegatives = 0;
        
        long astartTime = System.currentTimeMillis();
        for(int j=0;j<testData.x.length;j++){
            GMatrix vectorize = vectorize(nn.activate(testData.x[j]));

            if (compare(testData.y[j], vectorize)) {
                if (Math.abs(testData.y[j].getElement(0, 0)) == 1){
                    truePositives++;
                } else {
                    trueNegatives++;
                }
            } else {
                if (Math.abs(testData.y[j].getElement(0, 0)) == 1){
                    falseNegatives++;
                } else {
                    falsePositives++;
                }
                 System.err.println(testData.files.get(j).getAbsolutePath());
                 System.err.println("expected: ");
                 System.err.print(testData.y[j]);
                 System.err.println("actual: ");
                 System.err.println(vectorize);
                 System.err.println("-------------------------\n");
            }
        }
        double testTime = (System.currentTimeMillis() - astartTime) / 1000.;
        
        float p = truePositives/((float)(truePositives+falsePositives));// precision
        float r = truePositives/((float)(truePositives+falseNegatives));// recall
        float f1Score = 2*p*r/(p+r);
        System.out.println("Test time: " + testTime + " sec");
        System.out.println(((truePositives+trueNegatives) * 100. / testData.files.size()) + "% accuracy");
        System.out.println("Precision: "+p);
        System.out.println("Recall:    "+r);
        System.out.println("Harmonic:  "+f1Score);
    }
    
    private static void measureBulkActivationTime(TestData testData, NeuralNetwork nn){
        System.out.println("Reading movie matrix...");
        int n = testData.x[0].getNumCol();
        double[] allMovieX = new double[MANY_REPETITIONS*n];
        double[] tmp = new double[n];
        for (int l = 0; l < MANY_REPETITIONS; l++) {
            testData.x[l % testData.x.length].getRow(0, tmp);
            System.arraycopy(tmp, 0, allMovieX, l*n, n);
        }
        GMatrix x = new GMatrix(MANY_REPETITIONS, n, allMovieX);
        System.out.println("Movie matrix read.");
        
        int nums = 10;
        double durations = 0.;
        for (int k=0;k<nums;k++){
            long start = System.currentTimeMillis();
            nn.activate(x);
            long end = System.currentTimeMillis();
            durations += (double) (end - start) / 1000.;
        }
        System.out.println("Fps: " + (MANY_REPETITIONS *nums/ durations) + "; test duration: " + durations/nums + " sec");
    }

    private static boolean compare(GMatrix expected, GMatrix actual) {
        if (expected.getNumCol() != actual.getNumCol())
            return false;
        if (expected.getNumRow() != actual.getNumRow())
            return false;
        for (int i = 0; i < actual.getNumRow(); i++)
            for (int j = 0; j < actual.getNumCol(); j++)
                if (Math.round(expected.getElement(i, j)) != Math.round(actual.getElement(i, j)))
                    return false;
        return true;
    }

    public static String size(GMatrix mtx) {
        return mtx.getNumRow() + "x" + mtx.getNumCol();
    }
}
