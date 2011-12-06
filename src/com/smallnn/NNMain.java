package com.smallnn;

import static com.smallnn.AlgebraUtil.featureNormalize;
import static com.smallnn.AlgebraUtil.vectorize;
import static com.smallnn.input.TrainDataUtil.TEST_DIR;
import static com.smallnn.input.TrainDataUtil.getMixedData;
import static com.smallnn.input.TrainDataUtil.listSetFiles;
import static com.smallnn.input.TrainDataUtil.readData;
import static com.smallnn.output.ChartUtil.showChart;

import java.io.File;
import java.util.List;

import javax.vecmath.GMatrix;

import com.smallnn.AlgebraUtil.Normalized;
import com.smallnn.input.TrainDataUtil.Data;

public class NNMain {

    public static void main(String[] args) throws Exception {
        Data data = getMixedData();

        SingleLayerNetwork nn = new SingleLayerNetwork(30, data.y.getNumCol(), data.x.getNumCol());
        long time = System.currentTimeMillis();
        Normalized nrm = featureNormalize(data.x);
        
        int numberOfSteps = 2;
        /* Regularization cost: 1 */
        double[] trainingCosts = nn.train(nrm.values, data.y, 1., numberOfSteps);
        long traintime = System.currentTimeMillis() - time;
        System.out.println("Training time: " +  (traintime / 1000) + " sec; " + (traintime / numberOfSteps) + " sec per train");

        int accuracy = 0;
        time = System.currentTimeMillis();
        List<File> testSet = listSetFiles(TEST_DIR);
        for (File f : testSet) {
            
            double[] y_test = new double[nn.classes];
            y_test[0]   = f.getAbsolutePath().contains("non-bar") ? 0. : 1.;
            y_test[1]   = f.getAbsolutePath().contains("non-bar") ? 1. : 0.;
            GMatrix y = new GMatrix(1, nn.classes, y_test);
            
            GMatrix x_norm = featureNormalize(readData(f), nrm);
            nn.activate(x_norm);
            GMatrix vectorize = vectorize(nn.A3);

            if (compare(y, vectorize)){
                accuracy++;
            } else {
                System.err.println(f.getAbsolutePath());
                System.err.print(y);
                System.err.println("-------------------------\n");
            }
            
        }

        System.out.println("Test time: " +  ((System.currentTimeMillis() - time) / 1000) + " sec");
        System.out.println(((((float)accuracy)/testSet.size())*100)+"% accuracy");
        showChart(trainingCosts);
    }

    private static boolean compare(GMatrix expected, GMatrix actual) {
        if (expected.getNumCol() != actual.getNumCol())
            return false;
        if (expected.getNumRow() != actual.getNumRow())
            return false;
        for(int i=0;i<actual.getNumRow();i++)
            for(int j=0;j<actual.getNumCol();j++)
                if (Math.round(expected.getElement(i, j)) != Math.round(actual.getElement(i, j)))
                    return false;
        return true;
    }
    
    public static String size(GMatrix mtx) {
        return mtx.getNumRow() + "x" + mtx.getNumCol();
    }
}
