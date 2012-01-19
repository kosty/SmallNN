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
package com.smallnn.input;

import static com.smallnn.input.FileUtil.dirFilter;
import static com.smallnn.input.FileUtil.find;
import static com.smallnn.input.FileUtil.pngFilter;
import static com.smallnn.input.FileUtil.tildeExpand;
import static com.smallnn.input.ImageUtil.imageToDoubleArray;
import static com.smallnn.input.ImageUtil.zealousSubsCrop;

import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.imageio.ImageIO;
import javax.vecmath.GMatrix;

import com.smallnn.input.ImageUtil.Resolution;

public class TrainDataUtil {

    public static final String TRAINING_DIR = "training";
    public static final String VALIDATION_DIR = "verify";
    public static final String TEST_DIR = "test";

    public static List<File> getSubDirs(String path) {
        File[] listFiles = tildeExpand(path).listFiles(dirFilter);
        return Arrays.asList(listFiles);
    }
    
    public static Data getSingleClassCroppedData(String kitPath, Resolution orig, int boostLevel) throws Exception {
        List<File> trainingFileSet = interleaveTrainingClasses(kitPath, TRAINING_DIR, orig.width + "x" + orig.height);
        int inputLayerSize = (orig.width/2) * (orig.height/2);
        int m = trainingFileSet.size()*boostLevel;
        
        double[] x = new double[m * inputLayerSize];
        double[] y = new double[m];
        
        for (int i = 0; i < m; i++) {
            File f = trainingFileSet.get(i % trainingFileSet.size());
            double[] image = imageToDoubleArray(zealousSubsCrop(ImageIO.read(f)));
            System.arraycopy(image, 0, x, i * inputLayerSize, inputLayerSize);

            y[i] = i%2 == 0 ? 1. : 0.;
        }
        
        return new Data(new GMatrix(m, inputLayerSize, x), new GMatrix(m, 1, y));
    }
    
    public static Data getSingleClassData(String kitPath, Resolution orig, int boostLevel) throws Exception {
        List<File> trainingFileSet = interleaveTrainingClasses(kitPath, TRAINING_DIR, orig.width + "x" + orig.height);
        int inputLayerSize = orig.width * orig.height;
        int m = trainingFileSet.size()*boostLevel;
        
        double[] x = new double[m * inputLayerSize];
        double[] y = new double[m];
        
        for (int i = 0; i < m; i++) {
            File f = trainingFileSet.get(i % trainingFileSet.size());
            double[] image = imageToDoubleArray(ImageIO.read(f));
            System.arraycopy(image, 0, x, i * inputLayerSize, inputLayerSize);

            y[i] = i%2 == 0 ? 1. : 0.;
        }
        
        return new Data(new GMatrix(m, inputLayerSize, x), new GMatrix(m, 1, y));
    }

    public static Data getData(String kitPath, Resolution orig, Resolution fin, int boostLevel) throws Exception {
        int classes = getClassesNumber(kitPath);
        List<File> trainingFileSet = interleaveTrainingClasses(kitPath, TRAINING_DIR, orig.width + "x" + orig.height);
        int inputLayerSize = fin.width * fin.height;

        int m = trainingFileSet.size() * boostLevel;
        double[] x = new double[m * inputLayerSize];
        double[] y = new double[m * classes];

        for (int i = 0; i < m; i++) {
            File f = trainingFileSet.get(i % trainingFileSet.size());
            double[] image = imageToDoubleArray(ImageIO.read(f));
            System.arraycopy(image, 0, x, i * inputLayerSize, inputLayerSize);

            double[] expectedClass = new double[classes];
            Arrays.fill(expectedClass, 0.);
            expectedClass[i % classes] = 1.;
            System.arraycopy(expectedClass, 0, y, i * classes, classes);
        }
        return new Data(new GMatrix(m, inputLayerSize, x), new GMatrix(m, classes, y));
    }

    public static TestData getTestData(String kitPath, Resolution orig, Resolution fin)
            throws Exception {
        List<File> files = interleaveTrainingClasses(kitPath, TEST_DIR, orig.width + "x" + orig.height);
        GMatrix[] x_norms = new GMatrix[files.size()];
        GMatrix[] ys = new GMatrix[files.size()];
        int classes = TrainDataUtil.getClassesNumber(kitPath);
        for (int i = 0; i < files.size(); i++) {
            BufferedImage img = ImageIO.read(files.get(i));

            double[] image = imageToDoubleArray(img);
            x_norms[i] = new GMatrix(1, fin.width * fin.height, image);

            double[] y_test = new double[classes];
            Arrays.fill(y_test, 0.);
            y_test[i % classes] = 1.;
            ys[i] = new GMatrix(1, classes, y_test);
        }
        return new TestData(x_norms, ys, files);
    }

    public static int getClassesNumber(String pathPrefix) {
        List<File> testClasses = getSubDirs(pathPrefix + "/" + TEST_DIR);
        List<File> trainingClasses = getSubDirs(pathPrefix + "/" + TRAINING_DIR);
        assert testClasses.size() == trainingClasses.size();
        return trainingClasses.size();
    }

    public static List<File> interleaveTrainingClasses(String kitPath, String set, String resolution) {
        List<File> trainingClasses = getSubDirs(kitPath + "/" + set);
        List<File> refClass = find(new File(trainingClasses.get(0), resolution), pngFilter);
        System.out.println("Reference class of training set: " + trainingClasses.get(0).getAbsolutePath());

        List<List<File>> otherClasses = new ArrayList<List<File>>();
        for (int i = 1; i < trainingClasses.size(); i++) {
            otherClasses.add(find(new File(trainingClasses.get(i), resolution), pngFilter));
        }
        List<File> interleaved = new ArrayList<File>();
        for (int i = 0; i < refClass.size(); i++) {
            interleaved.add(refClass.get(i));
            for (List<File> otherClass : otherClasses)
                interleaved.add(otherClass.get(i % otherClass.size()));
        }
        return interleaved;
    }

    public static class Data {

        public final GMatrix x;
        public final GMatrix y;

        public Data(GMatrix x, GMatrix y) {
            this.x = x;
            this.y = y;
        }
    }

    public static class TestData {
        public final GMatrix[] x;
        public final GMatrix[] y;
        public final List<File> files;

        public TestData(GMatrix[] x, GMatrix[] y, List<File> files) {
            this.x = x;
            this.y = y;
            this.files = files;
        }
    }

    public static Data randomData() {
        int m = 50;
        int inputLayerSize = 2;
        int classes = 2;
        double[] x = new double[m * inputLayerSize];
        double[] y = new double[m * classes];
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
        return new Data(new GMatrix(m, inputLayerSize, x), new GMatrix(m, classes, y));
    }

    public static Data blackWhiteData() {
        double[] white = new double[] { 16777215., 16777215., 16777215., 16777215., 16777215., 16777215., 16777215.,
                16777215. };
        double[] black = new double[] { 0., 0., 0., 0., 0., 0., 0., 0. };
        int m = 50;
        int inputLayerSize = white.length;
        int classes = 2;
        double[] x = new double[m * inputLayerSize];
        double[] y = new double[m * classes];
        for (int i = 0; i < m; i++) {
            double[] tmparr = i % 2 == 0 ? white : black;
            for (int j = 0; j < inputLayerSize; j++) {
                x[i * inputLayerSize + j] = tmparr[j];
            }
            y[i * classes] = i % 2 == 0 ? 1. : 0.;
            y[i * classes + 1] = i % 2 == 0 ? 0. : 1.;
        }

        return new Data(new GMatrix(m, inputLayerSize, x), new GMatrix(m, classes, y));
    }

    public static Data interlieveTrainData(List<BufferedImage> subsFrames, List<BufferedImage> nonsubsFrames)
            throws IOException {
        double[] subsY = new double[] { 1., 0. };
        double[] nonSubsY = new double[] { 0., 1. };
        int size = subsFrames.size();
        int imageSize = subsFrames.get(0).getWidth() * subsFrames.get(0).getHeight();
        int totalEntries = 4 * size;
        double[] x = new double[totalEntries * imageSize];
        double[] y = new double[totalEntries * 2];
        for (int i = 0; i < size; i++) {
            double[] subsImage = imageToDoubleArray(subsFrames.get(i));
            double[] nonsubsImage = imageToDoubleArray(nonsubsFrames.get(i));

            int index = i * 4;

            System.arraycopy(subsImage, 0, x, index * imageSize, imageSize);
            System.arraycopy(subsY, 0, y, index * 2, 2);

            System.arraycopy(subsImage, 0, x, (index+1) * imageSize, imageSize);
            System.arraycopy(subsY, 0, y, (index+1) * 2, 2);
            
            System.arraycopy(subsImage, 0, x, (index+2) * imageSize, imageSize);
            System.arraycopy(subsY, 0, y, (index+2) * 2, 2);
            
            System.arraycopy(nonsubsImage, 0, x, (index + 3) * imageSize, imageSize);
            System.arraycopy(nonSubsY, 0, y, (index + 3) * 2, 2);
        }
        return new Data(new GMatrix(totalEntries, imageSize, x), new GMatrix(totalEntries, 2, y));
    }

    public static GMatrix readData(BufferedImage img) {
        int imageSize = img.getWidth() * img.getHeight();
        double[] x = imageToDoubleArray(img);
        return new GMatrix(1, imageSize, x);
    }

}
