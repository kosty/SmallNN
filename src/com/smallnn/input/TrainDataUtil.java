package com.smallnn.input;

import static com.smallnn.input.FileUtil.find;
import static com.smallnn.input.FileUtil.pngFilter;
import static com.smallnn.input.FileUtil.tildeExpand;
import static com.smallnn.input.ImageUtil.HEIGHT;
import static com.smallnn.input.ImageUtil.IMAGE_SIZE;
import static com.smallnn.input.ImageUtil.SMALL_HEIGHT;
import static com.smallnn.input.ImageUtil.SMALL_WIDTH;
import static com.smallnn.input.ImageUtil.WIDTH;
import static com.smallnn.input.ImageUtil.imageToGrayscale;
import static com.smallnn.input.ImageUtil.resize;

import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

import javax.imageio.ImageIO;
import javax.vecmath.GMatrix;

public class TrainDataUtil {
    
    public static final String dataPrefix = "~/References/barsDetect";
    public static final String TRAINING_DIR = "training";
    public static final String VALIDATION_DIR = "verify";
    public static final String TEST_DIR = "test";
    
    public static File getDataDir(String set, String cluster, String resolution){
        return tildeExpand(dataPrefix+"/"+set+"/"+cluster+"/"+resolution);
    }

    public static File getDataDir(String set){
        return tildeExpand(dataPrefix+"/"+set);
    }
    
    public static double[] longToDoubleArray(long[] arr){
        double[] result = new double[arr.length];
        for(int i=0;i<arr.length;i++){
            result[i] = (double)arr[i];
        }
        return result;
    }
    
    public static Data getMixedData() throws Exception{
        List<File> trainingFileSet = interleaveEntries();
        int inputLayerSize = SMALL_WIDTH*SMALL_HEIGHT;
        int classes = 2;
        int m=trainingFileSet.size();
        double[] x = new double[m*inputLayerSize];
        double[] y = new double[m*classes];
        
        for(int i=0;i<m;i++){
            File trainFile = trainingFileSet.get(i%trainingFileSet.size());
            double[] tmparr = longToDoubleArray(imageToGrayscale(resize(ImageIO.read(trainFile), SMALL_WIDTH, SMALL_HEIGHT))); 
            for (int j=0;j<inputLayerSize;j++){
                x[i*inputLayerSize+j] = tmparr[j];
            }
            y[i*classes]   = i%2 == 0 ? 1. : 0.;
            y[i*classes+1] = i%2 == 0 ? 0. : 1.;
        }
        
        return new Data(new GMatrix(m, inputLayerSize, x), new GMatrix(m, classes, y));
    }
    
    public static GMatrix readData(File f) throws Exception{
        double[] data = longToDoubleArray(imageToGrayscale(resize(ImageIO.read(f), SMALL_WIDTH, SMALL_HEIGHT)));
        return new GMatrix(1, SMALL_WIDTH*SMALL_HEIGHT, data);
    }
    
    
    public static List<File> interleaveEntries(){
        List<File> bars = find(getDataDir(TRAINING_DIR, "bars", "72x60"), pngFilter);
        List<File> nonBars = find(getDataDir(TRAINING_DIR, "non-bars", "72x60"), pngFilter);
        List<File> interleaved = new ArrayList<File>();
        for(int i=0;i<bars.size();i++){
            interleaved.add(bars.get(i));
            interleaved.add(nonBars.get(i%nonBars.size()));
        }
        return interleaved;
    }
    
    public static List<File> listSetFiles(String set) throws Exception {
        return find(getDataDir(set), pngFilter);
    }
    
    public static class Data {
        
        public final GMatrix x;
        public final GMatrix y;

        public Data(GMatrix x, GMatrix y) {
            this.x = x;
            this.y = y;
        }
    }
    
    
    public static Data randomInit() {
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
    
    public static Data blackWhiteInit(){
        double[] white = new double[]{16777215., 16777215., 16777215., 16777215., 16777215., 16777215., 16777215., 16777215.}; 
        double[] black = new double[]{0., 0., 0., 0., 0., 0., 0., 0.}; 
        int m = 50;
        int inputLayerSize = white.length;
        int classes =2;
        double[] x = new double[m * inputLayerSize];
        double[] y = new double[m * classes];
        for(int i=0;i<m;i++){
            double[] tmparr = i%2 == 0 ? white : black; 
            for (int j=0;j<inputLayerSize;j++){
                x[i*inputLayerSize+j] = tmparr[j];
            }
            y[i*classes]   = i%2 == 0 ? 1. : 0.;
            y[i*classes+1] = i%2 == 0 ? 0. : 1.;
        }
        
        return new Data(new GMatrix(m, inputLayerSize, x), new GMatrix(m, classes, y));
    }
    
    public static Data controlledImageSetInit() throws Exception{
        List<File> bars = find(getDataDir(TRAINING_DIR, "bars", "72x60"), pngFilter);
        List<File> nonBars = find(getDataDir(TRAINING_DIR, "non-bars", "72x60"), pngFilter);
        double[] bars0 = longToDoubleArray(imageToGrayscale(resize(ImageIO.read(bars.get(0)), SMALL_HEIGHT, SMALL_WIDTH))); //bars
        double[] nonBars1 =longToDoubleArray(imageToGrayscale(resize(ImageIO.read(nonBars.get(0)), SMALL_HEIGHT, SMALL_WIDTH))); //non-bars
        
        double[] bars2 = longToDoubleArray(imageToGrayscale(resize(ImageIO.read(bars.get(1)), SMALL_HEIGHT, SMALL_WIDTH))); //bars
        double[] nonBars3 =longToDoubleArray(imageToGrayscale(resize(ImageIO.read(nonBars.get(1)), SMALL_HEIGHT, SMALL_WIDTH))); //non-bars
        
        double[] bars4 = longToDoubleArray(imageToGrayscale(resize(ImageIO.read(bars.get(2)), SMALL_HEIGHT, SMALL_WIDTH))); //bars
        double[] nonBars5 =longToDoubleArray(imageToGrayscale(resize(ImageIO.read(nonBars.get(2)), SMALL_HEIGHT, SMALL_WIDTH))); //non-bars
        
        double[] bars6 = longToDoubleArray(imageToGrayscale(resize(ImageIO.read(bars.get(3)), SMALL_HEIGHT, SMALL_WIDTH))); //bars
        double[] nonBars7 =longToDoubleArray(imageToGrayscale(resize(ImageIO.read(nonBars.get(3)), SMALL_HEIGHT, SMALL_WIDTH))); //non-bars
        
        double[] bars8 = longToDoubleArray(imageToGrayscale(resize(ImageIO.read(bars.get(4)), SMALL_HEIGHT, SMALL_WIDTH))); //bars
        double[] nonBars9 =longToDoubleArray(imageToGrayscale(resize(ImageIO.read(nonBars.get(4)), SMALL_HEIGHT, SMALL_WIDTH))); //non-bars
        
        double[][] imagesArray = new double[][]{bars0, nonBars1, bars2, nonBars3, bars4, nonBars5, bars6, nonBars7, bars8, nonBars9};
        int inputLayerSize = bars0.length;
        int m=imagesArray.length*40;
        int classes = 2;
        double[] x = new double[m*inputLayerSize];
        double[] y = new double[m*classes];
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
        return new Data(new GMatrix(m, inputLayerSize, x), new GMatrix(m, classes, y));
    }

    public static Data fullImageInit() throws Exception {
        List<File> filenames = interleaveEntries();
        int classes = 2;
        double[] x = new double[filenames.size() * IMAGE_SIZE];
        double[] y = new double[filenames.size() * classes];
        for (int l = 0; l < filenames.size(); l++) {
            StringBuilder sb = new StringBuilder();
            int cnt = 0;

            BufferedImage img = ImageIO.read(filenames.get(l));
            sb.append(filenames.get(l).getName());
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
                if (filenames.get(l).getAbsolutePath().toLowerCase().contains("non-bars")) {
                    y[i] = 0.;
                    y[i + 1] = 1.;
                } else {
                    y[i] = 1.;
                    y[i + 1] = 0.;
                }
            }
        }
        return new Data(new GMatrix(filenames.size(), IMAGE_SIZE, x), new GMatrix(filenames.size(), classes, y));
    }

    
}
